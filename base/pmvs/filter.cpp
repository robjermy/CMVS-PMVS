#include "tinycthread.h"
#include <numeric>
#include <ctime>
#include <time.h>
#include "../numeric/mylapack.hpp"
#include "findMatch.hpp"
#include "filter.hpp"

using namespace Patch;
using namespace PMVS3;
using namespace std;

CFilter::CFilter(CFindMatch& findMatch) : _fm(findMatch) {
}

void CFilter::init(void) {
}

void CFilter::run(void) {
  setDepthMapsVGridsVPGridsAddPatchV(0);

  filterOutside();
  setDepthMapsVGridsVPGridsAddPatchV(1);

  filterExact();
  setDepthMapsVGridsVPGridsAddPatchV(1);
  
  filterNeighbor(1);
  setDepthMapsVGridsVPGridsAddPatchV(1);
  
  filterSmallGroups();
  setDepthMapsVGridsVPGridsAddPatchV(1);
}

void CFilter::filterOutside(void) {
  time_t tv;
  time(&tv); 
  time_t curtime = tv;
  cerr << "FilterOutside" << endl;
  //??? notice (1) here to avoid removing _fix=1
  _fm._pos.collectPatches(1);

  const int psize = (int)_fm._pos._ppatches.size();  
  _gains.resize(psize);
  
  cerr << "mainbody: " << flush;
  
  _fm._count = 0;
  vector<thrd_t> threads(_fm._CPU);
  for (int i = 0; i < _fm._CPU; ++i)
    thrd_create(&threads[i], &filterOutsideThreadTmp, (void*)this);
  for (int i = 0; i < _fm._CPU; ++i)
    thrd_join(threads[i], NULL);
  cerr << endl;

  // delete patches with positive _gains
  int count = 0;

  double ave = 0.0f;
  double ave2 = 0.0f;
  int denom = 0;  
  
  for (int p = 0; p < psize; ++p) {
    ave += _gains[p];
    ave2 += _gains[p] * _gains[p];
    ++denom;
    
    if (_gains[p] < 0.0) {
      _fm._pos.removePatch(_fm._pos._ppatches[p]);
      count++;
    }
  }

  if (denom == 0)
    denom = 1;
  ave /= denom;
  ave2 /= denom;
  ave2 = sqrt(max(0.0, ave2 - ave * ave));
  cerr << "Gain (ave/var): " << ave << ' ' << ave2 << endl;
  
  time(&tv);
  cerr << (int)_fm._pos._ppatches.size() << " -> "
       << (int)_fm._pos._ppatches.size() - count << " ("
       << 100 * ((int)_fm._pos._ppatches.size() - count) / (float)_fm._pos._ppatches.size()
       << "%)\t" << (tv - curtime) / CLOCKS_PER_SEC << " secs" << endl;
}

float CFilter::computeGain(const Patch::CPatch& patch, const int lock) {
  float gain = patch.score2(_fm._nccThreshold);

  const int size = (int)patch._images.size();  
  for (int i = 0; i < size; ++i) {
    const int& index = patch._images[i];
    if (_fm._tnum <= index)
      continue;
    
    const int& ix = patch._grids[i][0];      const int& iy = patch._grids[i][1];
    const int index2 = iy * _fm._pos._gwidths[index] + ix;
    
    float maxpressure = 0.0f;
    if (lock)
      _fm._imageLocks[index].rdlock();
    
    for (int j = 0; j < (int)_fm._pos._pgrids[index][index2].size(); ++j) {
      if (!_fm.isNeighbor(patch, *_fm._pos._pgrids[index][index2][j],
                           _fm._neighborThreshold1))
        maxpressure = max(maxpressure, _fm._pos._pgrids[index][index2][j]->_ncc -
                          _fm._nccThreshold);
    }
    if (lock)
      _fm._imageLocks[index].unlock();
    
    gain -= maxpressure;
  }
  
  const int vsize = (int)patch._vimages.size();
  for (int i = 0; i < vsize; ++i) {
    const int& index = patch._vimages[i];
    if (_fm._tnum <= index)
      continue;
    
    const float pdepth = _fm._pss.computeDepth(index, patch._coord);    
    
    const int& ix = patch._vgrids[i][0];      const int& iy = patch._vgrids[i][1];
    const int index2 = iy * _fm._pos._gwidths[index] + ix;
    float maxpressure = 0.0f;      

    if (lock)
      _fm._imageLocks[index].rdlock();
    
    for (int j = 0; j < (int)_fm._pos._pgrids[index][index2].size(); ++j) {
      const float bdepth = _fm._pss.computeDepth(index, _fm._pos._pgrids[index][index2][j]->_coord);
      if (pdepth < bdepth &&
          !_fm.isNeighbor(patch, *_fm._pos._pgrids[index][index2][j],
                           _fm._neighborThreshold1)) {
        maxpressure = max(maxpressure,
                          _fm._pos._pgrids[index][index2][j]->_ncc -
                          _fm._nccThreshold);
      }
    }
    if (lock)
      _fm._imageLocks[index].unlock();
    
    gain -= maxpressure;
  }
  return gain;
}

void CFilter::filterOutsideThread(void) {
  mtx_lock(&_fm._lock);
  const int id = _fm._count++;
  mtx_unlock(&_fm._lock);

  const int size = (int)_fm._pos._ppatches.size();  
  const int itmp = (int)ceil(size / (float)_fm._CPU);
  const int begin = id * itmp;
  const int end = min(size, (id + 1) * itmp);
  
  for (int p = begin; p < end; ++p) {
    PPatch& ppatch = _fm._pos._ppatches[p];
    _gains[p] = ppatch->score2(_fm._nccThreshold);
    
    const int size = (int)ppatch->_images.size();  
    for (int i = 0; i < size; ++i) {
      const int& index = ppatch->_images[i];
      if (_fm._tnum <= index)
        continue;
      
      const int& ix = ppatch->_grids[i][0];      const int& iy = ppatch->_grids[i][1];
      const int index2 = iy * _fm._pos._gwidths[index] + ix;
      
      float maxpressure = 0.0f;
      for (int j = 0; j < (int)_fm._pos._pgrids[index][index2].size(); ++j) {
	if (!_fm.isNeighbor(*ppatch, *_fm._pos._pgrids[index][index2][j],
                             _fm._neighborThreshold1))
	  maxpressure = max(maxpressure, _fm._pos._pgrids[index][index2][j]->_ncc -
			    _fm._nccThreshold);
      }
      
      _gains[p] -= maxpressure;
    }

    const int vsize = (int)ppatch->_vimages.size();
    for (int i = 0; i < vsize; ++i) {
      const int& index = ppatch->_vimages[i];
      if (_fm._tnum <= index)
        continue;
      
      const float pdepth = _fm._pss.computeDepth(index, ppatch->_coord);    
      
      const int& ix = ppatch->_vgrids[i][0];      const int& iy = ppatch->_vgrids[i][1];
      const int index2 = iy * _fm._pos._gwidths[index] + ix;
      float maxpressure = 0.0f;      
      
      for (int j = 0; j < (int)_fm._pos._pgrids[index][index2].size(); ++j) {
        const float bdepth = _fm._pss.computeDepth(index, _fm._pos._pgrids[index][index2][j]->_coord);
	if (pdepth < bdepth &&
            !_fm.isNeighbor(*ppatch, *_fm._pos._pgrids[index][index2][j],
                             _fm._neighborThreshold1)) {
	  maxpressure = max(maxpressure,
                            _fm._pos._pgrids[index][index2][j]->_ncc -
			    _fm._nccThreshold);
        }
      }
      _gains[p] -= maxpressure;
    }
  }
}

int CFilter::filterOutsideThreadTmp(void* arg) {
  ((CFilter*)arg)->filterOutsideThread();
  return 0;
}

void CFilter::filterExact(void) {
  time_t tv;
  time(&tv); 
  time_t curtime = tv;
  cerr << "Filter Exact: " << flush;

  //??? cannot use (1) because we use patch._id to set newimages,....
  _fm._pos.collectPatches();
  const int psize = (int)_fm._pos._ppatches.size();
  
  // dis associate images
  _newimages.clear();     _newgrids.clear();
  _removeimages.clear();  _removegrids.clear();
  _newimages.resize(psize);     _newgrids.resize(psize);
  _removeimages.resize(psize);  _removegrids.resize(psize);

  _fm._count = 0;
  vector<thrd_t> threads0(_fm._CPU);
  for (int i = 0; i < _fm._CPU; ++i)
    thrd_create(&threads0[i], &filterExactThreadTmp, (void*)this);  
  for (int i = 0; i < _fm._CPU; ++i)
    thrd_join(threads0[i], NULL);
  cerr << endl;

  //----------------------------------------------------------------------
  for (int p = 0; p < psize; ++p) {
    if (_fm._pos._ppatches[p]->_fix)
      continue;
    
    for (int i = 0; i < (int)_removeimages[p].size(); ++i) {
      const int index = _removeimages[p][i];
      if (_fm._tnum <= index) {
	cerr << "MUST NOT COME HERE" << endl;        exit (1);
      }
      const int ix = _removegrids[p][i][0];      const int iy = _removegrids[p][i][1];
      const int index2 = iy * _fm._pos._gwidths[index] + ix;

      _fm._pos._pgrids[index][index2].
	erase(remove(_fm._pos._pgrids[index][index2].begin(),
		     _fm._pos._pgrids[index][index2].end(),
		     _fm._pos._ppatches[p]),
	      _fm._pos._pgrids[index][index2].end());
    }
  }
  
  _fm._debug = 1;
  
  int count = 0;
  for (int p = 0; p < psize; ++p) {
    if (_fm._pos._ppatches[p]->_fix)
      continue;
    
    CPatch& patch = *_fm._pos._ppatches[p];

    // This should be images in targetting images. Has to come before the next for-loop.
    patch._timages = (int)_newimages[p].size();
    
    for (int i = 0; i < (int)patch._images.size(); ++i) {
      const int& index = patch._images[i];
      if (_fm._tnum <= index) {
	_newimages[p].push_back(patch._images[i]);
	_newgrids[p].push_back(patch._grids[i]);
      }
    }

    patch._images.swap(_newimages[p]);
    patch._grids.swap(_newgrids[p]);

    if (_fm._minImageNumThreshold <= (int)patch._images.size()) {
      _fm._optim.setRefImage(patch, 0);
      _fm._pos.setGrids(patch);
    }

    if ((int)patch._images.size() < _fm._minImageNumThreshold) {
      _fm._pos.removePatch(_fm._pos._ppatches[p]);
      count++;
    }
  }
  time(&tv); 
  cerr << (int)_fm._pos._ppatches.size() << " -> "
       << (int)_fm._pos._ppatches.size() - count << " ("
       << 100 * ((int)_fm._pos._ppatches.size() - count) / (float)_fm._pos._ppatches.size()
       << "%)\t" << (tv - curtime) / CLOCKS_PER_SEC << " secs" << endl;
}

void CFilter::filterExactThread(void) {
  const int psize = (int)_fm._pos._ppatches.size();
  vector<vector<int> > newimages, removeimages;
  vector<vector<TVec2<int> > > newgrids, removegrids;
  newimages.resize(psize);  removeimages.resize(psize);
  newgrids.resize(psize);   removegrids.resize(psize);

  while (1) {
    mtx_lock(&_fm._lock);
    const int image = _fm._count++;
    mtx_unlock(&_fm._lock);

    if (_fm._tnum <= image)
      break;
    
    cerr << '*' << flush;
    
    const int& w = _fm._pos._gwidths[image];
    const int& h = _fm._pos._gheights[image];
    int index = -1;
    for (int y = 0; y < h; ++y) {
      for (int x = 0; x < w; ++x) {
        ++index;
	for (int i = 0; i < (int)_fm._pos._pgrids[image][index].size(); ++i) {
	  const CPatch& patch = *_fm._pos._pgrids[image][index][i];
          if (patch._fix)
            continue;
          
	  int safe = 0;

	  if (_fm._pos.isVisible(patch, image, x, y, _fm._neighborThreshold1, 0))
	    safe = 1;
	  // use 4 neighbors?
	  else if (0 < x && _fm._pos.isVisible(patch, image, x - 1, y, _fm._neighborThreshold1, 0))
	    safe = 1;
	  else if (x < w - 1 && _fm._pos.isVisible(patch, image, x + 1, y, _fm._neighborThreshold1, 0))
	    safe = 1;
	  else if (0 < y && _fm._pos.isVisible(patch, image, x, y - 1, _fm._neighborThreshold1, 0))
	    safe = 1;
	  else if (y < h - 1 && _fm._pos.isVisible(patch, image, x, y + 1, _fm._neighborThreshold1, 0))
	    safe = 1;
	
	  if (safe) {
	    newimages[patch._id].push_back(image);
	    newgrids[patch._id].push_back(TVec2<int>(x, y));
	  }
	  else {
	    removeimages[patch._id].push_back(image);
	    removegrids[patch._id].push_back(TVec2<int>(x, y));
	  }
	}
      }
    }
  }

  mtx_lock(&_fm._lock);
  for (int p = 0; p < psize; ++p) {
    _newimages[p].insert(_newimages[p].end(),
			  newimages[p].begin(), newimages[p].end());
    _newgrids[p].insert(_newgrids[p].end(),
			  newgrids[p].begin(), newgrids[p].end());
    _removeimages[p].insert(_removeimages[p].end(),
			  removeimages[p].begin(), removeimages[p].end());
    _removegrids[p].insert(_removegrids[p].end(),
			  removegrids[p].begin(), removegrids[p].end());
  }
  mtx_unlock(&_fm._lock);
}

int CFilter::filterExactThreadTmp(void* arg) {
  ((CFilter*)arg)->filterExactThread();
  return 0;
}

void CFilter::filterNeighborThread(void) {
  const int size = (int)_fm._pos._ppatches.size();  
  while (1) {
    int jtmp = -1;
    mtx_lock(&_fm._lock);
    if (!_fm._jobs.empty()) {
      jtmp = _fm._jobs.front();
      _fm._jobs.pop_front();
    }
    mtx_unlock(&_fm._lock);
    if (jtmp == -1)
      break;

    const int begin = _fm._junit * jtmp;
    const int end = min(size, _fm._junit * (jtmp + 1));

    for (int p = begin; p < end; ++p) {
      PPatch& ppatch = _fm._pos._ppatches[p];
      if (_rejects[p])
        continue;
      
      vector<PPatch> neighbors;
      //_fm._pos.findNeighbors(*ppatch, neighbors, 0, 4, 2);
      _fm._pos.findNeighbors(*ppatch, neighbors, 0, 4, 2, 1);
      
      //?? new filter
      if ((int)neighbors.size() < 6)
        //if ((int)neighbors.size() < 8)
        _rejects[p] = _time + 1;
      else {
        // Fit a quadratic surface
        if (filterQuad(*ppatch, neighbors))
          _rejects[p] = _time + 1;
      }
    }
  }

  /*
  mtx_lock(&_fm._lock);
  const int id = _fm._count++;
  mtx_unlock(&_fm._lock);

  const int size = (int)_fm._pos._ppatches.size();  
  const int itmp = (int)ceil(size / (float)_fm._CPU);
  const int begin = id * itmp;
  const int end = min(size, (id + 1) * itmp);

  for (int p = begin; p < end; ++p) {
    PPatch& ppatch = _fm._pos._ppatches[p];
    if (_rejects[p])
      continue;

    vector<PPatch> neighbors;
    _fm._pos.findNeighbors(*ppatch, neighbors, 0, 4, 2);

    //?? new filter
    if ((int)neighbors.size() < 6)
    //if ((int)neighbors.size() < 8)
      _rejects[p] = _time + 1;
    else {
      // Fit a quadratic surface
      if (filterQuad(*ppatch, neighbors))
        _rejects[p] = _time + 1;
    }
  }
  */
}

int CFilter::filterQuad(const Patch::CPatch& patch,
                        const std::vector<PPatch>& neighbors) const {
  vector<vector<float> > A;
  vector<float> b, x;

  Vec4f xdir, ydir;
  ortho(patch._normal, xdir, ydir);

  const int nsize = (int)neighbors.size();

  float h = 0.0f;
  for (int n = 0; n < nsize; ++n)
    h += norm(neighbors[n]->_coord - patch._coord);
  h /= nsize;
  
  A.resize(nsize);
  b.resize(nsize);

  vector<float> fxs, fys, fzs;
  fxs.resize(nsize);
  fys.resize(nsize);
  fzs.resize(nsize);
  for (int n = 0; n < nsize; ++n) {
    A[n].resize(5);
    Vec4f diff = neighbors[n]->_coord - patch._coord;
    fxs[n] = diff * xdir / h;
    fys[n] = diff * ydir / h;
    fzs[n] = diff * patch._normal;

    A[n][0] = fxs[n] * fxs[n];
    A[n][1] = fys[n] * fys[n];
    A[n][2] = fxs[n] * fys[n];
    A[n][3] = fxs[n];
    A[n][4] = fys[n];
    b[n] = fzs[n];
  }
  x.resize(5);
  Cmylapack::lls(A, b, x);
  
  // Compute residual divided by _dscale
  const int inum = min(_fm._tau, (int)patch._images.size());
  float unit = 0.0;
  //for (int i = 0; i < (int)patch._images.size(); ++i)
  for (int i = 0; i < inum; ++i)
    unit += _fm._optim.getUnit(patch._images[i], patch._coord);
  //unit /= (int)patch._images.size();
  unit /= inum;
  
  float residual = 0.0f;
  for (int n = 0; n < nsize; ++n) {
    const float res =
      x[0] * (fxs[n] * fxs[n]) +
      x[1] * (fys[n] * fys[n]) +
      x[2] * (fxs[n] * fys[n]) +
      x[3] * fxs[n] +
      x[4] * fys[n] - fzs[n];
    //residual += fabs(res) / neighbors[n]->_dscale;
    residual += fabs(res) / unit;
  }
  
  residual /= (nsize - 5);

  if (residual < _fm._quadThreshold)
    return 0;
  else
    return 1;
}

int CFilter::filterNeighborThreadTmp(void* arg) {
  ((CFilter*)arg)->filterNeighborThread();
  return 0;
}
  
void CFilter::filterNeighbor(const int times) {
  time_t tv;
  time(&tv); 
  time_t curtime = tv;
  cerr << "FilterNeighbor:\t" << flush;

  //??? notice (1) to avoid removing _fix=1
  _fm._pos.collectPatches(1);
  if (_fm._pos._ppatches.empty())
    return;
  
  _rejects.resize((int)_fm._pos._ppatches.size());
  fill(_rejects.begin(), _rejects.end(), 0);

  // Lapack is not thread-safe? Sometimes, the code gets stuck here.
  int count = 0;
  for (_time = 0; _time < times; ++_time) {
    _fm._count = 0;

    _fm._jobs.clear();
    const int jtmp = (int)ceil(_fm._pos._ppatches.size() /
                               (float)_fm._junit);
    for (int j = 0; j < jtmp; ++j)
      _fm._jobs.push_back(j);
    
    vector<thrd_t> threads(_fm._CPU);
    for (int i = 0; i < _fm._CPU; ++i)
      thrd_create(&threads[i], &filterNeighborThreadTmp, (void*)this);
    for (int i = 0; i < _fm._CPU; ++i)
      thrd_join(threads[i], NULL);
    
    vector<PPatch>::iterator bpatch = _fm._pos._ppatches.begin();
    vector<PPatch>::iterator epatch = _fm._pos._ppatches.end();
    vector<int>::iterator breject = _rejects.begin();
    
    while (bpatch != epatch) {
      if ((*breject) == _time + 1) {
        count++;
        _fm._pos.removePatch(*bpatch);;
      }

      ++bpatch;
      ++breject;
    }
  }
  time(&tv);
  cerr << (int)_fm._pos._ppatches.size() << " -> "
       << (int)_fm._pos._ppatches.size() - count << " ("
       << 100 * ((int)_fm._pos._ppatches.size() - count) / (float)_fm._pos._ppatches.size()
       << "%)\t" << (tv - curtime) / CLOCKS_PER_SEC << " secs" << endl;
}

//----------------------------------------------------------------------
// Take out small connected components
//----------------------------------------------------------------------
void CFilter::filterSmallGroups(void) {
  time_t tv;
  time(&tv); 
  time_t curtime = tv;
  cerr << "FilterGroups:\t" << flush;
  _fm._pos.collectPatches();
  if (_fm._pos._ppatches.empty())
    return;
  
  const int psize = (int)_fm._pos._ppatches.size();
  vector<int> label;
  label.resize(psize);
  fill(label.begin(), label.end(), -1);

  list<int> untouch;
  vector<PPatch>::iterator bpatch = _fm._pos._ppatches.begin();
  for (int p = 0; p < psize; ++p, ++bpatch) {
    untouch.push_back(p);
    (*bpatch)->_flag = p;
  }
  
  int id = -1;
  while (!untouch.empty()) {
    const int pid = untouch.front();
    untouch.pop_front();

    if (label[pid] != -1)
      continue;

    label[pid] = ++id;
    list<int> ltmp;
    ltmp.push_back(pid);

    while (!ltmp.empty()) {
      const int ptmp = ltmp.front();
      ltmp.pop_front();

      filterSmallGroupsSub(ptmp, id, label, ltmp);
    }
  }  
  id++;
  
  vector<int> size;
  size.resize(id);
  vector<int>::iterator bite = label.begin();
  vector<int>::iterator eite = label.end();
  while (bite != eite) {
    ++size[*bite];
    ++bite;
  }
  
  const int threshold = max(20, psize / 10000);
  cerr << threshold << endl;
  
  bite = size.begin();
  eite = size.end();
  while (bite != eite) {
    if (*bite < threshold)
      *bite = 0;
    else
      *bite = 1;
    ++bite;
  }

  int count = 0;

  bite = label.begin();
  eite = label.end();
  bpatch = _fm._pos._ppatches.begin();
  while (bite != eite) {
    if ((*bpatch)->_fix) {
      ++bite;
      ++bpatch;
      continue;
    }
    
    if (size[*bite] == 0) {
      _fm._pos.removePatch(*bpatch);
      count++;
    }
    ++bite;
    ++bpatch;
  }
  time(&tv);
  cerr << (int)_fm._pos._ppatches.size() << " -> "
       << (int)_fm._pos._ppatches.size() - count << " ("
       << 100 * ((int)_fm._pos._ppatches.size() - count) / (float)_fm._pos._ppatches.size()
       << "%)\t" << (tv - curtime)/CLOCKS_PER_SEC << " secs" << endl;
}

void CFilter::filterSmallGroupsSub(const int pid, const int id,
                                   std::vector<int>& label,
                                   std::list<int>& ltmp) const {  
  // find neighbors of ptmp and set their ids
  const CPatch& patch = *_fm._pos._ppatches[pid];
  
  const int index = patch._images[0];
  const int ix = patch._grids[0][0];
  const int iy = patch._grids[0][1];
  const int gwidth = _fm._pos._gwidths[index];
  const int gheight = _fm._pos._gheights[index];
  
  for (int y = -1; y <= 1; ++y) {
    const int iytmp = iy + y;
    if (iytmp < 0 || gheight <= iytmp)
      continue;
    for (int x = -1; x <= 1; ++x) {
      const int ixtmp = ix + x;
      if (ixtmp < 0 || gwidth <= ixtmp)
	continue;
      
      //if (1 < abs(x) + abs(y))
      //continue;

      const int index2 = iytmp * gwidth + ixtmp;
      vector<PPatch>::iterator bgrid = _fm._pos._pgrids[index][index2].begin();
      vector<PPatch>::iterator egrid = _fm._pos._pgrids[index][index2].end();
      while (bgrid != egrid) {
        const int itmp = (*bgrid)->_flag;
	if (label[itmp] != -1) {
          ++bgrid;
	  continue;
        }
	
	if (_fm.isNeighbor(patch, **bgrid,
			    _fm._neighborThreshold2)) {
	  label[itmp] = id;
	  ltmp.push_back(itmp);
	}
        ++bgrid;
      }
      bgrid = _fm._pos._vpgrids[index][index2].begin();
      egrid = _fm._pos._vpgrids[index][index2].end();
      while (bgrid != egrid) {
        const int itmp = (*bgrid)->_flag;
	if (label[itmp] != -1) {
          ++bgrid;
	  continue;
        }
	
	if (_fm.isNeighbor(patch, **bgrid,
			    _fm._neighborThreshold2)) {
	  label[itmp] = id;
	  ltmp.push_back(itmp);
	}
        ++bgrid;
      }
    }
  }
}

void CFilter::setDepthMaps(void) {
  // initialize
  for (int index = 0; index < _fm._tnum; ++index) {
    fill(_fm._pos._dpgrids[index].begin(), _fm._pos._dpgrids[index].end(),
         _fm._pos._MAXDEPTH);
  }
  
  _fm._count = 0;
  vector<thrd_t> threads(_fm._CPU);
  for (int i = 0; i < _fm._CPU; ++i)
    thrd_create(&threads[i], &setDepthMapsThreadTmp, (void*)this);
  for (int i = 0; i < _fm._CPU; ++i)
    thrd_join(threads[i], NULL);
}

int CFilter::setDepthMapsThreadTmp(void* arg) {
  ((CFilter*)arg)->setDepthMapsThread();
  return 0;
}

void CFilter::setDepthMapsThread(void) {
  while (1) {
    mtx_lock(&_fm._lock);
    const int index = _fm._count++;
    mtx_unlock(&_fm._lock);

    if (_fm._tnum <= index)
      break;

    const int gwidth = _fm._pos._gwidths[index];
    const int gheight = _fm._pos._gheights[index];

    vector<PPatch>::iterator bpatch = _fm._pos._ppatches.begin();
    vector<PPatch>::iterator epatch = _fm._pos._ppatches.end();

    while (bpatch != epatch) {
      PPatch& ppatch = *bpatch;
      const Vec3f icoord =
        _fm._pss.project(index, ppatch->_coord, _fm._level);
      
      const float fx = icoord[0] / _fm._csize;
      const int xs[2] = {(int)floor(fx), (int)ceil(fx)};
      const float fy = icoord[1] / _fm._csize;
      const int ys[2] = {(int)floor(fy), (int)ceil(fy)};
    
      const float depth =
        _fm._pss._photos[index].OpticalAxis() * ppatch->_coord;

      for (int j = 0; j < 2; ++j) {
	for (int i = 0; i < 2; ++i) {
	  if (xs[i] < 0 || gwidth <= xs[i] || ys[j] < 0 || gheight <= ys[j])
	    continue;
          const int index2 = ys[j] * gwidth + xs[i];
	
	  if (_fm._pos._dpgrids[index][index2] == _fm._pos._MAXDEPTH)
	    _fm._pos._dpgrids[index][index2] = ppatch;
	  else {
	    const float dtmp = _fm._pss._photos[index].OpticalAxis() *
	      _fm._pos._dpgrids[index][index2]->_coord;
	    
	    if (depth < dtmp)
	      _fm._pos._dpgrids[index][index2] = ppatch;
	  }
	}
      }
      ++bpatch;
    }
  }
}  

void CFilter::setDepthMapsVGridsVPGridsAddPatchV(const int additive) {
  _fm._pos.collectPatches();
  setDepthMaps();

  // clear _vpgrids
  for (int index = 0; index < _fm._tnum; ++index) {
    vector<vector<PPatch> >::iterator bvvp = _fm._pos._vpgrids[index].begin();
    vector<vector<PPatch> >::iterator evvp = _fm._pos._vpgrids[index].end();
    while (bvvp != evvp) {
      (*bvvp).clear();
      ++bvvp;
    }
  }

  if (additive == 0) {
    // initialization
    vector<PPatch>::iterator bpatch = _fm._pos._ppatches.begin();
    vector<PPatch>::iterator epatch = _fm._pos._ppatches.end();
    while (bpatch != epatch) {
      (*bpatch)->_vimages.clear();
      (*bpatch)->_vgrids.clear();
      ++bpatch;
    }
  }
    
  _fm._count = 0;
  vector<thrd_t> threads0(_fm._CPU);
  for (int i = 0; i < _fm._CPU; ++i)
    thrd_create(&threads0[i], &setVGridsVPGridsThreadTmp, (void*)this);  
  for (int i = 0; i < _fm._CPU; ++i)
    thrd_join(threads0[i], NULL);
  
  _fm._count = 0;
  vector<thrd_t> threads1(_fm._CPU);
  for (int i = 0; i < _fm._CPU; ++i)
    thrd_create(&threads1[i], &addPatchVThreadTmp, (void*)this);
  for (int i = 0; i < _fm._CPU; ++i)
    thrd_join(threads1[i], NULL);
}

int CFilter::setVGridsVPGridsThreadTmp(void* arg) {
  ((CFilter*)arg)->setVGridsVPGridsThread();
  return 0;
}

int CFilter::addPatchVThreadTmp(void* arg) {
  ((CFilter*)arg)->addPatchVThread();
  return 0;
}

void CFilter::setVGridsVPGridsThread(void) {
  const int noj = 1000;
  const int size = (int)_fm._pos._ppatches.size();    
  const int job = max(1, size / (noj - 1));
  
  while (1) {
    mtx_lock(&_fm._lock);
    const int id = _fm._count++;
    mtx_unlock(&_fm._lock);

    const int begin = id * job;
    const int end = min(size, (id + 1) * job);
    
    if (size <= begin)
      break;
    
    // add patches to _vpgrids
    for (int p = begin; p < end; ++p) {
      PPatch& ppatch = _fm._pos._ppatches[p];
      _fm._pos.setVImagesVGrids(ppatch);
    }
  }
}

void CFilter::addPatchVThread(void) {
  while (1) {
    mtx_lock(&_fm._lock);
    const int index = _fm._count++;
    mtx_unlock(&_fm._lock);
    
    if (_fm._tnum <= index)
      break;
    
    vector<PPatch>::iterator bpatch = _fm._pos._ppatches.begin();
    vector<PPatch>::iterator epatch = _fm._pos._ppatches.end();
    while (bpatch != epatch) {
      PPatch& ppatch = *bpatch;
      vector<int>::iterator bimage = ppatch->_vimages.begin();
      vector<int>::iterator eimage = ppatch->_vimages.end();
      vector<Vec2i>::iterator bgrid = ppatch->_vgrids.begin();
      while (bimage != eimage) {
	if (*bimage == index) {
	  const int& ix = (*bgrid)[0];
	  const int& iy = (*bgrid)[1];
          const int index2 = iy * _fm._pos._gwidths[index] + ix;
          _fm._pos._vpgrids[index][index2].push_back(ppatch);
	  break;
	}
        ++bimage;
        ++bgrid;
      }
      ++bpatch;
    }
  }
}
