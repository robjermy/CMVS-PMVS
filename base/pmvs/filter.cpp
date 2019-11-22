#include <numeric>
#include <ctime>
#include <time.h>
#include "../numeric/mylapack.hpp"
#include "findMatch.hpp"
#include "filter.hpp"

PMVS3::CFilter::CFilter(CFindMatch& findMatch) : _fm(findMatch) {}

void PMVS3::CFilter::init(void) {}

void PMVS3::CFilter::run(void) {
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

void PMVS3::CFilter::filterOutside(void) {
  time_t tv;
  time(&tv);
  time_t curtime = tv;
  std::cerr << "FilterOutside" << std::endl;
  //??? notice (1) here to avoid removing _fix=1
  _fm._pos.collectPatches(1);

  const int psize = (int)_fm._pos._ppatches.size();
  _gains.resize(psize);

  std::cerr << "mainbody: " << std::flush;

  _fm.Count() = 0;
  std::vector<std::thread> threads(_fm._CPU);
  for (int i = 0; i < _fm._CPU; ++i) {
    threads[i] = std::thread(&CFilter::filterOutsideThread, this);
  }

  for (int i = 0; i < _fm._CPU; ++i) {
    if (threads[i].joinable()) {
      threads[i].join();
    }
  }
  std::cerr << std::endl;

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

  if (denom == 0) {
    denom = 1;
  }
  ave /= denom;
  ave2 /= denom;
  ave2 = sqrt(std::max(0.0, ave2 - ave * ave));
  std::cerr << "Gain (ave/var): " << ave << ' ' << ave2 << std::endl;

  time(&tv);
  std::cerr << (int)_fm._pos._ppatches.size() << " -> "
       << (int)_fm._pos._ppatches.size() - count << " ("
       << 100 * ((int)_fm._pos._ppatches.size() - count) / (float)_fm._pos._ppatches.size()
       << "%)\t" << (tv - curtime) / CLOCKS_PER_SEC << " secs" << std::endl;
}

float PMVS3::CFilter::computeGain(const Patch::CPatch& patch, const int lock) {
  float gain = patch.score2(_fm._nccThreshold);

  const int size = (int)patch._images.size();
  for (int i = 0; i < size; ++i) {
    const int& index = patch._images[i];
    if (_fm.NumTargetImages() <= index) continue;

    const int& ix = patch._grids[i][0];
    const int& iy = patch._grids[i][1];
    const int index2 = iy * _fm._pos._gwidths[index] + ix;

    float maxpressure = 0.0f;
    if (lock) {
      _fm.LockSharedImage(index);
    }

    for (int j = 0; j < (int)_fm._pos._pgrids[index][index2].size(); ++j) {
      if (!_fm.isNeighbor(patch, *_fm._pos._pgrids[index][index2][j], _fm._neighborThreshold1)) {
        maxpressure = std::max(maxpressure, _fm._pos._pgrids[index][index2][j]->_ncc - _fm._nccThreshold);
      }
    }
    if (lock) {
      _fm.UnlockSharedImage(index);
    }

    gain -= maxpressure;
  }

  const int vsize = (int)patch._vimages.size();
  for (int i = 0; i < vsize; ++i) {
    const int& index = patch._vimages[i];
    if (_fm.NumTargetImages() <= index) continue;

    const float pdepth = _fm._pss.computeDepth(index, patch._coord);

    const int& ix = patch._vgrids[i][0];
    const int& iy = patch._vgrids[i][1];
    const int index2 = iy * _fm._pos._gwidths[index] + ix;
    float maxpressure = 0.0f;

    if (lock) {
      _fm.LockSharedImage(index);
    }

    for (int j = 0; j < (int)_fm._pos._pgrids[index][index2].size(); ++j) {
      const float bdepth = _fm._pss.computeDepth(index, _fm._pos._pgrids[index][index2][j]->_coord);
      if (pdepth < bdepth && !_fm.isNeighbor(patch, *_fm._pos._pgrids[index][index2][j], _fm._neighborThreshold1)) {
        maxpressure = std::max(maxpressure, _fm._pos._pgrids[index][index2][j]->_ncc - _fm._nccThreshold);
      }
    }
    if (lock) {
      _fm.UnlockSharedImage(index);
    }

    gain -= maxpressure;
  }
  return gain;
}

void PMVS3::CFilter::filterOutsideThread(void) {
  _fm.Lock();
  const int id = _fm.Count()++;
  _fm.Unlock();

  const int size = (int)_fm._pos._ppatches.size();
  const int itmp = (int)ceil(size / (float)_fm._CPU);
  const int begin = id * itmp;
  const int end = std::min(size, (id + 1) * itmp);

  for (int p = begin; p < end; ++p) {
    Patch::PPatch& ppatch = _fm._pos._ppatches[p];
    _gains[p] = ppatch->score2(_fm._nccThreshold);

    const int size = (int)ppatch->_images.size();
    for (int i = 0; i < size; ++i) {
      const int& index = ppatch->_images[i];
      if (_fm.NumTargetImages() <= index) continue;

      const int& ix = ppatch->_grids[i][0];
      const int& iy = ppatch->_grids[i][1];
      const int index2 = iy * _fm._pos._gwidths[index] + ix;

      float maxpressure = 0.0f;
      for (int j = 0; j < (int)_fm._pos._pgrids[index][index2].size(); ++j) {
	      if (!_fm.isNeighbor(*ppatch, *_fm._pos._pgrids[index][index2][j], _fm._neighborThreshold1))
	        maxpressure = std::max(maxpressure, _fm._pos._pgrids[index][index2][j]->_ncc - _fm._nccThreshold);
      }

      _gains[p] -= maxpressure;
    }

    const int vsize = (int)ppatch->_vimages.size();
    for (int i = 0; i < vsize; ++i) {
      const int& index = ppatch->_vimages[i];
      if (_fm.NumTargetImages() <= index) continue;

      const float pdepth = _fm._pss.computeDepth(index, ppatch->_coord);

      const int& ix = ppatch->_vgrids[i][0];
      const int& iy = ppatch->_vgrids[i][1];
      const int index2 = iy * _fm._pos._gwidths[index] + ix;
      float maxpressure = 0.0f;

      for (int j = 0; j < (int)_fm._pos._pgrids[index][index2].size(); ++j) {
        const float bdepth = _fm._pss.computeDepth(index, _fm._pos._pgrids[index][index2][j]->_coord);
        if (pdepth < bdepth && !_fm.isNeighbor(*ppatch, *_fm._pos._pgrids[index][index2][j], _fm._neighborThreshold1)) {
	        maxpressure = std::max(maxpressure, _fm._pos._pgrids[index][index2][j]->_ncc - _fm._nccThreshold);
        }
      }
      _gains[p] -= maxpressure;
    }
  }
}

void PMVS3::CFilter::filterExact(void) {
  time_t tv;
  time(&tv);
  time_t curtime = tv;
  std::cerr << "Filter Exact: " << std::flush;

  //??? cannot use (1) because we use patch._id to set newimages,....
  _fm._pos.collectPatches();
  const int psize = (int)_fm._pos._ppatches.size();

  // dis associate images
  _newimages.clear();
  _newgrids.clear();
  _removeimages.clear();
  _removegrids.clear();

  _newimages.resize(psize);
  _newgrids.resize(psize);
  _removeimages.resize(psize);
  _removegrids.resize(psize);

  _fm.Count() = 0;
  std::vector<std::thread> threads0(_fm._CPU);
  for (int i = 0; i < _fm._CPU; ++i) {
    threads0[i] = std::thread(&CFilter::filterExactThread, this);
  }

  for (int i = 0; i < _fm._CPU; ++i) {
    if (threads0[i].joinable()) {
      threads0[i].join();
    }
  }
  std::cerr << std::endl;

  //----------------------------------------------------------------------
  for (int p = 0; p < psize; ++p) {
    if (_fm._pos._ppatches[p]->_fix) continue;

    for (int i = 0; i < (int)_removeimages[p].size(); ++i) {
      const int index = _removeimages[p][i];
      if (_fm.NumTargetImages() <= index) {
	      std::cerr << "MUST NOT COME HERE" << std::endl;
        exit (1);
      }
      const int ix = _removegrids[p][i][0];
      const int iy = _removegrids[p][i][1];
      const int index2 = iy * _fm._pos._gwidths[index] + ix;

      _fm._pos._pgrids[index][index2].erase(remove(_fm._pos._pgrids[index][index2].begin(), _fm._pos._pgrids[index][index2].end(), _fm._pos._ppatches[p]), _fm._pos._pgrids[index][index2].end());
    }
  }

  _fm._debug = 1;

  int count = 0;
  for (int p = 0; p < psize; ++p) {
    if (_fm._pos._ppatches[p]->_fix) continue;

    Patch::CPatch& patch = *_fm._pos._ppatches[p];

    // This should be images in targetting images. Has to come before the next for-loop.
    patch._timages = (int)_newimages[p].size();

    for (int i = 0; i < (int)patch._images.size(); ++i) {
      const int& index = patch._images[i];
      if (_fm.NumTargetImages() <= index) {
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
  std::cerr << (int)_fm._pos._ppatches.size() << " -> "
       << (int)_fm._pos._ppatches.size() - count << " ("
       << 100 * ((int)_fm._pos._ppatches.size() - count) / (float)_fm._pos._ppatches.size()
       << "%)\t" << (tv - curtime) / CLOCKS_PER_SEC << " secs" << std::endl;
}

void PMVS3::CFilter::filterExactThread(void) {
  const int psize = (int)_fm._pos._ppatches.size();
  std::vector<std::vector<int> > newimages, removeimages;
  std::vector<std::vector<TVec2<int> > > newgrids, removegrids;
  newimages.resize(psize);
  removeimages.resize(psize);
  newgrids.resize(psize);
  removegrids.resize(psize);

  while (1) {
    _fm.Lock();
    const int image = _fm.Count()++;
    _fm.Unlock();

    if (_fm.NumTargetImages() <= image) break;

    std::cerr << '*' << std::flush;

    const int& w = _fm._pos._gwidths[image];
    const int& h = _fm._pos._gheights[image];
    int index = -1;
    for (int y = 0; y < h; ++y) {
      for (int x = 0; x < w; ++x) {
        ++index;
        for (int i = 0; i < (int)_fm._pos._pgrids[image][index].size(); ++i) {
          const Patch::CPatch& patch = *_fm._pos._pgrids[image][index][i];
          if (patch._fix) continue;

          int safe = 0;

          if (_fm._pos.isVisible(patch, image, x, y, _fm._neighborThreshold1, 0)) {
            safe = 1;
          } else if (0 < x && _fm._pos.isVisible(patch, image, x - 1, y, _fm._neighborThreshold1, 0)) { // use 4 neighbors?
            safe = 1;
          } else if (x < w - 1 && _fm._pos.isVisible(patch, image, x + 1, y, _fm._neighborThreshold1, 0)) {
            safe = 1;
          } else if (0 < y && _fm._pos.isVisible(patch, image, x, y - 1, _fm._neighborThreshold1, 0)) {
            safe = 1;
          } else if (y < h - 1 && _fm._pos.isVisible(patch, image, x, y + 1, _fm._neighborThreshold1, 0)) {
            safe = 1;
          }
          if (safe) {
            newimages[patch._id].push_back(image);
            newgrids[patch._id].push_back(TVec2<int>(x, y));
          } else {
            removeimages[patch._id].push_back(image);
            removegrids[patch._id].push_back(TVec2<int>(x, y));
          }
        }
      }
    }
  }

  _fm.Lock();
  for (int p = 0; p < psize; ++p) {
    _newimages[p].insert(_newimages[p].end(), newimages[p].begin(), newimages[p].end());
    _newgrids[p].insert(_newgrids[p].end(), newgrids[p].begin(), newgrids[p].end());
    _removeimages[p].insert(_removeimages[p].end(), removeimages[p].begin(), removeimages[p].end());
    _removegrids[p].insert(_removegrids[p].end(), removegrids[p].begin(), removegrids[p].end());
  }
  _fm.Unlock();
}

void PMVS3::CFilter::filterNeighborThread(void) {
  const int size = (int)_fm._pos._ppatches.size();
  while (1) {
    int jtmp = -1;
    _fm.Lock();
    if (!_fm.Jobs().empty()) {
      jtmp = _fm.Jobs().front();
      _fm.Jobs().pop_front();
    }
    _fm.Unlock();
    if (jtmp == -1) break;

    const int begin = _fm.JUnit() * jtmp;
    const int end = std::min(size, _fm.JUnit() * (jtmp + 1));

    for (int p = begin; p < end; ++p) {
      Patch::PPatch& ppatch = _fm._pos._ppatches[p];
      if (_rejects[p]) continue;

      std::vector<Patch::PPatch> neighbors;
      //_fm._pos.findNeighbors(*ppatch, neighbors, 0, 4, 2);
      _fm._pos.findNeighbors(*ppatch, neighbors, 0, 4, 2, 1);

      //?? new filter
      if ((int)neighbors.size() < 6) {
        //if ((int)neighbors.size() < 8)
        _rejects[p] = _time + 1;
      } else {
        // Fit a quadratic surface
        if (filterQuad(*ppatch, neighbors)) {
          _rejects[p] = _time + 1;
        }
      }
    }
  }
}

int PMVS3::CFilter::filterQuad(const Patch::CPatch& patch, const std::vector<Patch::PPatch>& neighbors) const {
  std::vector<std::vector<float> > A;
  std::vector<float> b, x;

  Vec4f xdir, ydir;
  ortho(patch._normal, xdir, ydir);

  const int nsize = (int)neighbors.size();

  float h = 0.0f;
  for (int n = 0; n < nsize; ++n) {
    h += norm(neighbors[n]->_coord - patch._coord);
  }
  h /= nsize;

  A.resize(nsize);
  b.resize(nsize);

  std::vector<float> fxs, fys, fzs;
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
  const int inum = std::min(_fm._tau, (int)patch._images.size());
  float unit = 0.0;
  //for (int i = 0; i < (int)patch._images.size(); ++i)
  for (int i = 0; i < inum; ++i) {
    unit += _fm._optim.getUnit(patch._images[i], patch._coord);
  }
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

  if (residual < _fm._quadThreshold) {
    return 0;
  } else {
    return 1;
  }
}

void PMVS3::CFilter::filterNeighbor(const int times) {
  time_t tv;
  time(&tv);
  time_t curtime = tv;
  std::cerr << "FilterNeighbor:\t" << std::flush;

  //??? notice (1) to avoid removing _fix=1
  _fm._pos.collectPatches(1);
  if (_fm._pos._ppatches.empty()) return;

  _rejects.resize((int)_fm._pos._ppatches.size());
  fill(_rejects.begin(), _rejects.end(), 0);

  // Lapack is not thread-safe? Sometimes, the code gets stuck here.
  int count = 0;
  for (_time = 0; _time < times; ++_time) {
    _fm.Count() = 0;

    _fm.Jobs().clear();
    const int jtmp = (int)ceil(_fm._pos._ppatches.size() / (float)_fm.JUnit());
    for (int j = 0; j < jtmp; ++j) {
      _fm.Jobs().push_back(j);
    }

    std::vector<std::thread> threads(_fm._CPU);
    for (int i = 0; i < _fm._CPU; ++i) {
      threads[i] = std::thread(&CFilter::filterNeighborThread, this);
    }

    for (int i = 0; i < _fm._CPU; ++i) {
      if (threads[i].joinable()) {
        threads[i].join();
      }
    }

    auto bpatch = _fm._pos._ppatches.begin();
    auto epatch = _fm._pos._ppatches.end();
    auto breject = _rejects.begin();

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
  std::cerr << (int)_fm._pos._ppatches.size() << " -> "
       << (int)_fm._pos._ppatches.size() - count << " ("
       << 100 * ((int)_fm._pos._ppatches.size() - count) / (float)_fm._pos._ppatches.size()
       << "%)\t" << (tv - curtime) / CLOCKS_PER_SEC << " secs" << std::endl;
}

//----------------------------------------------------------------------
// Take out small connected components
//----------------------------------------------------------------------
void PMVS3::CFilter::filterSmallGroups(void) {
  time_t tv;
  time(&tv);
  time_t curtime = tv;
  std::cerr << "FilterGroups:\t" << std::flush;
  _fm._pos.collectPatches();
  if (_fm._pos._ppatches.empty()) return;

  const int psize = (int)_fm._pos._ppatches.size();
  std::vector<int> label;
  label.resize(psize);
  std::fill(label.begin(), label.end(), -1);

  std::list<int> untouch;
  auto bpatch = _fm._pos._ppatches.begin();
  for (int p = 0; p < psize; ++p, ++bpatch) {
    untouch.push_back(p);
    (*bpatch)->_flag = p;
  }

  int id = -1;
  while (!untouch.empty()) {
    const int pid = untouch.front();
    untouch.pop_front();

    if (label[pid] != -1) continue;

    label[pid] = ++id;
    std::list<int> ltmp;
    ltmp.push_back(pid);

    while (!ltmp.empty()) {
      const int ptmp = ltmp.front();
      ltmp.pop_front();

      filterSmallGroupsSub(ptmp, id, label, ltmp);
    }
  }
  id++;

  std::vector<int> size;
  size.resize(id);
  auto bite = label.begin();
  auto eite = label.end();
  while (bite != eite) {
    ++size[*bite];
    ++bite;
  }

  const int threshold = std::max(20, psize / 10000);
  std::cerr << threshold << std::endl;

  bite = size.begin();
  eite = size.end();
  while (bite != eite) {
    if (*bite < threshold) {
      *bite = 0;
    } else {
      *bite = 1;
    }
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
  std::cerr << (int)_fm._pos._ppatches.size() << " -> "
       << (int)_fm._pos._ppatches.size() - count << " ("
       << 100 * ((int)_fm._pos._ppatches.size() - count) / (float)_fm._pos._ppatches.size()
       << "%)\t" << (tv - curtime)/CLOCKS_PER_SEC << " secs" << std::endl;
}

void PMVS3::CFilter::filterSmallGroupsSub(const int pid, const int id, std::vector<int>& label, std::list<int>& ltmp) const {
  // find neighbors of ptmp and set their ids
  const Patch::CPatch& patch = *_fm._pos._ppatches[pid];

  const int index = patch._images[0];
  const int ix = patch._grids[0][0];
  const int iy = patch._grids[0][1];
  const int gwidth = _fm._pos._gwidths[index];
  const int gheight = _fm._pos._gheights[index];

  for (int y = -1; y <= 1; ++y) {
    const int iytmp = iy + y;
    if (iytmp < 0 || gheight <= iytmp) continue;

    for (int x = -1; x <= 1; ++x) {
      const int ixtmp = ix + x;
      if (ixtmp < 0 || gwidth <= ixtmp) continue;

      const int index2 = iytmp * gwidth + ixtmp;
      auto bgrid = _fm._pos._pgrids[index][index2].begin();
      auto egrid = _fm._pos._pgrids[index][index2].end();
      while (bgrid != egrid) {
        const int itmp = (*bgrid)->_flag;
        if (label[itmp] != -1) {
          ++bgrid;
          continue;
        }

        if (_fm.isNeighbor(patch, **bgrid, _fm._neighborThreshold2)) {
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

        if (_fm.isNeighbor(patch, **bgrid, _fm._neighborThreshold2)) {
          label[itmp] = id;
          ltmp.push_back(itmp);
        }
        ++bgrid;
      }
    }
  }
}

void PMVS3::CFilter::setDepthMaps(void) {
  // initialize
  for (int index = 0; index < _fm.NumTargetImages(); ++index) {
    std::fill(_fm._pos._dpgrids[index].begin(), _fm._pos._dpgrids[index].end(), _fm._pos._MAXDEPTH);
  }

  _fm.Count() = 0;
  std::vector<std::thread> threads(_fm._CPU);
  for (int i = 0; i < _fm._CPU; ++i) {
    threads[i] = std::thread(&CFilter::setDepthMapsThread, this);
  }

  for (int i = 0; i < _fm._CPU; ++i) {
    if (threads[i].joinable()) {
      threads[i].join();
    }
  }
}


void PMVS3::CFilter::setDepthMapsThread(void) {
  while (1) {
    _fm.Lock();
    const int index = _fm.Count()++;
    _fm.Unlock();

    if (_fm.NumTargetImages() <= index) break;

    const int gwidth = _fm._pos._gwidths[index];
    const int gheight = _fm._pos._gheights[index];

    auto bpatch = _fm._pos._ppatches.begin();
    auto epatch = _fm._pos._ppatches.end();

    while (bpatch != epatch) {
      Patch::PPatch& ppatch = *bpatch;
      const Vec3f icoord = _fm._pss.project(index, ppatch->_coord, _fm._level);

      const float fx = icoord[0] / _fm._csize;
      const int xs[2] = {(int)floor(fx), (int)ceil(fx)};
      const float fy = icoord[1] / _fm._csize;
      const int ys[2] = {(int)floor(fy), (int)ceil(fy)};

      const float depth = _fm._pss._photos[index].OpticalAxis() * ppatch->_coord;

      for (int j = 0; j < 2; ++j) {
        for (int i = 0; i < 2; ++i) {
          if (xs[i] < 0 || gwidth <= xs[i] || ys[j] < 0 || gheight <= ys[j]) continue;

          const int index2 = ys[j] * gwidth + xs[i];

          if (_fm._pos._dpgrids[index][index2] == _fm._pos._MAXDEPTH) {
            _fm._pos._dpgrids[index][index2] = ppatch;
          } else {
            const float dtmp = _fm._pss._photos[index].OpticalAxis() * _fm._pos._dpgrids[index][index2]->_coord;

            if (depth < dtmp) {
              _fm._pos._dpgrids[index][index2] = ppatch;
            }
          }
        }
      }
      ++bpatch;
    }
  }
}

void PMVS3::CFilter::setDepthMapsVGridsVPGridsAddPatchV(const int additive) {
  _fm._pos.collectPatches();
  setDepthMaps();

  // clear _vpgrids
  for (int index = 0; index < _fm.NumTargetImages(); ++index) {
    auto bvvp = _fm._pos._vpgrids[index].begin();
    auto evvp = _fm._pos._vpgrids[index].end();
    while (bvvp != evvp) {
      (*bvvp).clear();
      ++bvvp;
    }
  }

  if (additive == 0) {
    // initialization
    auto bpatch = _fm._pos._ppatches.begin();
    auto epatch = _fm._pos._ppatches.end();
    while (bpatch != epatch) {
      (*bpatch)->_vimages.clear();
      (*bpatch)->_vgrids.clear();
      ++bpatch;
    }
  }

  _fm.Count() = 0;
  // std::vector<thrd_t> threads0(_fm._CPU);
  std::vector<std::thread> threads0(_fm._CPU);
  for (int i = 0; i < _fm._CPU; ++i) {
    threads0[i] = std::thread(&CFilter::setVGridsVPGridsThread, this);
  }

  for (int i = 0; i < _fm._CPU; ++i) {
    if (threads0[i].joinable()) {
      threads0[i].join();
    }
  }

  _fm.Count() = 0;
  std::vector<std::thread> threads1(_fm._CPU);
  for (int i = 0; i < _fm._CPU; ++i) {
    threads1[i] = std::thread(&CFilter::addPatchVThread, this);
  }

  for (int i = 0; i < _fm._CPU; ++i) {
    if (threads1[i].joinable()) {
      threads1[i].join();
    }
  }
}

void PMVS3::CFilter::setVGridsVPGridsThread(void) {
  const int noj = 1000;
  const int size = (int)_fm._pos._ppatches.size();
  const int job = std::max(1, size / (noj - 1));

  while (1) {
    _fm.Lock();
    const int id = _fm.Count()++;
    _fm.Unlock();

    const int begin = id * job;
    const int end = std::min(size, (id + 1) * job);

    if (size <= begin) break;

    // add patches to _vpgrids
    for (int p = begin; p < end; ++p) {
      Patch::PPatch& ppatch = _fm._pos._ppatches[p];
      _fm._pos.setVImagesVGrids(ppatch);
    }
  }
}

void PMVS3::CFilter::addPatchVThread(void) {
  while (1) {
    _fm.Lock();
    const int index = _fm.Count()++;
    _fm.Unlock();

    if (_fm.NumTargetImages() <= index) break;

    auto bpatch = _fm._pos._ppatches.begin();
    auto epatch = _fm._pos._ppatches.end();
    while (bpatch != epatch) {
      Patch::PPatch& ppatch = *bpatch;
      auto bimage = ppatch->_vimages.begin();
      auto eimage = ppatch->_vimages.end();
      auto bgrid = ppatch->_vgrids.begin();

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
