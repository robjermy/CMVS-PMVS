#define _USE_MATH_DEFINES
#include <cmath>

#include <iomanip>
#include <limits>
#include <string>
#include "patchOrganizerS.hpp"
#include "findMatch.hpp"

Patch::PPatch PMVS3::CPatchOrganizerS::_MAXDEPTH(new Patch::CPatch());
Patch::PPatch PMVS3::CPatchOrganizerS::_BACKGROUND(new Patch::CPatch());

PMVS3::CPatchOrganizerS::CPatchOrganizerS(CFindMatch& findMatch) : _fm(findMatch) {}

// change the contents of _images from images to indexes
void PMVS3::CPatchOrganizerS::image2index(Patch::CPatch& patch) {
  // first image has to be target image
  std::vector<int> newimages;
  for (int i = 0; i < (int)patch._images.size(); ++i) {
    const int index = _fm._pss.image2index(patch._images[i]);
    if (index != -1) {
      newimages.push_back(index);
    }
  }

  patch._images.swap(newimages);

  // make sure that the reference image is the targeting image
  int exist = -1;
  for (int j = 0; j < (int)patch._images.size(); ++j) {
    if (patch._images[j] < _fm._tnum) {
      exist = j;
      break;
    }
  }
  if (exist == -1) {
    patch._images.clear();
  } else if (exist != 0) {
    std::swap(patch._images[0], patch._images[exist]);
  }
}

// change the contents of _images from indexes to images
void PMVS3::CPatchOrganizerS::index2image(Patch::CPatch& patch) {
  for (int i = 0; i < (int)patch._images.size(); ++i) {
    patch._images[i] = _fm._pss._images[patch._images[i]];
  }

  for (int i = 0; i < (int)patch._vimages.size(); ++i) {
    patch._vimages[i] = _fm._pss._images[patch._vimages[i]];
  }
}

void PMVS3::CPatchOrganizerS::init(void) {
  _pgrids.clear();
  _pgrids.resize(_fm._tnum);

  _vpgrids.clear();
  _vpgrids.resize(_fm._tnum);

  _dpgrids.clear();
  _dpgrids.resize(_fm._tnum);

  _counts.clear();
  _counts.resize(_fm._tnum);

  _gwidths.clear();
  _gwidths.resize(_fm._num);

  _gheights.clear();
  _gheights.resize(_fm._num);

  for (int index = 0; index < _fm._num; ++index) {
    const int gwidth = (_fm._pss.getWidth(index, _fm._level) + _fm._csize - 1) / _fm._csize;
    const int gheight = (_fm._pss.getHeight(index, _fm._level) + _fm._csize - 1) / _fm._csize;
    _gwidths[index] = gwidth;
    _gheights[index] = gheight;

    if (index < _fm._tnum) {
      _pgrids[index].resize(gwidth * gheight);
      _vpgrids[index].resize(gwidth * gheight);
      _dpgrids[index].resize(gwidth * gheight);
      _counts[index].resize(gwidth * gheight);
      std::fill(_dpgrids[index].begin(), _dpgrids[index].end(), _MAXDEPTH);
    }
  }
}

void PMVS3::CPatchOrganizerS::writePatches2(const std::string prefix, bool bExportPLY, bool bExportPatch, bool bExportPSet) {
  collectPatches(1);

  if (bExportPLY) {
    char buffer[1024];
    sprintf(buffer, "%s.ply", prefix.c_str());
    writePLY(_ppatches, buffer);
  }

  if (bExportPatch) {
    char buffer[1024];
    sprintf(buffer, "%s.patch", prefix.c_str());
    std::ofstream ofstr;
    ofstr.open(buffer);

    // Set output to full precision to be able to fully round-trip
    // the text representation with what we have in memory.
    ofstr << std::setprecision(std::numeric_limits<double>::max_digits10);

    ofstr << "PATCHES" << std::endl
          << (int)_ppatches.size() << std::endl;
    for (int p = 0; p < (int)_ppatches.size(); ++p) {
      Patch::CPatch patch = *_ppatches[p];
      index2image(patch);
      ofstr << patch << "\n";
    }
    ofstr.close();
  }

  if (bExportPSet) {
    char buffer[1024];
    sprintf(buffer, "%s.pset", prefix.c_str());
    std::ofstream ofstr;
    ofstr.open(buffer);
    for (int p = 0; p < (int)_ppatches.size(); ++p)
      ofstr << _ppatches[p]->_coord[0] << ' '
            << _ppatches[p]->_coord[1] << ' '
            << _ppatches[p]->_coord[2] << ' '
            << _ppatches[p]->_normal[0] << ' '
            << _ppatches[p]->_normal[1] << ' '
            << _ppatches[p]->_normal[2] << "\n";
    ofstr.close();
  }
}

void PMVS3::CPatchOrganizerS::readPatches(void) {
  // Read-in existing reconstructed points. set _fix to one for non-targeting images
  for (int i = 0; i < _fm._tnum; ++i) {
    const int image = _fm._images[i];
    char buffer[1024];
    sprintf(buffer, "%smodels/%08d.patc%d", _fm._prefix.c_str(), image, _fm._level);
    std::ifstream ifstr;
    ifstr.open(buffer);
    if (!ifstr.is_open()) continue;

    int pnum;
    std::string header;

    ifstr >> header >> pnum;
    std::cerr << image << ' ' << pnum << " patches" << std::endl;

    for (int p = 0; p < pnum; ++p) {
      Patch::PPatch ppatch(new Patch::CPatch());
      ifstr >> *ppatch;
      ppatch->_fix = 0;
      ppatch->_vimages.clear();

      image2index(*ppatch);
      if (ppatch->_images.empty()) continue;

      // _vimages must be targeting images
#ifdef DEBUG
      for (int j = 0; j < (int)ppatch->_vimages.size(); ++j) {
        if (_fm._tnum <= ppatch->_vimages[j]) {
          cerr << "Impossible in readPatches. _vimages must be targeting images" << std::endl
               << "for patches stored in targeting images, if visdata2 have been consistent" << std::endl;
          exit (1);
        }
      }
#endif
      setGrids(*ppatch);
      addPatch(ppatch);
    }

    ifstr.close();
  }

  // For patches in non-targeting images
  for (int i = _fm._tnum; i < _fm._num; ++i) {
    const int image = _fm._images[i];
    char buffer[1024];
    sprintf(buffer, "%smodels/%08d.patc%d", _fm._prefix.c_str(), image, _fm._level);
    std::ifstream ifstr;
    ifstr.open(buffer);
    if (!ifstr.is_open()) continue;

    int pnum;
    std::string header;
    ifstr >> header >> pnum;
    std::cerr << image << ' ' << pnum << " patches" << std::endl;

    for (int p = 0; p < pnum; ++p) {
      Patch::PPatch ppatch(new Patch::CPatch());
      ifstr >> *ppatch;
      ppatch->_fix = 1;
      ppatch->_vimages.clear();

      image2index(*ppatch);
      if (ppatch->_images.empty()) continue;

      setGrids(*ppatch);
      addPatch(ppatch);
    }

    ifstr.close();
  }
}

void PMVS3::CPatchOrganizerS::collectPatches(const int target) {
  _ppatches.clear();

  for (int index = 0; index < _fm._tnum; ++index) {
    for (int i = 0; i < (int)_pgrids[index].size(); ++i) {
      auto begin = _pgrids[index][i].begin();
      while (begin != _pgrids[index][i].end()) {
        (*begin)->_id = -1;
        begin++;
      }
    }
  }

  int count = 0;
  for (int index = 0; index < _fm._tnum; ++index) {
    for (int i = 0; i < (int)_pgrids[index].size(); ++i) {
      auto begin = _pgrids[index][i].begin();
      while (begin != _pgrids[index][i].end()) {
        if ((*begin)->_id == -1) {
          (*begin)->_id = count++;

          if (target == 0 || (*begin)->_fix == 0) {
            _ppatches.push_back(*begin);
          }
        }
        ++begin;
      }
    }
  }
}

void PMVS3::CPatchOrganizerS::collectPatches(std::priority_queue<Patch::PPatch, std::vector<Patch::PPatch>, P_compare>& pqpatches) {
  for (int index = 0; index < _fm._tnum; ++index) {
    for (int i = 0; i < (int)_pgrids[index].size(); ++i) {
      auto begin = _pgrids[index][i].begin();
      while (begin != _pgrids[index][i].end()) {
        if ((*begin)->_flag == 0) {
          (*begin)->_flag = 1;
          pqpatches.push(*begin);
        }
        ++begin;
      }
    }
  }
}

void PMVS3::CPatchOrganizerS::collectPatches(const int index, std::priority_queue<Patch::PPatch, std::vector<Patch::PPatch>, P_compare>& pqpatches) {
  _fm._imageLocks[index].wrlock();
  for (int i = 0; i < (int)_pgrids[index].size(); ++i) {
    auto begin = _pgrids[index][i].begin();
    auto end = _pgrids[index][i].end();

    while (begin != end) {
      if ((*begin)->_images[0] == index && (*begin)->_flag == 0) {
        (*begin)->_flag = 1;
        pqpatches.push(*begin);
      }
      ++begin;
    }
  }
  _fm._imageLocks[index].unlock();
}

// Should be used only for writing
void PMVS3::CPatchOrganizerS::collectNonFixPatches(const int index, std::vector<Patch::PPatch>& ppatches) {
  _fm._imageLocks[index].wrlock();;
  for (int i = 0; i < (int)_pgrids[index].size(); ++i) {
    auto begin = _pgrids[index][i].begin();
    auto end = _pgrids[index][i].end();

    while (begin != end) {
      if ((*begin)->_images[0] == index && (*begin)->_fix == 0) {
        ppatches.push_back(*begin);
      }
      ++begin;
    }
  }
  _fm._imageLocks[index].unlock();
}

void PMVS3::CPatchOrganizerS::clearFlags(void) {
  auto bppatch = _ppatches.begin();
  auto eppatch = _ppatches.end();

  while (bppatch != eppatch) {
    (*bppatch)->_flag = 0;
    ++bppatch;
  }
}

void PMVS3::CPatchOrganizerS::clearCounts(void) {
  for (int index = 0; index < _fm._tnum; ++index) {
    auto begin = _counts[index].begin();
    auto end = _counts[index].end();
    while (begin != end) {
      *begin = (unsigned char)0;
      ++begin;
    }
  }
}

void PMVS3::CPatchOrganizerS::addPatch(Patch::PPatch& ppatch) {
  // First handle _vimages
  auto bimage = ppatch->_images.begin();
  auto eimage = ppatch->_images.end();
  auto bgrid = ppatch->_grids.begin();
  while (bimage != eimage) {
    const int index = *bimage;
    if (_fm._tnum <= index) {
      ++bimage;
      ++bgrid;
      continue;
    }

    const int index2 = (*bgrid)[1] * _gwidths[index] + (*bgrid)[0];
    _fm._imageLocks[index].wrlock();
    _pgrids[index][index2].push_back(ppatch);
	  _fm._imageLocks[index].unlock();

    ++bimage;
    ++bgrid;
  }

  // If depth, set vimages
  if (_fm._depth == 0) return;

  bimage = ppatch->_vimages.begin();
  eimage = ppatch->_vimages.end();
  bgrid = ppatch->_vgrids.begin();

  while (bimage != eimage) {
    const int index = *bimage;
    const int index2 = (*bgrid)[1] * _gwidths[index] + (*bgrid)[0];
	  _fm._imageLocks[index].wrlock();
    _vpgrids[index][index2].push_back(ppatch);
	  _fm._imageLocks[index].unlock();

    ++bimage;
    ++bgrid;
  }

  updateDepthMaps(ppatch);
}

void PMVS3::CPatchOrganizerS::updateDepthMaps(Patch::PPatch& ppatch) {
  for (int image = 0; image < _fm._tnum; ++image) {
    const Vec3f icoord = _fm._pss.project(image, ppatch->_coord, _fm._level);

    const float fx = icoord[0] / _fm._csize;
    const int xs[2] = {(int)floor(fx), (int)ceil(fx)};
    const float fy = icoord[1] / _fm._csize;
    const int ys[2] = {(int)floor(fy), (int)ceil(fy)};

    const float depth = _fm._pss._photos[image].OpticalAxis() * ppatch->_coord;

	  _fm._imageLocks[image].wrlock();
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
	      if (xs[i] < 0 || _gwidths[image] <= xs[i] || ys[j] < 0 || _gheights[image] <= ys[j]) continue;

        const int index = ys[j] * _gwidths[image] + xs[i];
	      if (_dpgrids[image][index] == _MAXDEPTH) {
          _dpgrids[image][index] = ppatch;
        } else {
          const float dtmp = _fm._pss._photos[image].OpticalAxis() * _dpgrids[image][index]->_coord;

          if (depth < dtmp) {
            _dpgrids[image][index] = ppatch;
          }
        }
      }
    }
    _fm._imageLocks[image].unlock();
  }
}

void PMVS3::CPatchOrganizerS::setGridsImages(Patch::CPatch& patch, const std::vector<int>& images) const {
  patch._images.clear();
  patch._grids.clear();
  std::vector<int>::const_iterator bimage = images.begin();
  std::vector<int>::const_iterator eimage = images.end();

  while (bimage != eimage) {
    const Vec3f icoord = _fm._pss.project(*bimage, patch._coord, _fm._level);
    const int ix = ((int)floor(icoord[0] + 0.5f)) / _fm._csize;
    const int iy = ((int)floor(icoord[1] + 0.5f)) / _fm._csize;
    if (0 <= ix && ix < _gwidths[*bimage] && 0 <= iy && iy < _gheights[*bimage]) {
      patch._images.push_back(*bimage);
      patch._grids.push_back(Vec2i(ix, iy));
    }
    ++bimage;
  }
}

void PMVS3::CPatchOrganizerS::setGrids(Patch::PPatch& ppatch) const {
  setGrids(*ppatch);
}

void PMVS3::CPatchOrganizerS::setGrids(Patch::CPatch& patch) const {
  patch._grids.clear();
  for (int i = 0; i < (int)patch._images.size(); ++i) {
    const int image = patch._images[i];
    Vec3f icoord = _fm._pss.project(image, patch._coord, _fm._level);
    const int ix = ((int)floor(icoord[0] + 0.5f)) / _fm._csize;
    const int iy = ((int)floor(icoord[1] + 0.5f)) / _fm._csize;
    patch._grids.push_back(TVec2<int>(ix, iy));
  }
}

void PMVS3::CPatchOrganizerS::setVImagesVGrids(Patch::PPatch& ppatch) {
  setVImagesVGrids(*ppatch);
}

void PMVS3::CPatchOrganizerS::setVImagesVGrids(Patch::CPatch& patch) {
  std::vector<int> used;
  used.resize(_fm._tnum);
  fill(used.begin(), used.end(), 0);

  auto bimage = patch._images.begin();
  auto eimage = patch._images.end();
  while (bimage != eimage) {
    if ((*bimage) < _fm._tnum) {
      used[*(bimage)] = 1;
    }
    ++bimage;
  }

  bimage = patch._vimages.begin();
  eimage = patch._vimages.end();
  while (bimage != eimage) {
    used[*(bimage++)] = 1;
  }

  for (int image = 0; image < _fm._tnum; ++image) {
    if (used[image]) continue;

    int ix, iy;
    if (isVisible0(patch, image, ix, iy, _fm._neighborThreshold, 1) == 0) continue;
    if (_fm._pss.getEdge(patch._coord, image, _fm._level) == 0) continue;

    patch._vimages.push_back(image);
    patch._vgrids.push_back(TVec2<int>(ix, iy));
  }
}

void PMVS3::CPatchOrganizerS::removePatch(const Patch::PPatch& ppatch) {
  for (int i = 0; i < (int)ppatch->_images.size(); ++i) {
    const int image = ppatch->_images[i];
    if (_fm._tnum <= image) continue;

    const int& ix = ppatch->_grids[i][0];
    const int& iy = ppatch->_grids[i][1];
    const int index = iy * _gwidths[image] + ix;
    _pgrids[image][index].erase(remove(_pgrids[image][index].begin(), _pgrids[image][index].end(), ppatch), _pgrids[image][index].end());
  }

  for (int i = 0; i < (int)ppatch->_vimages.size(); ++i) {
    const int image = ppatch->_vimages[i];
#ifdef DEBUG
    if (_fm._tnum <= image) {
      cerr << "Impossible in removePatch. _vimages must be targetting images" << std::endl;
      exit (1);
    }
#endif

    const int& ix = ppatch->_vgrids[i][0];
    const int& iy = ppatch->_vgrids[i][1];
    const int index = iy * _gwidths[image] + ix;
    _vpgrids[image][index].erase(remove(_vpgrids[image][index].begin(), _vpgrids[image][index].end(), ppatch), _vpgrids[image][index].end());
  }
}

int PMVS3::CPatchOrganizerS::isVisible0(const Patch::CPatch& patch, const int image, int& ix, int& iy, const float strict, const int lock) {
  const Vec3f icoord = _fm._pss.project(image, patch._coord, _fm._level);
  ix = ((int)floor(icoord[0] + 0.5f)) / _fm._csize;
  iy = ((int)floor(icoord[1] + 0.5f)) / _fm._csize;

  return isVisible(patch, image, ix, iy, strict, lock);
}

int PMVS3::CPatchOrganizerS::isVisible(const Patch::CPatch& patch, const int image, const int& ix, const int& iy, const float strict, const int lock) {
  const int& gwidth = _gwidths[image];
  const int& gheight = _gheights[image];

  if (ix < 0 || gwidth <= ix || iy < 0 || gheight <= iy) return 0;
  if (_fm._depth == 0) return 1;

  int ans = 0;
  Patch::PPatch dppatch = _MAXDEPTH;
  const int index = iy * gwidth + ix;

  if (lock) {
    _fm._imageLocks[image].rdlock();
  }

  if (_dpgrids[image][index] == _MAXDEPTH) {
    ans = 1;
  } else {
    dppatch = _dpgrids[image][index];
  }

  if (lock) {
    _fm._imageLocks[image].unlock();
  }

  if (ans == 1) {
    return 1;
  }

  Vec4f ray = patch._coord - _fm._pss._photos[image].OpticalCenter();
  unitize(ray);
  const float diff = ray * (patch._coord - dppatch->_coord);
  const float factor = std::min(2.0, 2.0 + ray * patch._normal);

  if (diff < _fm._optim.getUnit(image, patch._coord) * _fm._csize * strict * factor) {
    return 1;
  } else {
    return 0;
  }
}

void PMVS3::CPatchOrganizerS::findNeighbors(const Patch::CPatch& patch, std::vector<Patch::PPatch>& neighbors, const int lock, const float scale, const int margin, const int skipvis) {
  const float radius = 1.5 * margin * _fm._expand.computeRadius(patch);

  std::vector<int>::const_iterator bimage = patch._images.begin();
  std::vector<int>::const_iterator eimage = patch._images.end();
  std::vector<Vec2i>::const_iterator bgrid = patch._grids.begin();

#ifdef DEBUG
  if (patch._images.empty()) {
    cerr << "Empty patches in findCloses" << std::endl;
    exit (1);
  }
#endif
  float unit = 0.0f;
  for (int i = 0; i < (int)patch._images.size(); ++i) {
    unit += _fm._optim.getUnit(patch._images[i], patch._coord);
  }
  unit /= (int)patch._images.size();
  unit *= _fm._csize;

  while (bimage != eimage) {
    if (_fm._tnum <= *bimage) {
      ++bimage;
      ++bgrid;
      continue;
    }

    const int image = *bimage;
    const int& ix = (*bgrid)[0];
    const int& iy = (*bgrid)[1];
    if (lock) {
      _fm._imageLocks[image].rdlock();
    }

    for (int j = -margin; j <= margin; ++j) {
      const int ytmp = iy + j;
      if (ytmp < 0 || _fm._pos._gheights[image] <= ytmp) continue;

      for (int i = -margin; i <= margin; ++i) {
        const int xtmp = ix + i;

        if (xtmp < 0 || _fm._pos._gwidths[image] <= xtmp) continue;

        const int index = ytmp * _fm._pos._gwidths[image] + xtmp;
        std::vector<Patch::PPatch>::const_iterator bpatch = _fm._pos._pgrids[image][index].begin();
        std::vector<Patch::PPatch>::const_iterator epatch = _fm._pos._pgrids[image][index].end();
        while (bpatch != epatch) {
          if (_fm.isNeighborRadius(patch, **bpatch, unit, _fm._neighborThreshold * scale, radius))
            neighbors.push_back(*bpatch);
          ++bpatch;
        }
        bpatch = _fm._pos._vpgrids[image][index].begin();
        epatch = _fm._pos._vpgrids[image][index].end();
        while (bpatch != epatch) {
          if (_fm.isNeighborRadius(patch, **bpatch, unit, _fm._neighborThreshold * scale, radius)) {
            neighbors.push_back(*bpatch);
          }
          ++bpatch;
        }
      }
    }
    if (lock) {
      _fm._imageLocks[image].unlock();
    }

    ++bimage;
    ++bgrid;
  }

  if (skipvis == 0) {
    bimage = patch._vimages.begin();
    eimage = patch._vimages.end();
    bgrid = patch._vgrids.begin();

    while (bimage != eimage) {
      const int image = *bimage;
      const int& ix = (*bgrid)[0];
      const int& iy = (*bgrid)[1];

      if (lock) {
        _fm._imageLocks[image].rdlock();
      }

      for (int j = -margin; j <= margin; ++j) {
        const int ytmp = iy + j;
        if (ytmp < 0 || _fm._pos._gheights[image] <= ytmp) continue;

        for (int i = -margin; i <= margin; ++i) {
          const int xtmp = ix + i;
          if (xtmp < 0 || _fm._pos._gwidths[image] <= xtmp) continue;

          const int index = ytmp * _fm._pos._gwidths[image] + xtmp;
          std::vector<Patch::PPatch>::const_iterator bpatch = _fm._pos._pgrids[image][index].begin();
          std::vector<Patch::PPatch>::const_iterator epatch = _fm._pos._pgrids[image][index].end();

          while (bpatch != epatch) {
            if (_fm.isNeighborRadius(patch, **bpatch, unit, _fm._neighborThreshold * scale, radius)) {
              neighbors.push_back(*bpatch);
            }
            ++bpatch;
          }

          bpatch = _fm._pos._vpgrids[image][index].begin();
          epatch = _fm._pos._vpgrids[image][index].end();
          while (bpatch != epatch) {
            if (_fm.isNeighborRadius(patch, **bpatch, unit, _fm._neighborThreshold * scale, radius)) {
              neighbors.push_back(*bpatch);
            }
            ++bpatch;
          }
        }
      }
      if (lock) {
        _fm._imageLocks[image].unlock();
      }
      ++bimage;
      ++bgrid;
    }
  }

  std::sort(neighbors.begin(), neighbors.end());
  neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
}

float PMVS3::CPatchOrganizerS::computeUnit(const Patch::CPatch& patch) const {
  float unit = 0.0f;
  for (int i = 0; i < (int)patch._images.size(); ++i) {
    unit += _fm._optim.getUnit(patch._images[i], patch._coord);
  }
  unit /= (int)patch._images.size();
  unit *= _fm._csize;
  return unit;
}

void PMVS3::CPatchOrganizerS::setScales(Patch::CPatch& patch) const {
  const float unit = _fm._optim.getUnit(patch._images[0], patch._coord);
  const float unit2 = 2.0f * unit;
  Vec4f ray = patch._coord - _fm._pss._photos[patch._images[0]].OpticalCenter();
  unitize(ray);

  const int inum = std::min(_fm._tau, (int)patch._images.size());

  // First compute, how many pixel difference per unit along vertical
  //for (int i = 1; i < (int)patch._images.size(); ++i) {
  for (int i = 1; i < inum; ++i) {
    Vec3f diff = _fm._pss.project(patch._images[i], patch._coord, _fm._level) - _fm._pss.project(patch._images[i], patch._coord - ray * unit2, _fm._level);
    patch._dscale += norm(diff);
  }

  // set _dscale to the vertical distance where average pixel move is half pixel
  //patch._dscale /= (int)patch._images.size() - 1;
  patch._dscale /= inum - 1;
  patch._dscale = unit2 / patch._dscale;

  patch._ascale = atan(patch._dscale / (unit * _fm._wsize / 2.0f));
}

// write out results
void PMVS3::CPatchOrganizerS::writePLY(const std::vector<Patch::PPatch>& patches, const std::string filename) {
  std::ofstream ofstr;
  ofstr.open(filename.c_str());

  // Set output to full precision to be able to fully round-trip
  // the text representation with what we have in memory.
  ofstr << std::setprecision(std::numeric_limits<double>::max_digits10);

  ofstr << "ply" << '\n'
       << "format ascii 1.0" << '\n'
       << "element vertex " << (int)patches.size() << '\n'
       << "property float x" << '\n'
       << "property float y" << '\n'
       << "property float z" << '\n'
       << "property float nx" << '\n'
       << "property float ny" << '\n'
       << "property float nz" << '\n'
       << "property uchar diffuse_red" << '\n'
       << "property uchar diffuse_green" << '\n'
       << "property uchar diffuse_blue" << '\n'
       << "property float quality" << '\n'
       << "end_header" << '\n';

  std::vector<Patch::PPatch>::const_iterator bpatch = patches.begin();
  std::vector<Patch::PPatch>::const_iterator bend = patches.end();

  while (bpatch != bend) {
    // Get color
    Vec3i color;

    const int mode = 0;
    // 0: color from images
    // 1: fix
    // 2: angle
    if (mode == 0) {
      int denom = 0;
      Vec3f colorf;
      for (int i = 0; i < (int)(*bpatch)->_images.size(); ++i) {
        const int image = (*bpatch)->_images[i];
        colorf += _fm._pss.getColor((*bpatch)->_coord, image, _fm._level);
        denom++;
      }

      colorf /= denom;
      color[0] = std::min(255,(int)floor(colorf[0] + 0.5f));
      color[1] = std::min(255,(int)floor(colorf[1] + 0.5f));
      color[2] = std::min(255,(int)floor(colorf[2] + 0.5f));
    }
    else if (mode == 1) {
      if ((*bpatch)->_tmp == 1.0f) {
        color[0] = 255;
        color[1] = 0;
        color[2] = 0;
      } else {
        color[0] = 255;
        color[1] = 255;
        color[2] = 255;
      }
    } else if (mode == 2) {
      float angle = 0.0f;
      auto bimage = (*bpatch)->_images.begin();
      auto eimage = (*bpatch)->_images.end();

      while (bimage != eimage) {
        const int index = *bimage;
        Vec4f ray = _fm._pss._photos[index].OpticalCenter() - (*bpatch)->_coord;
        ray[3] = 0.0f;
        unitize(ray);

        angle += acos(ray * (*bpatch)->_normal);
        ++bimage;
      }

      angle = angle / (M_PI / 2.0f);
      float r, g, b;
      Image::CImage::gray2rgb(angle, r, g, b);
      color[0] = (int)(r * 255.0f);
      color[1] = (int)(g * 255.0f);
      color[2] = (int)(b * 255.0f);
    }

    ofstr << (*bpatch)->_coord[0] << ' '
          << (*bpatch)->_coord[1] << ' '
          << (*bpatch)->_coord[2] << ' '
          << (*bpatch)->_normal[0] << ' '
          << (*bpatch)->_normal[1] << ' '
          << (*bpatch)->_normal[2] << ' '
          << color[0] << ' ' << color[1] << ' ' << color[2] << ' '
          << (*bpatch)->_ncc << '\n';
      ++bpatch;
  }
  ofstr.close();
}

void PMVS3::CPatchOrganizerS::writePLY(const std::vector<Patch::PPatch>& patches, const std::string filename, const std::vector<Vec3i>& colors) {
  std::ofstream ofstr;
  ofstr.open(filename.c_str());

  // Set output to full precision to be able to fully round-trip
  // the text representation with what we have in memory.
  ofstr << std::setprecision(std::numeric_limits<double>::max_digits10);

  ofstr << "ply" << '\n'
       << "format ascii 1.0" << '\n'
       << "element vertex " << (int)patches.size() << '\n'
       << "property float x" << '\n'
       << "property float y" << '\n'
       << "property float z" << '\n'
       << "property float nx" << '\n'
       << "property float ny" << '\n'
       << "property float nz" << '\n'
       << "property uchar diffuse_red" << '\n'
       << "property uchar diffuse_green" << '\n'
       << "property uchar diffuse_blue" << '\n'
       << "end_header" << '\n';

  std::vector<Patch::PPatch>::const_iterator bpatch = patches.begin();
  std::vector<Patch::PPatch>::const_iterator bend = patches.end();
  std::vector<Vec3i>::const_iterator colorb = colors.begin();

  while (bpatch != bend) {
    ofstr << (*bpatch)->_coord[0] << ' '
          << (*bpatch)->_coord[1] << ' '
          << (*bpatch)->_coord[2] << ' '
          << (*bpatch)->_normal[0] << ' '
          << (*bpatch)->_normal[1] << ' '
          << (*bpatch)->_normal[2] << ' '
          << *colorb << '\n';
    ++bpatch;
    ++colorb;
  }
  ofstr.close();
}
