#define _USE_MATH_DEFINES
#include <cmath>

#include <ctime>
#include <map>

#include "pmvs/detectFeatures.hpp"
#include "pmvs/findMatch.hpp"

PMVS3::CFindMatch::CFindMatch() : _pos(*this), _seed(*this), _expand(*this), _filter(*this), _optim(*this) {
  _debug = 0;
}

PMVS3::CFindMatch::~CFindMatch() {
  for (int image = 0; image < _num; ++image) {
    delete _imageLocks[image];
    delete _countLocks[image];
  }
  _imageLocks.clear();
  _countLocks.clear();
}

void PMVS3::CFindMatch::updateThreshold() {
  _nccThreshold -= 0.05f;
  _nccThresholdBefore -= 0.05f;

  _countThreshold1 = 2;
}

void PMVS3::CFindMatch::init(const SOption& option) {
  _timages = option._timages;
  _oimages = option._oimages;
  _images.clear();
  _images.insert(_images.end(), _timages.begin(), _timages.end());
  _images.insert(_images.end(), _oimages.begin(), _oimages.end());

  _tnum = (int)_timages.size();
  _num = (int)_images.size();

  _prefix = option._prefix;
  _level = option._level;
  _csize = option._csize;
  _nccThreshold = option._threshold;
  _wsize = option._wsize;
  _minImageNumThreshold = option._minImageNum;
  _CPU = option._CPU;
  _setEdge = option._setEdge;
  _sequenceThreshold = option._sequence;

  _junit = 100;
  // This initialization does not matter
  _visibleThreshold = 0.0f;
  _visibleThresholdLoose = 0.0f;

  //_tau = max(option._minImageNum * 2, min(_num, 5));
  _tau = std::min(option._minImageNum * 2, _num);

  _depth = 0;

  // set target images and other images
  _bindexes = option._bindexes;
  _visdata = option._visdata;
  _visdata2 = option._visdata2;

  //----------------------------------------------------------------------
  // mtx_init(&_lock, mtx_plain | mtx_recursive);
  for (int image = 0; image < _num; ++image) {
    _imageLocks.push_back(new std::shared_mutex());
    _countLocks.push_back(new std::shared_mutex());
  }
  // We set _level + 3, to use multi-resolutional texture grabbing
  _pss.init(_images, _prefix, _level + 3, _wsize, 1);

  if (_setEdge != 0.0f) {
    _pss.setEdge(_setEdge);
  }
  _pss.setDistances();

  // Detect features if not yet done
  CDetectFeatures df;
  const int fcsize = 16;
  df.run(_pss, _num, fcsize, _level, _CPU);

  // Initialize each core member. _po should be first
  _pos.init();
  _seed.init(df._points);
  _expand.init();
  _filter.init();
  _optim.init();
  //----------------------------------------------------------------------
  // Init thresholds
  _angleThreshold0 = 60.0f * M_PI / 180.0f;
  _angleThreshold1 = 60.0f * M_PI / 180.0f;

  _countThreshold0 = 2;
  _countThreshold1 = 4;
  _countThreshold2 = 2;

  _neighborThreshold0 = 0.5f;
  _neighborThreshold1 = 1.0f;
  _neighborThreshold2 = 1.0f;
  
  _maxAngleThreshold = option._maxAngleThreshold;
  _nccThresholdBefore = _nccThreshold - 0.3f;
  _quadThreshold = option._quadThreshold;
  _epThreshold = 2.0f;
}

int PMVS3::CFindMatch::insideBimages(const Vec4f& coord) const {
  for (int i = 0; i < (int)_bindexes.size(); ++i) {
    const int index = _bindexes[i];
    const Vec3f icoord = _pss.project(index, coord, _level);
    if (icoord[0] < 0.0 || _pss.getWidth(index, _level) - 1 < icoord[0] || icoord[1] < 0.0 || _pss.getHeight(index, _level) - 1 < icoord[1]) {
      return 0;
    }
  }
  return 1;
}

int PMVS3::CFindMatch::isNeighbor(const Patch::CPatch& lhs, const Patch::CPatch& rhs, const float neighborThreshold) const {
  const float hunit = (_optim.getUnit(lhs._images[0], lhs._coord) + _optim.getUnit(rhs._images[0], rhs._coord)) / 2.0 * _csize;
  return isNeighbor(lhs, rhs, hunit, neighborThreshold);
}

int PMVS3::CFindMatch::isNeighbor(const Patch::CPatch& lhs, const Patch::CPatch& rhs, const float hunit, const float neighborThreshold) const {
  if (lhs._normal * rhs._normal < cos(120.0 * M_PI / 180.0)) {
    return 0;
  }

  const Vec4f diff = rhs._coord - lhs._coord;
  const float vunit = lhs._dscale + rhs._dscale;

  const float f0 = lhs._normal * diff;
  const float f1 = rhs._normal * diff;
  float ftmp = (fabs(f0) + fabs(f1)) / 2.0;
  ftmp /= vunit;

  // this may loosen the isneighbor testing. need to tighten (decrease) threshold?
  const float hsize = norm(2 * diff - lhs._normal * f0 - rhs._normal * f1) / 2.0 / hunit;
  if (1.0 < hsize) {
    ftmp /= std::min(2.0f, hsize);
  }

  if (ftmp < neighborThreshold) {
    return 1;
  }

  return 0;
}

int PMVS3::CFindMatch::isNeighborRadius(
  const Patch::CPatch& lhs,
  const Patch::CPatch& rhs,
  const float hunit,
  const float neighborThreshold,
  const float radius
) const {
  if (lhs._normal * rhs._normal < cos(120.0 * M_PI / 180.0)) {
    return 0;
  }

  const Vec4f diff = rhs._coord - lhs._coord;
  const float vunit = lhs._dscale + rhs._dscale;

  const float f0 = lhs._normal * diff;
  const float f1 = rhs._normal * diff;
  float ftmp = (fabs(f0) + fabs(f1)) / 2.0;
  ftmp /= vunit;

  // this may loosen the isneighbor testing. need to tighten (decrease) threshold?
  const float hsize = norm(2 * diff - lhs._normal * f0 - rhs._normal * f1) / 2.0 / hunit;

  // radius check
  if (radius / hunit < hsize) return 0;

  if (1.0 < hsize) {
    ftmp /= std::min(2.0f, hsize);
  }

  if (ftmp < neighborThreshold) {
    return 1;
  } else {
    return 0;
  }
}

void PMVS3::CFindMatch::run(void) {
  time_t tv;
  time(&tv);
  time_t curtime = tv;

  // Seed generation
  _seed.run();
  _seed.clear();

  ++_depth;
  _pos.collectPatches();

  // Expansion
  const int TIME = 3;
  for (int t = 0; t < TIME; ++t) {
    _expand.run();
    _filter.run();

    updateThreshold();

    std::cout << "STATUS: ";
    for (int i = 0; i < (int)_optim._status.size(); ++i) {
      std::cout << _optim._status[i] << ' ';
      if (i % 10 == 9) {
        std::cout << std::endl;
      }
    }
    std::cout << std::endl;

    ++_depth;
  }
  time(&tv);
  std::cerr << "---- Total: " << (tv - curtime)/CLOCKS_PER_SEC << " secs ----" << std::endl;
}

void PMVS3::CFindMatch::write(const std::string prefix, bool bExportPLY, bool bExportPatch, bool bExportPSet) {
  _pos.writePatches2(prefix, bExportPLY, bExportPatch, bExportPSet);
}
