#define _USE_MATH_DEFINES
#include <cmath>

#include <algorithm>
#include <numeric>
#include "pmvs/findMatch.hpp"
#include "pmvs/optim.hpp"
#include <cstdio>

#include "nlopt.hpp"

PMVS3::COptim* PMVS3::COptim::_one = NULL;

PMVS3::COptim::COptim(CFindMatch& findMatch) : _fm(findMatch) {
  _one = this;
  _status.resize(35);
  fill(_status.begin(), _status.end(), 0);
}

void PMVS3::COptim::init(void) {
  _vect0T.resize(_fm._CPU);
  _centersT.resize(_fm._CPU);
  _raysT.resize(_fm._CPU);
  _indexesT.resize(_fm._CPU);
  _dscalesT.resize(_fm._CPU);
  _ascalesT.resize(_fm._CPU);
  _paramsT.resize(_fm._CPU);

  _texsT.resize(_fm._CPU);
  _weightsT.resize(_fm._CPU);

  for (int c = 0; c < _fm._CPU; ++c) {
    _texsT[c].resize(_fm.NumImages());
    _weightsT[c].resize(_fm.NumImages());
    for (int j = 0; j < _fm._tau; ++j) {
      _texsT[c][j].resize(3 * _fm._wsize * _fm._wsize);
    }
  }

  setAxesScales();
}

void PMVS3::COptim::setAxesScales(void) {
  _xaxes.resize(_fm.NumImages());
  _yaxes.resize(_fm.NumImages());
  _zaxes.resize(_fm.NumImages());
  for (int index = 0; index < _fm.NumImages(); ++index) {
    _zaxes[index] = Vec3f(_fm._pss._photos[index].OpticalAxis()[0], _fm._pss._photos[index].OpticalAxis()[1], _fm._pss._photos[index].OpticalAxis()[2]);
    _xaxes[index] = Vec3f(_fm._pss._photos[index].ProjectionMatrix()[0][0][0], _fm._pss._photos[index].ProjectionMatrix()[0][0][1], _fm._pss._photos[index].ProjectionMatrix()[0][0][2]);
    _yaxes[index] = cross(_zaxes[index], _xaxes[index]);
    unitize(_yaxes[index]);
    _xaxes[index] = cross(_yaxes[index], _zaxes[index]);
  }

  _ipscales.resize(_fm.NumImages());
  for (int index = 0; index < _fm.NumImages(); ++index) {
    const Vec4f xaxe(_xaxes[index][0], _xaxes[index][1], _xaxes[index][2], 0.0);
    const Vec4f yaxe(_yaxes[index][0], _yaxes[index][1], _yaxes[index][2], 0.0);

    const float fx = xaxe * _fm._pss._photos[index].ProjectionMatrix()[0][0];
    const float fy = yaxe * _fm._pss._photos[index].ProjectionMatrix()[0][1];
    _ipscales[index] = fx + fy;
  }
}

void PMVS3::COptim::collectImages(const int index, std::vector<int>& indexes) const{
  // Find images with constraints _angleThreshold, _visdata,
  // _sequenceThreshold, _targets. Results are sorted by
  // CphotoSet::_distances.
  indexes.clear();
  Vec4f ray0 = _fm._pss._photos[index].OpticalAxis();
  ray0[3] = 0.0f;

  std::vector<Vec2f> candidates;
  // Search for only related images
  for (int i = 0; i < (int)_fm._visdata2[index].size(); ++i) {
    const int indextmp = _fm._visdata2[index][i];

    if (_fm._sequenceThreshold != -1 && _fm._sequenceThreshold < abs(index - indextmp)) continue;

    Vec4f ray1 = _fm._pss._photos[indextmp].OpticalAxis();
    ray1[3] = 0.0f;

    if (ray0 * ray1 < cos(_fm._angleThreshold0)) continue;

    candidates.push_back(Vec2f(_fm._pss._distances[index][indextmp], indextmp));
  }

  sort(candidates.begin(), candidates.end(), Svec2cmp<float>());
  for (int i = 0; i < std::min(_fm._tau, (int)candidates.size()); ++i) {
    indexes.push_back((int)candidates[i][1]);
  }
}

int PMVS3::COptim::preProcess(Patch::CPatch& patch, const int id, const int seed) {
  addImages(patch);

  // Here define reference images, and sort images.
  // Something similar to constraintImages is done inside.
  constraintImages(patch, _fm._nccThresholdBefore, id);

  // Fix the reference image and sort the other  _tau - 1 images.
  sortImages(patch);

  // Pierre Moulon (it avoid crash in some case)
  if( (int)patch._images.size() > 0) {
    // setSscales should be here to avoid noisy output
    _fm._pos.setScales(patch);
  }

  // Check minimum number of images
  if ((int)patch._images.size() < _fm._minImageNumThreshold) return 1;

  const int flag = _fm._pss.checkAngles(patch._coord, patch._images, _fm._maxAngleThreshold, _fm._angleThreshold1, _fm._minImageNumThreshold);

  if (flag) {
    patch._images.clear();
    return 1;
  }

  return 0;
}

void PMVS3::COptim::filterImagesByAngle(Patch::CPatch& patch) {
  std::vector<int> newindexes;

  auto bimage = patch._images.begin();
  auto eimage = patch._images.end();

  while (bimage != eimage) {
    const int index = *bimage;
    Vec4f ray = _fm._pss._photos[index].OpticalCenter() - patch._coord;
    unitize(ray);
    if (ray * patch._normal < cos(_fm._angleThreshold1)) {
      // if reference image is deleted, over
      if (bimage == patch._images.begin()) {
        patch._images.clear();
        return;
      }
    }
    else {
      newindexes.push_back(index);
    }
    ++bimage;
  }

  patch._images.swap(newindexes);
}

int PMVS3::COptim::postProcess(Patch::CPatch& patch, const int id, const int seed) {
  if ((int)patch._images.size() < _fm._minImageNumThreshold) return 1;

  if (_fm._pss.getMask(patch._coord, _fm._level) == 0 || _fm.insideBimages(patch._coord) == 0) return 1;

  addImages(patch);

  constraintImages(patch, _fm._nccThreshold, id);
  filterImagesByAngle(patch);

  if ((int)patch._images.size() < _fm._minImageNumThreshold) return 1;

  _fm._pos.setGrids(patch);

  setRefImage(patch, id);
  constraintImages(patch, _fm._nccThreshold, id);

  if ((int)patch._images.size() < _fm._minImageNumThreshold) return 1;

  _fm._pos.setGrids(patch);

  // set _timages
  patch._timages = 0;
  std::vector<int>::const_iterator begin = patch._images.begin();
  std::vector<int>::const_iterator end = patch._images.end();
  while (begin != end) {
    if (*begin < _fm.NumTargetImages()) {
      ++patch._timages;
    }
    ++begin;
  }

  patch._tmp = patch.score2(_fm._nccThreshold);
  // Set vimages vgrids.
  if (_fm._depth) {
    _fm._pos.setVImagesVGrids(patch);

    if (2 <= _fm._depth && check(patch)) return 1;
  }
  return 0;
}

void PMVS3::COptim::constraintImages(Patch::CPatch& patch, const float nccThreshold, const int id) {
  std::vector<float> inccs;
  setINCCs(patch, inccs, patch._images, id, 0);

  //----------------------------------------------------------------------
  // Constraint images
  std::vector<int> newimages;
  newimages.push_back(patch._images[0]);
  for (int i = 1; i < (int)patch._images.size(); ++i) {
    if (inccs[i] < 1.0f - nccThreshold) {
      newimages.push_back(patch._images[i]);
    }
  }
  patch._images.swap(newimages);
}

void PMVS3::COptim::setRefImage(Patch::CPatch& patch, const int id) {
#ifdef DEBUG
  if (patch._images.empty()) {
    cerr << "empty images" << endl;    exit (1);
  }
#endif
  //----------------------------------------------------------------------
  // Set the reference image
  // Only for target images
  std::vector<int> indexes;
  std::vector<int>::const_iterator begin = patch._images.begin();
  std::vector<int>::const_iterator end = patch._images.end();
  while (begin != end) {
    if (*begin < _fm.NumTargetImages()) {
      indexes.push_back(*begin);
    }
    ++begin;
  }
  // To avoid segmentation error on alley dataset. (this code is necessary because of the use of filterExact)
  if (indexes.empty()) {
    patch._images.clear();
    return;
  }

  std::vector<std::vector<float> > inccs;
  setINCCs(patch, inccs, indexes, id, 1);

  int refindex = -1;
  float refncc = INT_MAX/2;
  for (int i = 0; i < (int)indexes.size(); ++i) {
    const float sum = accumulate(inccs[i].begin(), inccs[i].end(), 0.0f);
    if (sum < refncc) {
      refncc = sum;
      refindex = i;
    }
  }

  const int refIndex = indexes[refindex];
  for (int i = 0; i < (int)patch._images.size(); ++i) {
    if (patch._images[i] == refIndex) {
      const int itmp = patch._images[0];
      patch._images[0] = refIndex;
      patch._images[i] = itmp;
      break;
    }
  }
}

// When no sampling was done, this is used
void PMVS3::COptim::setRefConstraintImages(Patch::CPatch& patch, const float nccThreshold, const int id) {
  //----------------------------------------------------------------------
  // Set the reference image
  std::vector<std::vector<float> > inccs;
  setINCCs(patch, inccs, patch._images, id, 1);

  int refindex = -1;
  float refncc = INT_MAX/2;
  for (int i = 0; i < (int)patch._images.size(); ++i) {
    const float sum = accumulate(inccs[i].begin(), inccs[i].end(), 0.0f);
    if (sum < refncc) {
      refncc = sum;
      refindex = i;
    }
  }

  const float robustThreshold = robustincc(1.0f - nccThreshold);
  std::vector<int> newimages;
  newimages.push_back(patch._images[refindex]);
  for (int i = 0; i < (int)patch._images.size(); ++i) {
    if (i != refindex && inccs[refindex][i] < robustThreshold) {
      newimages.push_back(patch._images[i]);
    }
  }
  patch._images.swap(newimages);
}

void PMVS3::COptim::sortImages(Patch::CPatch& patch) const{
  const int newm = 1;
  if (newm == 1) {
    const float threshold = 1.0f - cos(10.0 * M_PI / 180.0);
    std::vector<int> indexes, indexes2;
    std::vector<float> units, units2;
    std::vector<Vec4f> rays, rays2;

    computeUnits(patch, indexes, units, rays);

    patch._images.clear();
    if (indexes.size() < 2) return;

    units[0] = 0.0f;

    while (!indexes.empty()) {
      auto ite = min_element(units.begin(), units.end());
      const int index = ite - units.begin();

      patch._images.push_back(indexes[index]);

      // Remove other images within 5 degrees
      indexes2.clear();
      units2.clear();
      rays2.clear();
      for (int j = 0; j < (int)rays.size(); ++j) {
        if (j == index) continue;

        indexes2.push_back(indexes[j]);
        rays2.push_back(rays[j]);
        const float ftmp = std::min(threshold, std::max(threshold / 2.0f, 1.0f - rays[index] * rays[j]));

        units2.push_back(units[j] * (threshold / ftmp));
      }
      indexes2.swap(indexes);
      units2.swap(units);
      rays2.swap(rays);
    }
  } else {
    //----------------------------------------------------------------------
    //Sort and grab the best _tau images. All the other images don't
    //matter.  First image is the reference and fixed
    const float threshold = cos(5.0 * M_PI / 180.0);
    std::vector<int> indexes, indexes2;
    std::vector<float> units, units2;
    std::vector<Vec4f> rays, rays2;

    computeUnits(patch, indexes, units, rays);

    patch._images.clear();
    if (indexes.size() < 2) return;

    units[0] = 0.0f;

    while (!indexes.empty()) {
      //for (int i = 0; i < size; ++i) {
      auto ite = min_element(units.begin(), units.end());
      const int index = ite - units.begin();

      patch._images.push_back(indexes[index]);

      // Remove other images within 5 degrees
      indexes2.clear();
      units2.clear();
      rays2.clear();
      for (int j = 0; j < (int)rays.size(); ++j) {
        if (rays[index] * rays[j] < threshold) {
          indexes2.push_back(indexes[j]);
          units2.push_back(units[j]);
          rays2.push_back(rays[j]);
        }
      }
      indexes2.swap(indexes);
      units2.swap(units);
      rays2.swap(rays);
    }
  }
}

int PMVS3::COptim::check(Patch::CPatch& patch) {
  const float gain = _fm._filter.computeGain(patch, 1);
  patch._tmp = gain;

  if (gain < 0.0) {
    patch._images.clear();
    return 1;
  }

  {
    std::vector<Patch::PPatch> neighbors;
    _fm._pos.findNeighbors(patch, neighbors, 1, 4, 2);
    // Only check when enough number of neighbors
    if (6 < (int)neighbors.size() && _fm._filter.filterQuad(patch, neighbors)) {
      patch._images.clear();
      return 1;
    }
  }

  return 0;
}

void PMVS3::COptim::removeImagesEdge(Patch::CPatch& patch) const{
  std::vector<int> newindexes;
  std::vector<int>::const_iterator bimage = patch._images.begin();
  std::vector<int>::const_iterator eimage = patch._images.end();
  while (bimage != eimage) {
    if (_fm._pss.getEdge(patch._coord, *bimage, _fm._level)) {
      newindexes.push_back(*bimage);
    }
    ++bimage;
  }
  patch._images.swap(newindexes);
}

void PMVS3::COptim::addImages(Patch::CPatch& patch) const{
  // take into account _edge
  std::vector<int> used;
  used.resize(_fm.NumImages());
  for (int index = 0; index < _fm.NumImages(); ++index) {
    used[index] = 0;
  }

  std::vector<int>::const_iterator bimage = patch._images.begin();
  std::vector<int>::const_iterator eimage = patch._images.end();
  while (bimage != eimage) {
    used[*bimage] = 1;
    ++bimage;
  }

  bimage = _fm._visdata2[patch._images[0]].begin();
  eimage = _fm._visdata2[patch._images[0]].end();

  const float athreshold = cos(_fm._angleThreshold0);
  while (bimage != eimage) {
    if (used[*bimage]) {
      ++bimage;
      continue;
    }

    const Vec3f icoord = _fm._pss.project(*bimage, patch._coord, _fm._level);
    if (icoord[0] < 0.0f || _fm._pss.getWidth(*bimage, _fm._level) - 1 <= icoord[0] || icoord[1] < 0.0f || _fm._pss.getHeight(*bimage, _fm._level) - 1 <= icoord[1]) {
      ++bimage;
      continue;
    }

    if (_fm._pss.getEdge(patch._coord, *bimage, _fm._level) == 0) {
      ++bimage;
      continue;
    }

    Vec4f ray = _fm._pss._photos[*bimage].OpticalCenter() - patch._coord;
    unitize(ray);
    const float ftmp = ray * patch._normal;

    if (athreshold <= ftmp) {
      patch._images.push_back(*bimage);
    }

    ++bimage;
  }
}

void PMVS3::COptim::computeUnits(const Patch::CPatch& patch, std::vector<float>& units) const{
  const int size = (int)patch._images.size();
  units.resize(size);

  std::vector<int>::const_iterator bimage = patch._images.begin();
  std::vector<int>::const_iterator eimage = patch._images.end();

  auto bfine = units.begin();

  while (bimage != eimage) {
    *bfine = INT_MAX/2;

    *bfine = getUnit(*bimage, patch._coord);
    Vec4f ray = _fm._pss._photos[*bimage].OpticalCenter() - patch._coord;
    unitize(ray);
    const float denom = ray * patch._normal;
    if (0.0 < denom) {
      *bfine /= denom;
    } else {
      *bfine = INT_MAX/2;
    }

    ++bimage;
    ++bfine;
  }
}

void PMVS3::COptim::computeUnits(const Patch::CPatch& patch, std::vector<int>& indexes, std::vector<float>& units, std::vector<Vec4f>& rays) const{
  std::vector<int>::const_iterator bimage = patch._images.begin();
  std::vector<int>::const_iterator eimage = patch._images.end();

  while (bimage != eimage) {
    Vec4f ray = _fm._pss._photos[*bimage].OpticalCenter() - patch._coord;
    unitize(ray);
    const float dot = ray * patch._normal;
    if (dot <= 0.0f) {
      ++bimage;
      continue;
    }

    const float scale = getUnit(*bimage, patch._coord);
    const float fine = scale / dot;

    indexes.push_back(*bimage);
    units.push_back(fine);
    rays.push_back(ray);
    ++bimage;
  }
}

void PMVS3::COptim::refinePatch(Patch::CPatch& patch, const int id, const int time) {
  if(!refinePatchBFGS(patch, id, 1000, 1)) {
	  std::cout << "refinePatchBFGS failed!" << std::endl;
  }

  if (patch._images.empty()) return;
}

//----------------------------------------------------------------------
// BFGS functions
//----------------------------------------------------------------------
double PMVS3::COptim::my_f(unsigned n, const double *x, double *grad, void *my_func_data)
{
  double xs[3] = {x[0], x[1], x[2]};
  const int id = *((int*)my_func_data);

  const float angle1 = xs[1] * _one->_ascalesT[id];
  const float angle2 = xs[2] * _one->_ascalesT[id];

  double ret = 0.0;

  //?????
  const double bias = 0.0f;//2.0 - exp(- angle1 * angle1 / sigma2) - exp(- angle2 * angle2 / sigma2);

  Vec4f coord, normal;
  _one->decode(coord, normal, xs, id);

  const int index = _one->_indexesT[id][0];
  Vec4f pxaxis, pyaxis;
  _one->getPAxes(index, coord, normal, pxaxis, pyaxis);

  const int size = std::min(_one->_fm._tau, (int)_one->_indexesT[id].size());
  const int mininum = std::min(_one->_fm._minImageNumThreshold, size);

  for (int i = 0; i < size; ++i) {
    int flag = _one->grabTex(coord, pxaxis, pyaxis, normal, _one->_indexesT[id][i], _one->_fm._wsize, _one->_texsT[id][i]);

    if (flag == 0) {
      _one->normalize(_one->_texsT[id][i]);
    }
  }

  const int pairwise = 0;
  if (pairwise) {
    double ans = 0.0f;
    int denom = 0;
    for (int i = 0; i < size; ++i) {
      for (int j = i+1; j < size; ++j) {
        if (_one->_texsT[id][i].empty() || _one->_texsT[id][j].empty()) continue;

        ans += robustincc(1.0 - _one->dot(_one->_texsT[id][i], _one->_texsT[id][j]));
        denom++;
      }
    }

    if (denom < mininum * (mininum - 1) / 2) {
      ret = 2.0f;
    } else {
      ret = ans / denom + bias;
    }
  } else {
    if (_one->_texsT[id][0].empty()) {
      return 2.0;
    }

    double ans = 0.0f;
    int denom = 0;
    for (int i = 1; i < size; ++i) {
      if (_one->_texsT[id][i].empty()) continue;

      ans += robustincc(1.0 - _one->dot(_one->_texsT[id][0], _one->_texsT[id][i]));
      denom++;
    }

    if (denom < mininum - 1) {
      ret = 2.0f;
    } else {
      ret = ans / denom + bias;
    }
  }

  return ret;
}

bool PMVS3::COptim::refinePatchBFGS(Patch::CPatch& patch, const int id, const int time, const int ncc)
{
  int idtmp = id;

  _centersT[id] = patch._coord;
  _raysT[id] = patch._coord - _fm._pss._photos[patch._images[0]].OpticalCenter();
  unitize(_raysT[id]);
  _indexesT[id] = patch._images;

  _dscalesT[id] = patch._dscale;
  _ascalesT[id] = M_PI / 48.0f;//patch._ascale;

  computeUnits(patch, _weightsT[id]);
  for (int i = 1; i < (int)_weightsT[id].size(); ++i) {
    _weightsT[id][i] = std::min(1.0f, _weightsT[id][0] / _weightsT[id][i]);
  }
  _weightsT[id][0] = 1.0f;

  double p[3];
  encode(patch._coord, patch._normal, p, id);

  double min_angle = -23.99999;	//(- M_PI / 2.0) / _one->_ascalesT[id];
  double max_angle = 23.99999;	//(M_PI / 2.0) / _one->_ascalesT[id];

  std::vector<double> lower_bounds(3);
  lower_bounds[0] = -HUGE_VAL;		// Not bound
  lower_bounds[1] = min_angle;
  lower_bounds[2] = min_angle;
  std::vector<double> upper_bounds(3);
  upper_bounds[0] = HUGE_VAL;		// Not bound
  upper_bounds[1] = max_angle;
  upper_bounds[2] = max_angle;

  bool success = false;

  try {
    // LN_NELDERMEAD: Corresponds to the N-Simplex-Algorithm of GSL, that was used originally here
    // LN_SBPLX
    // LN_COBYLA
    // LN_BOBYQA
    // LN_PRAXIS
    nlopt::opt opt(nlopt::LN_BOBYQA, 3);
    opt.set_min_objective(my_f, &idtmp);
    opt.set_xtol_rel(1.e-7);
    opt.set_maxeval(time);

    opt.set_lower_bounds(lower_bounds);
    opt.set_upper_bounds(upper_bounds);

    std::vector<double> x(3);
    for(int i = 0; i < 3; i++)
    {
      // NLOPT returns an error if x is not within the bounds
      x[i] = std::max(std::min(p[i], upper_bounds[i]), lower_bounds[i]);
    }

    double minf;
    nlopt::srand(1);
    nlopt::result result = opt.optimize(x, minf);

    p[0] = x[0];
    p[1] = x[1];
    p[2] = x[2];

    success = (result == nlopt::SUCCESS || result == nlopt::STOPVAL_REACHED || result == nlopt::FTOL_REACHED || result == nlopt::XTOL_REACHED);
  } catch(std::exception &e) {
	  success = false;
  }

  if (success) {
    decode(patch._coord, patch._normal, p, id);

    patch._ncc = 1.0 - unrobustincc(computeINCC(patch._coord, patch._normal, patch._images, id, 1));
  } else {
    return false;
  }

  return true;
}

void PMVS3::COptim::encode(const Vec4f& coord, double* const vect, const int id) const {
  vect[0] = (coord - _centersT[id]) * _raysT[id] / _dscalesT[id];
}

void PMVS3::COptim::encode(const Vec4f& coord, const Vec4f& normal, double* const vect, const int id) const {
  encode(coord, vect, id);

  const int image = _indexesT[id][0];
  const float fx = _xaxes[image] * proj(normal); // projects from 4D to 3D, divide by last value
  const float fy = _yaxes[image] * proj(normal);
  const float fz = _zaxes[image] * proj(normal);

  vect[2] = asin(std::max(-1.0f, std::min(1.0f, fy)));
  const float cosb = cos(vect[2]);

  if (cosb == 0.0) {
    vect[1] = 0.0;
  } else {
    const float sina = fx / cosb;
    const float cosa = - fz / cosb;
    vect[1] = acos(std::max(-1.0f, std::min(1.0f, cosa)));
    if (sina < 0.0) {
      vect[1] = - vect[1];
    }
  }

  vect[1] = vect[1] / _ascalesT[id];
  vect[2] = vect[2] / _ascalesT[id];
}

void PMVS3::COptim::decode(Vec4f& coord, Vec4f& normal, const double* const vect, const int id) const {
  decode(coord, vect, id);
  const int image = _indexesT[id][0];

  const float angle1 = vect[1] * _ascalesT[id];
  const float angle2 = vect[2] * _ascalesT[id];

  const float fx = sin(angle1) * cos(angle2);
  const float fy = sin(angle2);
  const float fz = - cos(angle1) * cos(angle2);

  Vec3f ftmp = _xaxes[image] * fx + _yaxes[image] * fy + _zaxes[image] * fz;
  normal = Vec4f(ftmp[0], ftmp[1], ftmp[2], 0.0f);
}

void PMVS3::COptim::decode(Vec4f& coord, const double* const vect, const int id) const {
  coord = _centersT[id] + _dscalesT[id] * vect[0] * _raysT[id];
}

void PMVS3::COptim::setINCCs(const Patch::CPatch& patch, std::vector<float> & inccs, const std::vector<int>& indexes, const int id, const int robust) {
  const int index = indexes[0];
  Vec4f pxaxis, pyaxis;
  getPAxes(index, patch._coord, patch._normal, pxaxis, pyaxis);

  std::vector<std::vector<float> >& texs = _texsT[id];

  const int size = (int)indexes.size();
  for (int i = 0; i < size; ++i) {
    const int flag = grabTex(patch._coord, pxaxis, pyaxis, patch._normal, indexes[i], _fm._wsize, texs[i]);
    if (flag == 0) {
      normalize(texs[i]);
    }
  }

  inccs.resize(size);
  if (texs[0].empty()) {
    fill(inccs.begin(), inccs.end(), 2.0f);
    return;
  }

  for (int i = 0; i < size; ++i) {
    if (i == 0) {
      inccs[i] = 0.0f;
    } else if (!texs[i].empty()) {
      if (robust == 0) {
        inccs[i] = 1.0f - dot(texs[0], texs[i]);
      } else {
        inccs[i] = robustincc(1.0f - dot(texs[0], texs[i]));
      }
    }
    else {
      inccs[i] = 2.0f;
    }
  }
}

void PMVS3::COptim::setINCCs(const Patch::CPatch& patch, std::vector<std::vector<float> >& inccs, const std::vector<int>& indexes, const int id, const int robust) {
  const int index = indexes[0];
  Vec4f pxaxis, pyaxis;
  getPAxes(index, patch._coord, patch._normal, pxaxis, pyaxis);

  std::vector<std::vector<float>>& texs = _texsT[id];

  const int size = (int)indexes.size();
  for (int i = 0; i < size; ++i) {
    const int flag = grabTex(patch._coord, pxaxis, pyaxis, patch._normal, indexes[i], _fm._wsize, texs[i]);

    if (flag == 0) {
      normalize(texs[i]);
    }
  }

  inccs.resize(size);
  for (int i = 0; i < size; ++i) {
    inccs[i].resize(size);
  }

  for (int i = 0; i < size; ++i) {
    inccs[i][i] = 0.0f;
    for (int j = i+1; j < size; ++j) {
      if (!texs[i].empty() && !texs[j].empty()) {
        if (robust == 0) {
          inccs[j][i] = inccs[i][j] = 1.0f - dot(texs[i], texs[j]);
        } else {
          inccs[j][i] = inccs[i][j] = robustincc(1.0f - dot(texs[i], texs[j]));
        }
      } else {
        inccs[j][i] = inccs[i][j] = 2.0f;
      }
    }
  }
}

int PMVS3::COptim::grabSafe(const int index, const int size, const Vec3f& center, const Vec3f& dx, const Vec3f& dy, const int level) const {
  const int margin = size / 2;

  const Vec3f tl = center - dx * margin - dy * margin;
  const Vec3f tr = center + dx * margin - dy * margin;

  const Vec3f bl = center - dx * margin + dy * margin;
  const Vec3f br = center + dx * margin + dy * margin;

  const float minx = std::min(tl[0], std::min(tr[0], std::min(bl[0], br[0])));
  const float maxx = std::max(tl[0], std::max(tr[0], std::max(bl[0], br[0])));
  const float miny = std::min(tl[1], std::min(tr[1], std::min(bl[1], br[1])));
  const float maxy = std::max(tl[1], std::max(tr[1], std::max(bl[1], br[1])));

  // 1 should be enough
  const int margin2 = 3;
  // ??? may need to change if we change interpolation method
  if (minx < margin2 || _fm._pss.getWidth(index, level) - 1 - margin2 <= maxx || miny < margin2 || _fm._pss.getHeight(index, level) - 1 - margin2 <= maxy) {
    return 0;
  }

  return 1;
}

// My own optimisaton
float MyPow2(int x) {
	const float answers[] = {0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
	return answers[x + 4];
}

static float Log2 = log(2.0f);

int PMVS3::COptim::grabTex(const Vec4f& coord, const Vec4f& pxaxis, const Vec4f& pyaxis, const Vec4f& pzaxis, const int index, const int size, std::vector<float>& tex) const {
  tex.clear();

  Vec4f ray = _fm._pss._photos[index].OpticalCenter() - coord;
  unitize(ray);
  const float weight = std::max(0.0f, ray * pzaxis);

  //???????
  if (weight < cos(_fm._angleThreshold1)) return 1;

  const int margin = size / 2;

  Vec3f center = _fm._pss.project(index, coord, _fm._level);
  Vec3f dx = _fm._pss.project(index, coord + pxaxis, _fm._level) - center;
  Vec3f dy = _fm._pss.project(index, coord + pyaxis, _fm._level) - center;

  const float ratio = (norm(dx) + norm(dy)) / 2.0f;
  int leveldif = (int)floor(log(ratio) / Log2 + 0.5f);

  // Upper limit is 2
  leveldif = std::max(-_fm._level, std::min(2, leveldif));

  const float scale = MyPow2(leveldif);
  const int newlevel = _fm._level + leveldif;

  center /= scale;
  dx /= scale;
  dy /= scale;

  if (grabSafe(index, size, center, dx, dy, newlevel) == 0) return 1;

  Vec3f left = center - dx * margin - dy * margin;

  tex.resize(3 * size * size);
  float* texp = &tex[0] - 1;
  for (int y = 0; y < size; ++y) {
    Vec3f vftmp = left;
    left += dy;
    for (int x = 0; x < size; ++x) {
      Vec3f color = _fm._pss.getColor(index, vftmp[0], vftmp[1], newlevel);
      *(++texp) = color[0];
      *(++texp) = color[1];
      *(++texp) = color[2];
      vftmp += dx;
    }
  }

  return 0;
}

double PMVS3::COptim::computeINCC(const Vec4f& coord, const Vec4f& normal, const std::vector<int>& indexes, const int id, const int robust) {
  if ((int)indexes.size() < 2) return 2.0;

  const int index = indexes[0];
  Vec4f pxaxis, pyaxis;
  getPAxes(index, coord, normal, pxaxis, pyaxis);

  return computeINCC(coord, normal, indexes, pxaxis, pyaxis, id, robust);
}

double PMVS3::COptim::computeINCC(const Vec4f& coord, const Vec4f& normal, const std::vector<int>& indexes, const Vec4f& pxaxis, const Vec4f& pyaxis, const int id, const int robust) {
  if ((int)indexes.size() < 2) return 2.0;

  const int size = std::min(_fm._tau, (int)indexes.size());
  std::vector<std::vector<float> >& texs = _texsT[id];

  for (int i = 0; i < size; ++i) {
    int flag;
    flag = grabTex(coord, pxaxis, pyaxis, normal, indexes[i], _fm._wsize, texs[i]);

    if (flag == 0) {
      normalize(texs[i]);
    }
  }

  if (texs[0].empty()) {
    return 2.0;
  }

  double score = 0.0;

  // pure pairwise of reference based
#ifdef PMVS_PAIRNCC
  float totalweight = 0.0;
  for (int i = 0; i < size; ++i) {
    for (int j = i+1; j < size; ++j) {
      if (!texs[i].empty() && !texs[j].empty()) {
        const float ftmp = _weightsT[id][i] * _weightsT[id][j];
        totalweight += ftmp;
        if (robust) {
          score += robustincc(1.0 - dot(texs[i], texs[j])) * ftmp;
        } else {
          score += (1.0 - dot(texs[i], texs[j])) * ftmp;
        }
      }
    }
  }

  if (totalweight == 0.0) {
    score = 2.0;
  } else {
    score /= totalweight;
  }
#else
  float totalweight = 0.0;
  for (int i = 1; i < size; ++i) {
    if (!texs[i].empty()) {
      totalweight += _weightsT[id][i];
      if (robust) {
        score += robustincc(1.0 - dot(texs[0], texs[i])) * _weightsT[id][i];
      } else {
        score += (1.0 - dot(texs[0], texs[i])) * _weightsT[id][i];
      }
    }
  }
  if (totalweight == 0.0) {
    score = 2.0;
  } else {
    score /= totalweight;
  }
#endif

  return score;
}

void PMVS3::COptim::lfunc(double* p, double* hx, int m, int n, void* adata) {
  int iflag;
  _one->func(n, m, p, hx, &iflag, adata);
}

void PMVS3::COptim::func(int m, int n, double* x, double* fvec, int* iflag, void* arg) {
  const int id = *((int*)arg);
  double xs[3] = {x[0], x[1], x[2]};

  for (int i = 0; i < m; ++i) {
    fvec[i] = 2.0;
  }

  Vec4f coord, normal;
  decode(coord, normal, xs, id);

  const int index = _indexesT[id][0];
  Vec4f pxaxis, pyaxis;
  getPAxes(index, coord, normal, pxaxis, pyaxis);

  const int size = std::min(_fm._tau, (int)_indexesT[id].size());

  for (int i = 0; i < size; ++i) {
    int flag;
    flag = grabTex(coord, pxaxis, pyaxis, normal, _indexesT[id][i], _fm._wsize, _texsT[id][i]);

    if (flag == 0) {
      normalize(_texsT[id][i]);
    }
  }

  int count = -1;
  for (int i = 0; i < size; ++i) {
    for (int j = i+1; j < size; ++j) {
      count++;
      if (_texsT[id][i].empty() || _texsT[id][j].empty()) continue;

      fvec[count] = robustincc(1.0 - dot(_texsT[id][i], _texsT[id][j]));
    }
  }
}

// Normalize only scale for each image
void PMVS3::COptim::normalize(std::vector<std::vector<float> >& texs, const int size) {
  // compute average rgb
  Vec3f ave;
  int denom = 0;

  std::vector<Vec3f> rgbs;
  rgbs.resize(size);
  for (int i = 0; i < size; ++i) {
    if (texs[i].empty()) continue;

    int count = 0;
    while (count < (int)texs[i].size()) {
      rgbs[i][0] += texs[i][count++];
      rgbs[i][1] += texs[i][count++];
      rgbs[i][2] += texs[i][count++];
    }
    rgbs[i] /= (int)texs[i].size() / 3;

    ave += rgbs[i];
    ++denom;
  }

  // overall average
  if (denom == 0) return;

  ave /= denom;

  // Scale all the colors
  for (int i = 0; i < size; ++i) {
    if (texs[i].empty()) continue;

    int count = 0;
    // compute scale
    Vec3f scale;
    for (int j = 0; j < 3; ++j) {
      if (rgbs[i][j] != 0.0f) {
        scale[j] = ave[j] / rgbs[i][j];
      }
    }

    while (count < (int)texs[i].size()) {
      texs[i][count++] *= scale[0];
      texs[i][count++] *= scale[1];
      texs[i][count++] *= scale[2];
    }
  }
}

void PMVS3::COptim::normalize(std::vector<float>& tex) {
  const int size = (int)tex.size();
  const int size3 = size / 3;
  Vec3f ave;

  float* texp = &tex[0] - 1;
  for (int i = 0; i < size3; ++i) {
    ave[0] += *(++texp);
    ave[1] += *(++texp);
    ave[2] += *(++texp);
  }

  ave /= size3;

  float ave2 = 0.0;
  texp = &tex[0] - 1;
  for (int i = 0; i < size3; ++i) {
    const float f0 = ave[0] - *(++texp);
    const float f1 = ave[1] - *(++texp);
    const float f2 = ave[2] - *(++texp);

    ave2 += f0 * f0 + f1 * f1 + f2 * f2;
  }

  ave2 = sqrt(ave2 / size);

  if (ave2 == 0.0f) {
    ave2 = 1.0f;
  }

  texp = &tex[0] - 1;
  for (int i = 0; i < size3; ++i) {
    *(++texp) -= ave[0];    *texp /= ave2;
    *(++texp) -= ave[1];    *texp /= ave2;
    *(++texp) -= ave[2];    *texp /= ave2;
  }
}

float PMVS3::COptim::dot(const std::vector<float>& tex0, const std::vector<float>& tex1) const{
#ifndef PMVS_WNCC
  // Pierre Moulon (use classic access to array, windows STL do not like begin()-1)
  const int size = (int)tex0.size();
  float ans = 0.0f;
  for (int i = 0; i < size; ++i) {
    ans += tex0[i] * tex1[i];
  }
  return ans / size;
#else
  const int size = (int)tex0.size();
  std::vector<float>::const_iterator i0 = tex0.begin();
  std::vector<float>::const_iterator i1 = tex1.begin();
  float ans = 0.0f;
  for (int i = 0; i < size; ++i, ++i0, ++i1) {
    ans += (*i0) * (*i1) * _template[i];
  }
  return ans;
#endif
}

float PMVS3::COptim::ssd(const std::vector<float>& tex0, const std::vector<float>& tex1) const{
  const float scale = 0.01;

#ifndef PMVS_WNCC
  // Pierre Moulon (use classic access to array, windows STL do not like begin()-1)
  const int size = (int)tex0.size();
  float ans = 0.0f;
  for(int i = 0; i < size; ++i) {
    ans += fabs( tex0[i] - tex1[i] );
  }

  return scale * ans / size;
#else
  const int size = (int)tex0.size();
  std::vector<float>::const_iterator i0 = tex0.begin();
  std::vector<float>::const_iterator i1 = tex1.begin();
  float ans = 0.0f;
  for (int i = 0; i < size; ++i, ++i0, ++i1) {
    const float ftmp = fabs((*i0) - (*i1));
    //ans += (*i0) * (*i1) * _template[i];
    ans += ftmp * _template[i];
  }
  return scale * ans;
#endif
}

float PMVS3::COptim::getUnit(const int index, const Vec4f& coord) const {
  const float fz = norm(coord - _fm._pss._photos[index].OpticalCenter());
  const float ftmp = _ipscales[index];
  if (ftmp == 0.0) {
    return 1.0;
  }

  return 2.0 * fz * (0x0001 << _fm._level) / ftmp;
}

// get x and y axis to collect textures given reference image and normal
void PMVS3::COptim::getPAxes(const int index, const Vec4f& coord, const Vec4f& normal, Vec4f& pxaxis, Vec4f& pyaxis) const{
  // yasu changed here for fpmvs
  const float pscale = getUnit(index, coord);

  Vec3f normal3(normal[0], normal[1], normal[2]);
  Vec3f yaxis3 = cross(normal3, _xaxes[index]);
  unitize(yaxis3);
  Vec3f xaxis3 = cross(yaxis3, normal3);
  pxaxis[0] = xaxis3[0];  pxaxis[1] = xaxis3[1];  pxaxis[2] = xaxis3[2];  pxaxis[3] = 0.0;
  pyaxis[0] = yaxis3[0];  pyaxis[1] = yaxis3[1];  pyaxis[2] = yaxis3[2];  pyaxis[3] = 0.0;

  pxaxis *= pscale;
  pyaxis *= pscale;
  const float xdis = norm(_fm._pss.project(index, coord + pxaxis, _fm._level) - _fm._pss.project(index, coord, _fm._level));
  const float ydis = norm(_fm._pss.project(index, coord + pyaxis, _fm._level) - _fm._pss.project(index, coord, _fm._level));
  pxaxis /= xdis;
  pyaxis /= ydis;
}

void PMVS3::COptim::setWeightsT(const Patch::CPatch& patch, const int id) {
  computeUnits(patch, _weightsT[id]);
  for (int i = 1; i < (int)_weightsT[id].size(); ++i) {
    _weightsT[id][i] = std::min(1.0f, _weightsT[id][0] / _weightsT[id][i]);
  }
  _weightsT[id][0] = 1.0f;
}
