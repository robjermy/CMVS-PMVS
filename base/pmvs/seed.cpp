#include <time.h>
#include <numeric>
#include <ctime>
#include <time.h>
#include <random>
#include <thread>
#include "seed.hpp"
#include "findMatch.hpp"

PMVS3::CSeed::CSeed(CFindMatch& findMatch) : _fm(findMatch) {}

void PMVS3::CSeed::init(const std::vector<std::vector<CPoint> >& points) {
  _ppoints.clear();
  _ppoints.resize(_fm._num);

  for (int index = 0; index < _fm._num; ++index) {
    const int gheight = _fm._pos._gheights[index];
    const int gwidth = _fm._pos._gwidths[index];
    _ppoints[index].resize(gwidth * gheight);
  }

  readPoints(points);
}

void PMVS3::CSeed::readPoints(const std::vector<std::vector<CPoint> >& points) {
  for (int index = 0; index < _fm._num; ++index) {
    for (int i = 0; i < (int)points[index].size(); ++i) {
      PPoint ppoint(new CPoint(points[index][i]));
      ppoint->_itmp = index;
      const int ix = ((int)floor(ppoint->_icoord[0] + 0.5f)) / _fm._csize;
      const int iy = ((int)floor(ppoint->_icoord[1] + 0.5f)) / _fm._csize;
      const int index2 = iy * _fm._pos._gwidths[index] + ix;
      _ppoints[index][index2].push_back(ppoint);
    }
  }
}

// Arbitrary seed for deterministic pseudorandomness.
static const unsigned int RANDOM_SEED = 42;

void PMVS3::CSeed::run(void) {
  _fm._count = 0;
  _fm._jobs.clear();
  _scounts.resize(_fm._CPU);
  _fcounts0.resize(_fm._CPU);
  _fcounts1.resize(_fm._CPU);
  _pcounts.resize(_fm._CPU);

  std::fill(_scounts.begin(), _scounts.end(), 0);
  std::fill(_fcounts0.begin(), _fcounts0.end(), 0);
  std::fill(_fcounts1.begin(), _fcounts1.end(), 0);
  std::fill(_pcounts.begin(), _pcounts.end(), 0);

  std::vector<int> vitmp;
  for (int i = 0; i < _fm._tnum; ++i) {
    vitmp.push_back(i);
  }

  std::mt19937 gen(RANDOM_SEED);
  std::shuffle(vitmp.begin(), vitmp.end(), gen);
  _fm._jobs.insert(_fm._jobs.end(), vitmp.begin(), vitmp.end());

  std::cerr << "adding seeds " << std::endl;

  _fm._pos.clearCounts();

  // If there already exists a patch, don't use
  for (int index = 0; index < (int)_fm._tnum; ++index) {
    for (int j = 0; j < (int)_fm._pos._pgrids[index].size(); ++j) {
      if (!_fm._pos._pgrids[index][j].empty()) {
        _fm._pos._counts[index][j] = _fm._countThreshold2;
      }
    }
  }

  time_t tv;
  time(&tv);
  time_t curtime = tv;
  std::vector<std::thread> threads(_fm._CPU);
  for (int i = 0; i < _fm._CPU; ++i) {
    threads[i] = std::thread(&CSeed::initialMatchThread, this);
  }

  for (int i = 0; i < _fm._CPU; ++i) {
    if (threads[i].joinable()) {
      threads[i].join();
    }
  }

  //----------------------------------------------------------------------
  std::cerr << "done" << std::endl;
  time(&tv);
  std::cerr << "---- Initial: " << (tv - curtime)/CLOCKS_PER_SEC << " secs ----" << std::endl;

  const int trial = accumulate(_scounts.begin(), _scounts.end(), 0);
  const int fail0 = accumulate(_fcounts0.begin(), _fcounts0.end(), 0);
  const int fail1 = accumulate(_fcounts1.begin(), _fcounts1.end(), 0);
  const int pass = accumulate(_pcounts.begin(), _pcounts.end(), 0);
  std::cerr << "Total pass fail0 fail1 refinepatch: "
       << trial << ' ' << pass << ' '
       << fail0 << ' ' << fail1 << ' ' << pass + fail1 << std::endl;
  std::cerr << "Total pass fail0 fail1 refinepatch: "
       << 100 * trial / (float)trial << ' '
       << 100 * pass / (float)trial << ' '
       << 100 * fail0 / (float)trial << ' '
       << 100 * fail1 / (float)trial << ' '
       << 100 * (pass + fail1) / (float)trial << std::endl;
}

void PMVS3::CSeed::initialMatchThread(void) {
  _fm._lock.lock();
  const int id = _fm._count++;
  _fm._lock.unlock();

  while (1) {
    int index = -1;
    _fm._lock.lock();
    if (!_fm._jobs.empty()) {
      index = _fm._jobs.front();
      _fm._jobs.pop_front();
    }
    _fm._lock.unlock();

    if (index == -1) break;

    initialMatch(index, id);
  }
}

void PMVS3::CSeed::clear(void) {
  std::vector<std::vector<std::vector<PPoint>>>().swap(_ppoints);
}

void PMVS3::CSeed::initialMatch(const int index, const int id) {
  std::vector<int> indexes;
  _fm._optim.collectImages(index, indexes);

  if (_fm._tau < (int)indexes.size()) {
    indexes.resize(_fm._tau);
  }

  if (indexes.empty()) return;

  int totalcount = 0;
  //======================================================================
  // for each feature point, starting from the optical center, keep on
  // matching until we find candidateThreshold patches
  const int gheight = _fm._pos._gheights[index];
  const int gwidth = _fm._pos._gwidths[index];

  int index2 = -1;
  for (int y = 0; y < gheight; ++y) {
    for (int x = 0; x < gwidth; ++x) {
      ++index2;
      if (!canAdd(index, x, y)) continue;

      for (int p = 0; p < (int)_ppoints[index][index2].size(); ++p) {
        // collect features that satisfies epipolar geometry
        // constraints and sort them according to the differences of
        // distances between two cameras.
        std::vector<PPoint> vcp;
        collectCandidates(index, indexes, *_ppoints[index][index2][p], vcp);

        int count = 0;
        Patch::CPatch bestpatch;
        //======================================================================
        for (int i = 0; i < (int)vcp.size(); ++i) {
          Patch::CPatch patch;
          patch._coord = vcp[i]->_coord;
          patch._normal = _fm._pss.Photo(index).OpticalCenter() - patch._coord;

          unitize(patch._normal);
          patch._normal[3] = 0.0;
          patch._flag = 0;

          ++_fm._pos._counts[index][index2];
          const int ix = ((int)floor(vcp[i]->_icoord[0] + 0.5f)) / _fm._csize;
          const int iy = ((int)floor(vcp[i]->_icoord[1] + 0.5f)) / _fm._csize;
          const int index3 = iy * _fm._pos._gwidths[vcp[i]->_itmp] + ix;
          if (vcp[i]->_itmp < _fm._tnum) {
            ++_fm._pos._counts[vcp[i]->_itmp][index3];
          }

          const int flag = initialMatchSub(index, vcp[i]->_itmp, id, patch);
          if (flag == 0) {
            ++count;
            if (bestpatch.score(_fm._nccThreshold) < patch.score(_fm._nccThreshold)) {
              bestpatch = patch;
            }

            if (_fm._countThreshold0 <= count) break;
          }
        }

        if (count != 0) {
          Patch::PPatch ppatch(new Patch::CPatch(bestpatch));
          _fm._pos.addPatch(ppatch);
          ++totalcount;
          break;
        }
      }
    }
  }

  std::cerr << '(' << index << ',' << totalcount << ')' << std::flush;
}

void PMVS3::CSeed::collectCells(const int index0, const int index1, const CPoint& p0, std::vector<Vec2i>& cells) {
  Vec3 point(p0._icoord[0], p0._icoord[1], p0._icoord[2]);
#ifdef DEBUG
  if (p0._icoord[2] != 1.0f) {
    cerr << "Impossible in collectCells" << endl;    exit (1);
  }
#endif

  Mat3 F;
  Image::setF(_fm._pss.Photo(index0), _fm._pss.Photo(index1), F, _fm._level);
  const int gwidth = _fm._pos._gwidths[index1];
  const int gheight = _fm._pos._gheights[index1];

  Vec3 line = transpose(F) * point;
  if (line[0] == 0.0 && line[1] == 0.0) {
    std::cerr << "Point right on top of the epipole?"
         << index0 << ' ' << index1 << std::endl;
    return;
  }

  // vertical
  if (fabs(line[0]) > fabs(line[1])) {
    for (int y = 0; y < gheight; ++y) {
      const float fy = (y + 0.5) * _fm._csize - 0.5f;
      float fx = (- line[1] * fy - line[2]) / line[0];
      fx = std::max((float)(INT_MIN + 3.0f), std::min((float)(INT_MAX - 3.0f), fx));

      const int ix = ((int)floor(fx + 0.5f)) / _fm._csize;
      if (0 <= ix && ix < gwidth) {
        cells.push_back(TVec2<int>(ix, y));
      }

      if (0 <= ix - 1 && ix - 1 < gwidth) {
        cells.push_back(TVec2<int>(ix - 1, y));
      }

      if (0 <= ix + 1 && ix + 1 < gwidth) {
        cells.push_back(TVec2<int>(ix + 1, y));
      }
    }
  } else {
    for (int x = 0; x < gwidth; ++x) {
      const float fx = (x + 0.5) * _fm._csize - 0.5f;
      float fy = (- line[0] * fx - line[2]) / line[1];
      fy = std::max((float)(INT_MIN + 3.0f), std::min((float)(INT_MAX - 3.0f), fy));

      const int iy = ((int)floor(fy + 0.5f)) / _fm._csize;
      if (0 <= iy && iy < gheight) {
        cells.push_back(TVec2<int>(x, iy));
      }

      if (0 <= iy - 1 && iy - 1 < gheight) {
        cells.push_back(TVec2<int>(x, iy - 1));
      }

      if (0 <= iy + 1 && iy + 1 < gheight) {
        cells.push_back(TVec2<int>(x, iy + 1));
      }
    }
  }
}

// make sorted array of feature points in images, that satisfy the
// epipolar geometry coming from point in image
void PMVS3::CSeed::collectCandidates(const int index, const std::vector<int>& indexes, const CPoint& point, std::vector<PPoint>& vcp) {
  const Vec3 p0(point._icoord[0], point._icoord[1], 1.0);
  for (int i = 0; i < (int)indexes.size(); ++i) {
    const int indexid = indexes[i];

    std::vector<TVec2<int> > cells;
    collectCells(index, indexid, point, cells);
    Mat3 F;
    Image::setF(_fm._pss.Photo(index), _fm._pss.Photo(indexid), F, _fm._level);

    for (int i = 0; i < (int)cells.size(); ++i) {
      const int x = cells[i][0];      const int y = cells[i][1];
      if (!canAdd(indexid, x, y)) continue;

      const int index2 = y * _fm._pos._gwidths[indexid] + x;

      auto begin = _ppoints[indexid][index2].begin();
      auto end = _ppoints[indexid][index2].end();
      while (begin != end) {
        CPoint& rhs = **begin;
        // ? use type to reject candidates?
        if (point._type != rhs._type) {
          ++begin;
          continue;
        }

        const Vec3 p1(rhs._icoord[0], rhs._icoord[1], 1.0);
        if (_fm._epThreshold <= Image::computeEPD(F, p0, p1)) {
          ++begin;
          continue;
        }
        vcp.push_back(*begin);
        ++begin;
      }
    }
  }

  // set distances to _response
  std::vector<PPoint> vcptmp;
  for (int i = 0; i < (int)vcp.size(); ++i) {
    unproject(index, vcp[i]->_itmp, point, *vcp[i], vcp[i]->_coord);

    if (_fm._pss.Photo(index).ProjectionMatrix()[_fm._level][2] * vcp[i]->_coord <= 0.0) continue;
    if (_fm._pss.getMask(vcp[i]->_coord, _fm._level) == 0 || _fm.insideBimages(vcp[i]->_coord) == 0) continue;

    //??? from the closest
    vcp[i]->_response = fabs(norm(vcp[i]->_coord - _fm._pss.Photo(index).OpticalCenter()) - norm(vcp[i]->_coord - _fm._pss.Photo(vcp[i]->_itmp).OpticalCenter()));

    vcptmp.push_back(vcp[i]);
  }
  vcptmp.swap(vcp);
  sort(vcp.begin(), vcp.end());
}

int PMVS3::CSeed::canAdd(const int index, const int x, const int y) {
  if (!_fm._pss.getMask(index, _fm._csize * x, _fm._csize * y, _fm._level)) return 0;
  if (_fm._tnum <= index) return 1;

  const int index2 = y * _fm._pos._gwidths[index] + x;

  // Check if _pgrids already contains something
  if (!_fm._pos._pgrids[index][index2].empty()) return 0;

  //??? critical
  if (_fm._countThreshold2 <= _fm._pos._counts[index][index2]) return 0;

  return 1;
}

void PMVS3::CSeed::unproject(const int index0, const int index1, const CPoint& p0, const CPoint& p1, Vec4f& coord) const{
  Mat4 A;
  A[0][0] = _fm._pss.Photo(index0).ProjectionMatrix()[_fm._level][0][0] - p0._icoord[0] * _fm._pss.Photo(index0).ProjectionMatrix()[_fm._level][2][0];
  A[0][1] = _fm._pss.Photo(index0).ProjectionMatrix()[_fm._level][0][1] - p0._icoord[0] * _fm._pss.Photo(index0).ProjectionMatrix()[_fm._level][2][1];
  A[0][2] = _fm._pss.Photo(index0).ProjectionMatrix()[_fm._level][0][2] - p0._icoord[0] * _fm._pss.Photo(index0).ProjectionMatrix()[_fm._level][2][2];
  A[1][0] = _fm._pss.Photo(index0).ProjectionMatrix()[_fm._level][1][0] - p0._icoord[1] * _fm._pss.Photo(index0).ProjectionMatrix()[_fm._level][2][0];
  A[1][1] = _fm._pss.Photo(index0).ProjectionMatrix()[_fm._level][1][1] - p0._icoord[1] * _fm._pss.Photo(index0).ProjectionMatrix()[_fm._level][2][1];
  A[1][2] = _fm._pss.Photo(index0).ProjectionMatrix()[_fm._level][1][2] - p0._icoord[1] * _fm._pss.Photo(index0).ProjectionMatrix()[_fm._level][2][2];
  A[2][0] = _fm._pss.Photo(index1).ProjectionMatrix()[_fm._level][0][0] - p1._icoord[0] * _fm._pss.Photo(index1).ProjectionMatrix()[_fm._level][2][0];
  A[2][1] = _fm._pss.Photo(index1).ProjectionMatrix()[_fm._level][0][1] - p1._icoord[0] * _fm._pss.Photo(index1).ProjectionMatrix()[_fm._level][2][1];
  A[2][2] = _fm._pss.Photo(index1).ProjectionMatrix()[_fm._level][0][2] - p1._icoord[0] * _fm._pss.Photo(index1).ProjectionMatrix()[_fm._level][2][2];
  A[3][0] = _fm._pss.Photo(index1).ProjectionMatrix()[_fm._level][1][0] - p1._icoord[1] * _fm._pss.Photo(index1).ProjectionMatrix()[_fm._level][2][0];
  A[3][1] = _fm._pss.Photo(index1).ProjectionMatrix()[_fm._level][1][1] - p1._icoord[1] * _fm._pss.Photo(index1).ProjectionMatrix()[_fm._level][2][1];
  A[3][2] = _fm._pss.Photo(index1).ProjectionMatrix()[_fm._level][1][2] - p1._icoord[1] * _fm._pss.Photo(index1).ProjectionMatrix()[_fm._level][2][2];

  Vec4 b;
  b[0] = p0._icoord[0] * _fm._pss.Photo(index0).ProjectionMatrix()[_fm._level][2][3] - _fm._pss.Photo(index0).ProjectionMatrix()[_fm._level][0][3];
  b[1] = p0._icoord[1] * _fm._pss.Photo(index0).ProjectionMatrix()[_fm._level][2][3] - _fm._pss.Photo(index0).ProjectionMatrix()[_fm._level][1][3];
  b[2] = p1._icoord[0] * _fm._pss.Photo(index1).ProjectionMatrix()[_fm._level][2][3] - _fm._pss.Photo(index1).ProjectionMatrix()[_fm._level][0][3];
  b[3] = p1._icoord[1] * _fm._pss.Photo(index1).ProjectionMatrix()[_fm._level][2][3] - _fm._pss.Photo(index1).ProjectionMatrix()[_fm._level][1][3];

  Mat4 AT = transpose(A);
  Mat4 ATA = AT * A;
  Vec4 ATb = AT * b;

  Mat3 ATA3;
  for (int y = 0; y < 3; ++y) {
    for (int x = 0; x < 3; ++x) {
      ATA3[y][x] = ATA[y][x];
    }
  }

  Vec3 ATb3;
  for (int y = 0; y < 3; ++y) {
    ATb3[y] = ATb[y];
  }

  Mat3 iATA3;
  invert(iATA3, ATA3);
  Vec3 ans = iATA3 * ATb3;
  for (int y = 0; y < 3; ++y) {
    coord[y] = ans[y];
  }
  coord[3] = 1.0f;
}

// starting with (index, indexs), set visible images by looking at correlation.
int PMVS3::CSeed::initialMatchSub(const int index0, const int index1, const int id, Patch::CPatch& patch) {
  //----------------------------------------------------------------------
  patch._images.clear();
  patch._images.push_back(index0);
  patch._images.push_back(index1);

  ++_scounts[id];

  //----------------------------------------------------------------------
  // We know that patch._coord is inside bimages and inside mask
  if (_fm._optim.preProcess(patch, id, 1)) {
    ++_fcounts0[id];
    return 1;
  }

  //----------------------------------------------------------------------
  _fm._optim.refinePatch(patch, id, 100);

  //----------------------------------------------------------------------
  if (_fm._optim.postProcess(patch, id, 1)) {
    ++_fcounts1[id];
    return 1;
  }

  ++_pcounts[id];
  //----------------------------------------------------------------------
  return 0;
}
