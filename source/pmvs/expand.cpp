#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#include <thread>
#include "time.h"

#include "pmvs/expand.hpp"
#include "pmvs/findMatch.hpp"


PMVS3::CExpand::CExpand(CFindMatch& findMatch) : _fm(findMatch) {}

void PMVS3::CExpand::init() {}

void PMVS3::CExpand::run() {
  _fm.Count() = 0;
  _fm.Jobs().clear();
  _ecounts.resize(_fm._CPU);
  _fcounts0.resize(_fm._CPU);
  _fcounts1.resize(_fm._CPU);
  _pcounts.resize(_fm._CPU);

  fill(_ecounts.begin(), _ecounts.end(), 0);
  fill(_fcounts0.begin(), _fcounts0.end(), 0);
  fill(_fcounts1.begin(), _fcounts1.end(), 0);
  fill(_pcounts.begin(), _pcounts.end(), 0);

  time_t starttime = time(NULL);

  _fm._pos.clearCounts();
  _fm._pos.clearFlags();

  if (!_queue.empty()) {
    std::cerr << "Queue is not empty in expand" << std::endl;
    exit (1);
  }
  // set queue
  _fm._pos.collectPatches(_queue);

  std::cerr << "Expanding patches..." << std::flush;
  std::vector<std::thread> threads(_fm._CPU);
  for (int c = 0; c < _fm._CPU; ++c) {
    threads[c] = std::thread(&CExpand::expandThread, this);
  }

  for (int c = 0; c < _fm._CPU; ++c) {
    // thrd_join(threads[c], NULL);
    if (threads[c].joinable()) {
      threads[c].join();
    }
  }

  std::cerr << std::endl << "---- EXPANSION: " << (time(NULL) - starttime) << " secs ----" << std::endl;

  const int trial = accumulate(_ecounts.begin(), _ecounts.end(), 0);
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

void PMVS3::CExpand::expandThread(void) {
  _fm.Lock();
  const int id = _fm.Count()++;
  _fm.Unlock();

  while (1) {
    Patch::PPatch ppatch;
    int empty = 0;
    _fm.Lock();
    if (_queue.empty()) {
      empty = 1;
    } else {
      ppatch = _queue.top();
      _queue.pop();
    }
    _fm.Unlock();

    if (empty) break;

    // For each direction;
    std::vector<std::vector<Vec4f> > canCoords;
    findEmptyBlocks(ppatch, canCoords);

    for (int i = 0; i < (int)canCoords.size(); ++i) {
      for (int j = 0; j < (int)canCoords[i].size(); ++j) {
        const int flag = expandSub(ppatch, id, canCoords[i][j]);
        // fail
        if (flag) {
          ppatch->_dflag |= (0x0001) << i;
        }
      }
    }
  }
}

void PMVS3::CExpand::findEmptyBlocks(const Patch::PPatch& ppatch, std::vector<std::vector<Vec4f> >& canCoords) {
  // dnum must be at most 8, because _dflag is char
  const int dnum = 6;
  const Patch::CPatch& patch = *ppatch;

  // Empty six directions
  Vec4f xdir, ydir;
  ortho(ppatch->_normal, xdir, ydir);

  // -1: not empty
  // pos: number of free _pgrids
  //
  // Check if each direction satisfies both of the following two constraints.
  // a. No neighbor
  // b. At least minImageNumThreshold _pgrids without any patches and few _counts
  std::vector<float> fill;
  fill.resize(dnum);
  std::fill(fill.begin(), fill.end(), 0.0f);

  //----------------------------------------------------------------------
  // We look at the effective resolution of each image at the patch.
  // We only use images with good effective resolution to determine
  // empty blocks, because lwo-resolution images can easily satisfy
  // the first condition (neighbors), and no expansion will occur.
  // ----------------------------------------------------------------------
  // Minimum number of images required to obtain high res results, and
  // explor empty blocks.
  const float radius = computeRadius(patch);
  const float radiuslow = radius / 6.0f;//2.0f;
  const float radiushigh = radius * 2.5f;//2.0f;//1.5f;

  std::vector<Patch::PPatch> neighbors;
  _fm._pos.findNeighbors(patch, neighbors, 1, 4.0f);//3.0f);

  auto bpatch = neighbors.begin();
  auto epatch = neighbors.end();
  while (bpatch != epatch) {
    const Vec4f diff = (*bpatch)->_coord - ppatch->_coord;
    Vec2f f2(diff * xdir, diff * ydir);
    const float len = norm(f2);
    if (len < radiuslow || radiushigh < len) {
      ++bpatch;
      continue;
    }

    f2 /= len;

    float angle = atan2(f2[1], f2[0]);
    if (angle < 0.0) {
      angle += 2 * M_PI;
    }

    const float findex = angle / (2 * M_PI / dnum);
    const int lindex = (int)floor(findex);
    const int hindex = lindex + 1;

    fill[lindex % dnum] += hindex - findex;
    fill[hindex % dnum] += findex - lindex;
    ++bpatch;
  }

  canCoords.resize(dnum);
  for (int i = 0; i < dnum; ++i) {
    if (0.0f < fill[i]) continue;

    // If already failed, don't try, because we fail again.
    if (ppatch->_dflag & (0x0001 << i)) continue;

    const float angle = 2 * M_PI * i / dnum;
    Vec4f canCoord = ppatch->_coord + cos(angle) * radius * xdir + sin(angle) * radius * ydir;
    canCoords[i].push_back(canCoord);
  }
}

float PMVS3::CExpand::computeRadius(const Patch::CPatch& patch) {
  const int minnum = 2;
  std::vector<float> units;
  _fm._optim.computeUnits(patch, units);
  std::vector<float> vftmp = units;
#ifdef DEBUG
  if ((int)units.size() < minnum) {
    cerr << "units size less than minnum: " << (int)units.size() << ' ' << minnum << endl;
    cout << (int)patch._images.size() << endl;
    exit (1);
  }
#endif
  std::nth_element(vftmp.begin(), vftmp.begin() + minnum - 1, vftmp.end());
  // Threshold is the second smallest value with some margin
  // ??? critical
  return (*(vftmp.begin() + minnum - 1)) * _fm._csize;
}

int PMVS3::CExpand::expandSub(const Patch::PPatch& orgppatch, const int id, const Vec4f& canCoord) {
  // Choose the closest one
  Patch::CPatch patch;
  patch._coord = canCoord;
  patch._normal = orgppatch->_normal;
  patch._flag = 1;

  _fm._pos.setGridsImages(patch, orgppatch->_images);
  if (patch._images.empty()) return 1;

  //-----------------------------------------------------------------
  // Check bimages and mask. Then, initialize possible visible images
  if (_fm._pss.getMask(patch._coord, _fm._level) == 0 || _fm.insideBimages(patch._coord) == 0) return 1;

  // Check _counts and maybe _pgrids
  const int flag = checkCounts(patch);
  if (flag) return 1;

  // Check edge
  _fm._optim.removeImagesEdge(patch);
  if (patch._images.empty()) return 1;

  ++_ecounts[id];
  //-----------------------------------------------------------------
  // Preprocess
  if (_fm._optim.preProcess(patch, id, 0)) {
    ++_fcounts0[id];
    return 1;
  }

  //-----------------------------------------------------------------
  _fm._optim.refinePatch(patch, id, 100);

  //-----------------------------------------------------------------
  if (_fm._optim.postProcess(patch, id, 0)) {
    ++_fcounts1[id];
    return 1;
  }
  ++_pcounts[id];

  //-----------------------------------------------------------------
  // Finally
  Patch::PPatch ppatch(new Patch::CPatch(patch));

  //patch._images = orgppatch->_images;
  const int add = updateCounts(patch);

  _fm._pos.addPatch(ppatch);

  if (add) {
    _fm.Lock();
    _queue.push(ppatch);
    _fm.Unlock();
  }

  return 0;
}

int PMVS3::CExpand::checkCounts(Patch::CPatch& patch) {
  int full = 0;  int empty = 0;

  auto begin = patch._images.begin();
  auto end = patch._images.end();
  auto begin2 = patch._grids.begin();

  while (begin != end) {
    const int index = *begin;
    if (_fm.NumTargetImages() <= index) {
      ++begin;
      ++begin2;
      continue;
    }

    const int ix = (*begin2)[0];    const int iy = (*begin2)[1];
    if (ix < 0 || _fm._pos._gwidths[index] <= ix || iy < 0 || _fm._pos._gheights[index] <= iy) {
      ++begin;
      ++begin2;
      continue;
    }

    const int index2 = iy * _fm._pos._gwidths[index] + ix;

    int flag = 0;
	  _fm.LockSharedImage(index);;
    if (!_fm._pos._pgrids[index][index2].empty()) {
      flag = 1;
    }
	  _fm.UnlockSharedImage(index);

    if (flag) {
      ++full;
      ++begin;
      ++begin2;
      continue;
    }

    //mtx_lock(&_fm._countLocks[index]);
	  _fm.LockSharedCount(index);
    if (_fm._countThreshold1 <= _fm._pos._counts[index][index2]) {
      ++full;
    } else {
      ++empty;
    }
	  _fm.UnlockSharedCount(index);
    ++begin;
    ++begin2;
  }

  //First expansion is expensive and make the condition strict
  if (_fm._depth <= 1) {
    if (empty < _fm._minImageNumThreshold && full != 0) {
      return 1;
    } else {
      return 0;
    }
  }
  else {
    if (empty < _fm._minImageNumThreshold - 1 && full != 0) {
      return 1;
    } else {
      return 0;
    }
  }
}

int PMVS3::CExpand::updateCounts(const Patch::CPatch& patch) {
  // Use _images and _vimages. Loosen when to set add = 1
  int full = 0;  int empty = 0;
  {
    std::vector<int>::const_iterator begin = patch._images.begin();
    std::vector<int>::const_iterator end = patch._images.end();
    std::vector<Vec2i>::const_iterator begin2 = patch._grids.begin();

    while (begin != end) {
      const int index = *begin;
      if (_fm.NumTargetImages() <= index) {
        ++begin;
        ++begin2;
        continue;
      }

      const int ix = (*begin2)[0];
      const int iy = (*begin2)[1];
      if (ix < 0 || _fm._pos._gwidths[index] <= ix || iy < 0 || _fm._pos._gheights[index] <= iy) {
        ++begin;
        ++begin2;
        continue;
      }

      const int index2 = iy * _fm._pos._gwidths[index] + ix;

	    _fm.LockSharedCount(index);
      if (_fm._countThreshold1 <= _fm._pos._counts[index][index2]) {
        ++full;
      } else {
        ++empty;
      }
      ++_fm._pos._counts[index][index2];
	    _fm.UnlockSharedCount(index);
      ++begin;
      ++begin2;
    }
  }

  {
    std::vector<int>::const_iterator begin = patch._vimages.begin();
    std::vector<int>::const_iterator end = patch._vimages.end();
    std::vector<Vec2i>::const_iterator begin2 = patch._vgrids.begin();

    while (begin != end) {
      const int index = *begin;
#ifdef DEBUG
      if (_fm.NumTargetImages() <= index) {
        cerr << "Impossible in updateCounts" << endl;
        exit (1);
      }
#endif

      const int ix = (*begin2)[0];
      const int iy = (*begin2)[1];
      if (ix < 0 || _fm._pos._gwidths[index] <= ix || iy < 0 || _fm._pos._gheights[index] <= iy) {
        ++begin;
        ++begin2;
        continue;
      }

      const int index2 = iy * _fm._pos._gwidths[index] + ix;

	    _fm.LockSharedCount(index);
      if (_fm._countThreshold1 <= _fm._pos._counts[index][index2]) {
        ++full;
      } else {
        ++empty;
      }
      ++_fm._pos._counts[index][index2];
	    _fm.UnlockSharedCount(index);
      ++begin;
      ++begin2;
    }
  }

  if (empty != 0) {
    return 1;
  }

  return 0;
}
