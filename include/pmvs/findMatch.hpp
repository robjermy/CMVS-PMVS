#pragma once

#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <queue>
#include "patch.hpp"

#include "image/photoSetS.hpp"
#include "patchOrganizerS.hpp"
#include "seed.hpp"
#include "expand.hpp"
#include "filter.hpp"
#include "optim.hpp"
#include "option.hpp"

#include <shared_mutex>
#include <thread>

namespace PMVS3 {
  class CFindMatch {
  public:
    CFindMatch();
    virtual ~CFindMatch();

    void init(const PMVS3::SOption& option);
    void run(void);
    void write(const std::string prefix, bool bExportPLY, bool bExportPatch, bool bExportPSet);

    int insideBimages(const Vec4f& coord) const;

    int isNeighborRadius(
      const Patch::CPatch& lhs,
      const Patch::CPatch& rhs,
      const float hunit,
      const float neighborThreshold,
      const float radius
    ) const;

    int isNeighbor(const Patch::CPatch& lhs, const Patch::CPatch& rhs, const float hunit, const float neighborThreshold) const;
    int isNeighbor(const Patch::CPatch& lhs, const Patch::CPatch& rhs, const float neighborThreshold) const;

    // Const getters
    const int NumImages() const { return _num; }
    const int NumTargetImages() const { return _tnum; }
    const int Count() const { return _count; }
    const std::list<int> Jobs() const { return _jobs; }
    const int JUnit() const { return _junit; }

    // Reference getters
    int& Count() { return _count; }
    std::list<int>& Jobs() { return _jobs; }

    // Thread locking and unlocking
    void Lock() { _lock.lock(); }
    void Unlock() { _lock.unlock(); }

    void LockImage(const int index) { _imageLocks[index]->lock(); }
    void LockSharedImage(const int index) { _imageLocks[index]->lock_shared(); }
    void UnlockImage(const int index) { _imageLocks[index]->unlock(); }
    void UnlockSharedImage(const int index) { _imageLocks[index]->unlock_shared(); }

    void LockCount(const int index) { _countLocks[index]->lock(); }
    void LockSharedCount(const int index) { _countLocks[index]->lock_shared(); }
    void UnlockCount(const int index) { _countLocks[index]->unlock(); }
    void UnlockSharedCount(const int index) { _countLocks[index]->unlock_shared(); }

    //----------------------------------------------------------------------
    // target images
    std::vector<int> _timages;
    // other images where patches are not computed
    std::vector<int> _oimages;
    // total images
    std::vector<int> _images;

    // prefix
    std::string _prefix;
    // level
    int _level;
    // cellsize
    int _csize;
    // nccThreshold
    float _nccThreshold;
    // windows size
    int _wsize;
    // mininum image num threshold
    int _minImageNumThreshold;
    // use edge detection or not
    float _setEdge;
    // bounding images
    std::vector<int> _bindexes;
    // visdata from SfM. _num x _num matrix
    std::vector<std::vector<int> > _visdata;
    // an array of relavant images
    std::vector<std::vector<int> > _visdata2;
    // sequence Threshold
    int _sequenceThreshold;
    // CPU
    int _CPU;
    // Threshold on filterQuad
    float _quadThreshold;

    // Maximum number of images used in the optimization
    int _tau;

    // If patches are dense or not, that is, if we use check(patch) after patch optimization
    int _depth;

    //----------------------------------------------------------------------
    // Thresholds
    //----------------------------------------------------------------------
    // For first feature matching. Images within this angle are used in
    // matching.
    float _angleThreshold0;
    // tigher angle
    float _angleThreshold1;

    // Number of success generation from each seed point
    int _countThreshold0;
    // Number of counts, expansion can be tried
    int _countThreshold1;

    // Number of trials for each cell in seed
    int _countThreshold2;

    // Parameter for isNeighbor in findemptyblocks
    float _neighborThreshold;
    // Parameter for isNeighbor in filterOutside
    float _neighborThreshold1;
    // Parameter for filterNeighbor
    float _neighborThreshold2;

    // ncc threshold before optim
    float _nccThresholdBefore;
    // Maximum angle of images must be at least as large as this
    float _maxAngleThreshold;

    // visibility consistency threshold
    float _visibleThreshold;
    float _visibleThresholdLoose;

    // Epipolar distance in seed generation
    float _epThreshold;

    //----------------------------------------------------------------------
    // Core members
    //----------------------------------------------------------------------
    // Images
    Image::CPhotoSetS _pss;
    // Patch organizer
    CPatchOrganizerS _pos;

    // Seed generator
    CSeed _seed;
    // Patch expansion
    CExpand _expand;
  public:
    // Patch filter
    CFilter _filter;
    // Patch optimizer
    COptim _optim;

    int _debug;

  protected: // variables
    // num of target images
    int _tnum;
    // num of total images
    int _num;
    // count
    int _count;


    //----------------------------------------------------------------------
    // For threads related
    //----------------------------------------------------------------------
    // General lock
    std::mutex _lock;
    // For each image
    std::vector<std::shared_mutex*> _imageLocks;
    std::vector<std::shared_mutex*> _countLocks;
    // jobs
    std::list<int> _jobs;
    // job unit
    int _junit;

  protected: // methods
    void init(void);
    void initTargets(void);
    void updateThreshold(void);
    void initImages(void);
  };
}
