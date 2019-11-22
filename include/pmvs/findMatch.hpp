#pragma once

#include <fstream>
#include <iostream>
#include <list>
#include <queue>
#include <shared_mutex>
#include <string>
#include <thread>
#include <vector>

#include "expand.hpp"
#include "filter.hpp"
#include "image/photoSetS.hpp"
#include "optim.hpp"
#include "option.hpp"
#include "patch.hpp"
#include "patchOrganizerS.hpp"
#include "seed.hpp"


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
    const std::list<int>& Jobs() const { return _jobs; }
    const int JUnit() const { return _junit; }
    const int CPU() const { return _CPU; }
    const Image::CPhotoSetS& PhotoSets() const { return _pss; }
    const CPatchOrganizerS& PatchOrganizer() const { return _pos; }
    const std::vector<int>& Images() const { return _images; }
    const std::string& Prefix() const { return _prefix; }
    const int Level() const { return _level; }
    const int CSize() const { return _csize; }
    const int WSize() const { return _wsize; }
    const int NCCThreshold() const { return _nccThreshold; }
    const int NCCThresholdBefore() const { return _nccThresholdBefore; }
    const int MinImageNumThreshold() const { return _minImageNumThreshold; }
    const int Tau() const { return _tau; }
    const int Depth() const { return _depth; }
    const int SequenceThreshold() const { return _sequenceThreshold; }
    const float QuadThreshold() const { return _quadThreshold; }
    const float AngleThreshold0() const { return _angleThreshold0; }
    const float AngleThreshold1() const { return _angleThreshold1; }
    const int CountThreshold0() const { return _countThreshold0; }
    const int CountThreshold1() const { return _countThreshold1; }
    const int CountThreshold2() const { return _countThreshold2; }
    const float NeighborThreshold0() const { return _neighborThreshold0; }
    const float NeighborThreshold1() const { return _neighborThreshold1; }
    const float NeighborThreshold2() const { return _neighborThreshold2; }
    const float MaxAngleThreshold() const { return _maxAngleThreshold; }
    const float EpipolarThreshold() const { return _epThreshold; }

    // Reference getters
    CExpand& Expand() { return _expand; }
    CFilter& Filter() { return _filter; }
    COptim& Optimizer() { return _optim; }

    int& Count() { return _count; }
    Image::CPhotoSetS& PhotoSets() { return _pss; }
    CPatchOrganizerS& PatchOrganizer() { return _pos; }
    std::list<int>& Jobs() { return _jobs; }
    int& Debug() { return _debug; }

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

    // bounding images
    std::vector<int> _bindexes;
    // visdata from SfM. _num x _num matrix
    std::vector<std::vector<int> > _visdata;
    // an array of relavant images
    std::vector<std::vector<int> > _visdata2;

  public:


  protected: // variables
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
    // Patch filter
    CFilter _filter;
    // Patch optimizer
    COptim _optim;


    // target images
    std::vector<int> _timages;
    // other images where patches are not computed
    std::vector<int> _oimages;
    // total images
    std::vector<int> _images;

    // num of target images
    int _tnum;
    // num of total images
    int _num;
    // count
    int _count;
    // CPU
    int _CPU;
    // prefix
    std::string _prefix;
    // level
    int _level;
    // cellsize
    int _csize;
    // windows size
    int _wsize;
    // mininum image num threshold
    int _minImageNumThreshold;
    // use edge detection or not
    float _setEdge;
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
    // sequence Threshold
    int _sequenceThreshold;
    // Threshold on filterQuad
    float _quadThreshold;

    // ncc threshold before optim
    float _nccThresholdBefore;
    // nccThreshold
    float _nccThreshold;

    // Number of success generation from each seed point
    int _countThreshold0;
    // Number of counts, expansion can be tried
    int _countThreshold1;
    // Number of trials for each cell in seed
    int _countThreshold2;

    // Parameter for isNeighbor in findemptyblocks
    float _neighborThreshold0;
    // Parameter for isNeighbor in filterOutside
    float _neighborThreshold1;
    // Parameter for filterNeighbor
    float _neighborThreshold2;
    
    // Maximum angle of images must be at least as large as this
    float _maxAngleThreshold;

    // visibility consistency threshold
    float _visibleThreshold;
    float _visibleThresholdLoose;

    // Epipolar distance in seed generation
    float _epThreshold;

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

    int _debug;

  protected: // methods
    void init(void);
    void initTargets(void);
    void updateThreshold(void);
    void initImages(void);
  };
}
