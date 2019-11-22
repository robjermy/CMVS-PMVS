#pragma once

/*
 * A main class to detect features
 */
#include <list>
#include <shared_mutex>
#include <string>

#include "image/photoSetS.hpp"
#include "point.hpp"

namespace Image {
  class CPhotoSetS;
};

namespace PMVS3 {
  class CDetectFeatures {
  public:
    CDetectFeatures(void);
    virtual ~CDetectFeatures();

    void run(
      const Image::CPhotoSetS& pss,
      const int num,
      const int csize,
      const int level,
      const int CPU = 1
    );

    std::vector<std::vector<CPoint>> _points;

  protected:
    const Image::CPhotoSetS* _ppss;
    int _csize;
    int _level;

    //----------------------------------------------------------------------
    // thread related
    //----------------------------------------------------------------------
    std::mutex _lock;
    int _CPU;

    std::list<int> _jobs;

    void runThread(void);
  };
}
