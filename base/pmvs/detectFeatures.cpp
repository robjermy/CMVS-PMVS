#include <iostream>
#include <fstream>
#include <thread>
#include "../image/image.hpp"
#include "detectFeatures.hpp"
#include "harris.hpp"
#include "dog.hpp"
#include "point.hpp"

PMVS3::CDetectFeatures::CDetectFeatures() {}

PMVS3::CDetectFeatures::~CDetectFeatures() {}

void PMVS3::CDetectFeatures::run(
  const Image::CPhotoSetS& pss,
  const int num,
  const int csize,
  const int level,
  const int CPU
) {
  _ppss = &pss;
  _csize = csize;
  _level = level;
  _CPU = CPU;

  _points.clear();
  _points.resize(num);

  for (int index = 0; index < num; ++index) {
    _jobs.push_back(index);
  }

  std::vector<std::thread> threads(_CPU);
  for (int i = 0; i < _CPU; ++i) {
    threads[i] = std::thread(&CDetectFeatures::runThread, this);
  }

  for (int i = 0; i < _CPU; ++i) {
    if (threads[i].joinable()) {
      threads[i].join();
    }
  }

  std::cerr << "done" << std::endl;
}

void PMVS3::CDetectFeatures::runThread(void) {
  while (1) {
    int index = -1;
    _lock.lock();
    {
      if (!_jobs.empty()) {
        index = _jobs.front();
        _jobs.pop_front();
      }
    }
    _lock.unlock();
    if (index == -1) break;

    const int image = _ppss->_images[index];
    std::cerr << image << ' ' << std::flush;

    //?????????????  May need file lock, because targetting images
    //should not overlap among multiple processors.
    char buffer[1024];
    sprintf(buffer, "%smodels/%08d.affin%d", _ppss->_prefix.c_str(), image, _level);
    std::ifstream ifstr;
    ifstr.open(buffer);
    if (ifstr.is_open()) {
      ifstr.close();
      continue;
    }
    ifstr.close();

    //----------------------------------------------------------------------
    // parameters
    // for harris
    const float sigma = 4.0f;
    // for DoG
    const float firstScale = 1.0f;
    const float lastScale = 3.0f;

    //----------------------------------------------------------------------
    // Harris
    {
      CHarris harris;
      std::multiset<CPoint> result;
      harris.run(
        _ppss->Photo(index).getImage(_level),
        _ppss->Photo(index).CImage::getMask(_level),
        _ppss->Photo(index).CImage::getEdge(_level),
        _ppss->Photo(index).getWidth(_level),
        _ppss->Photo(index).getHeight(_level), _csize, sigma, result
      );

      std::multiset<CPoint>::reverse_iterator rbegin = result.rbegin();
      while (rbegin != result.rend()) {
        _points[index].push_back(*rbegin);
        rbegin++;
      }
    }

    //----------------------------------------------------------------------
    // DoG
    {
      CDifferenceOfGaussians dog;
      std::multiset<CPoint> result;
      dog.run(
        _ppss->Photo(index).getImage(_level),
        _ppss->Photo(index).CImage::getMask(_level),
        _ppss->Photo(index).CImage::getEdge(_level),
        _ppss->Photo(index).getWidth(_level),
        _ppss->Photo(index).getHeight(_level),
        _csize, firstScale, lastScale, result
      );

      std::multiset<CPoint>::reverse_iterator rbegin = result.rbegin();
      while (rbegin != result.rend()) {
        _points[index].push_back(*rbegin);
        rbegin++;
      }
    }
  }
}
