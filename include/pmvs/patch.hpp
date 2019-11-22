#pragma once

#include <vector>
#include <iostream>
#include <memory>
#include "numeric/vec4.hpp"

namespace Patch {
  class CPatch {
  public:
    CPatch(void) {
      _ncc = -1.0;
      _timages = 0;
      _fix = 0;
      // dflag is initialized only once. if failed in one direction, we
      // never try that.
      _dflag = 0;

      // All non-class member variables need to be initialized so that
      // they aren't just uninitialized memory.
      _flag = 0;
      _id = 0;
      _dscale = 0;
      _ascale = 0;
      _tmp = 0;
    }

    //----------------------------------------------------------------------
    // saved information
    // 3D coordinates of the center of the patch
    Vec4f _coord;
    // patch outward normal vector
    Vec4f _normal;

    // associated image ids. first image id is the reference one. images
    // can be non-targetting image.
    std::vector<int> _images;
    std::vector<TVec2<int> > _grids;

    // visible images. _vimages must be targetting images.
    std::vector<int> _vimages;
    std::vector<TVec2<int> > _vgrids;

    //----------------------------------------------------------------------
    inline float score(const float threshold) const{
      return std::max(0.0f, _ncc - threshold) * (int)_images.size();
    }
    inline float score2(const float threshold) const{
      return std::max(0.0f, _ncc - threshold) * _timages;
    }

    // average ncc
    float _ncc;
    // number of targetting images in _images
    int _timages;

    // flat for expansion
    // 0: not yet tested
    // 1: done
    int _flag;

    // for directional flag
    unsigned char _dflag;

    // fixed patch or not
    char _fix;

    // id number in _ppatches
    int _id;

    // scaling factor corresponding to one pixel difference
    float _dscale;
    float _ascale;

    float _tmp;
  };

  typedef std::shared_ptr<CPatch> PPatch;

  struct Spatchcmp {
    bool operator()(const PPatch& lhs, const PPatch& rhs) {
      return lhs.get() < rhs.get();
    }
  };

  std::istream& operator>>(std::istream& istr, Patch::CPatch& rhs);
  std::ostream& operator<<(std::ostream& ostr, const Patch::CPatch& rhs);
}
