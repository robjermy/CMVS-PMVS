#pragma once

#include "../numeric/vec4.hpp"
#include "../numeric/mat3.hpp"

namespace PMVS3 {
  class CPoint {
  public:
    CPoint(void);
    virtual ~CPoint();

    Vec3f _icoord;
    float _response;

    // 0: Harris
    // 1: DoG
    int _type;

    // tempporary variable, used to store original imageid in initial match
    int _itmp;

    // 3D coordinate
    Vec4f _coord;

    bool operator<(const CPoint& rhs) const {
      return _response < rhs._response;
    }

    friend std::istream& operator>>(std::istream& istr, CPoint& rhs);
    friend std::ostream& operator<<(std::ostream& ostr, const CPoint& rhs);
  };

  bool SortCpoint(const CPoint& a, const CPoint& b);

  std::istream& operator>>(std::istream& istr, CPoint& rhs);
  std::ostream& operator<<(std::ostream& ostr, const CPoint& rhs);
}
