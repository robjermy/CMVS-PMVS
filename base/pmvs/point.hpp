#ifndef PMVS3_POINT_H
#define PMVS3_POINT_H

#include "../numeric/vec4.hpp"
#include "../numeric/mat3.hpp"

namespace PMVS3 {
class CPoint {
 public:
  CPoint(void);
  virtual ~CPoint();
  
  Vec3f m_icoord;
  float m_response;

  // 0: Harris
  // 1: DoG
  int m_type;

  // tempporary variable, used to store original imageid in initial match
  int m_itmp;

  // 3D coordinate
  Vec4f m_coord;
  
  bool operator < (const CPoint& rhs) const {
    return m_response < rhs.m_response;
  }

  friend std::istream& operator >>(std::istream& istr, CPoint& rhs);
  friend std::ostream& operator <<(std::ostream& ostr, const CPoint& rhs);
};

bool SortCpoint(const CPoint& a, const CPoint& b);

std::istream& operator >>(std::istream& istr, CPoint& rhs);
std::ostream& operator <<(std::ostream& ostr, const CPoint& rhs);
};

#endif //PMVS3_POINT_H
