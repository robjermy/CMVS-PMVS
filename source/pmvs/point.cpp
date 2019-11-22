#include "pmvs/point.hpp"
#include <iostream>

PMVS3::CPoint::CPoint(void) {
  _response = -1.0;
  _type = -1;
}

PMVS3::CPoint::~CPoint() {
}

std::istream& PMVS3::operator >>(std::istream& istr, CPoint& rhs) {
  std::string header;
  char str[1024];
  istr >> str;
  header = std::string(str);
  istr >> rhs._icoord[0] >> rhs._icoord[1] >> rhs._response >> rhs._type;
  rhs._icoord[2] = 1.0f;
  return istr;
}

std::ostream& PMVS3::operator <<(std::ostream& ostr, const CPoint& rhs) {
  ostr << "POINT0" << std::endl
       << rhs._icoord[0] << ' ' << rhs._icoord[1] << ' ' << rhs._response << ' '
       << rhs._type;
  return ostr;
}

bool SortCpoint(const PMVS3::CPoint& a, const PMVS3::CPoint& b) {
  return a._response < b._response;
}
