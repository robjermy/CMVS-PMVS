#include "point.hpp"
#include <iostream>

using namespace PMVS3;
using namespace std;

CPoint::CPoint(void) {
  _response = -1.0;
  _type = -1;
}

CPoint::~CPoint() {
}

std::istream& PMVS3::operator >>(std::istream& istr, CPoint& rhs) {
  string header;
  char str[1024];
  istr >> str;
  header = string(str);
  istr >> rhs._icoord[0] >> rhs._icoord[1] >> rhs._response >> rhs._type;
  rhs._icoord[2] = 1.0f;
  return istr;
}

std::ostream& PMVS3::operator <<(std::ostream& ostr, const CPoint& rhs) {
  ostr << "POINT0" << endl
       << rhs._icoord[0] << ' ' << rhs._icoord[1] << ' ' << rhs._response << ' '
       << rhs._type;
  return ostr;
}

bool SortCpoint(const CPoint& a, const CPoint& b)
{
    return a._response < b._response;
}
