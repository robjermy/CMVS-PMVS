#include <string>

#include "numeric/vec4.hpp"
#include "pmvs/patch.hpp"

std::istream& Patch::operator>>(std::istream& istr, CPatch& rhs) {
  std::string header;
  int itmp;
  istr >> header >> rhs._coord >> rhs._normal >> rhs._ncc
       >> rhs._dscale >> rhs._ascale;

  if (header == "PATCHA") {
    int type;    Vec4f dir;
    istr >> type >> dir;
  }

  istr >> itmp;
  rhs._images.resize(itmp);
  for (int i = 0; i < itmp; ++i)
    istr >> rhs._images[i];

  istr >> itmp;
  rhs._vimages.resize(itmp);
  for (int i = 0; i < itmp; ++i)
    istr >> rhs._vimages[i];

  return istr;
}

std::ostream& Patch::operator<<(std::ostream& ostr, const CPatch& rhs) {
  ostr << "PATCHS" << std::endl
       << rhs._coord << std::endl
       << rhs._normal << std::endl
       << rhs._ncc << ' '
       << rhs._dscale << ' '
       << rhs._ascale << std::endl
       << (int)rhs._images.size() << std::endl;
  for (int i = 0; i < (int)rhs._images.size(); ++i)
    ostr << rhs._images[i] << ' ';
  ostr << std::endl;

  ostr << (int)rhs._vimages.size() << std::endl;
  for (int i = 0; i < (int)rhs._vimages.size(); ++i)
    ostr << rhs._vimages[i] << ' ';
  ostr << std::endl;

  return ostr;
}
