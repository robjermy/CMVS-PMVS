#include <iostream>
#include <string>
#include "../base/pmvs/findMatch.hpp"
#include "../base/pmvs/option.hpp"

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " prefix option_file [Optional export]" << std::endl
         << std::endl
         << "--------------------------------------------------" << std::endl
         << "level       1    csize    2" << std::endl
         << "threshold   0.7  wsize    7" << std::endl
         << "minImageNum 3    CPU      4" << std::endl
         << "useVisData  0    sequence -1" << std::endl
         << "quad        2.5  maxAngle 10.0" << std::endl
         << "--------------------------------------------------" << std::endl
         << "2 ways to specify targetting images" << std::endl
         << "timages  5  1 3 5 7 9 (enumeration)" << std::endl
         << "        -1  0 24 (range specification)" << std::endl
         << "--------------------------------------------------" << std::endl
         << "4 ways to specify other images" << std::endl
         << "oimages  5  0 2 4 6 8 (enumeration)" << std::endl
         << "        -1  24 48 (range specification)" << std::endl
         << std::endl
         << "[Optional export] PATCH PSET" << std::endl
         << " i.e export patch and pset: prefix option_file PATCH PSET"
         << " i.e export patch only: prefix option_file PATCH" << std::endl;
    exit (1);
  }

  for(int i=0; i < argc; ++i)
  {
	  std::cout << std::endl << argv[i];
  }
  std::cout << std::endl;

  PMVS3::SOption option;
  option.init(argv[1], argv[2]);

  PMVS3::CFindMatch findMatch;
  findMatch.init(option);
  findMatch.run();

  bool bExportPLY = true;
  bool bExportPatch = false;
  bool bExportPSet = false;

  for (int i=3; i < argc; ++i) {
    std::string option(argv[i]);
    if (option == "PATCH") {
      bExportPatch = true;
    }

    if (option == "PSET") {
      bExportPSet = true;
    }
  }

  char buffer[1024];
  sprintf(buffer, "%smodels/%s", argv[1], argv[2]);
  findMatch.write(buffer, bExportPLY, bExportPatch, bExportPSet);
}
