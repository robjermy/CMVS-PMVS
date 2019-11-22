#define _USE_MATH_DEFINES
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>

#include "pmvs/option.hpp"

PMVS3::SOption::SOption(void) {
  _level = 1;
  _csize = 2;
  _threshold = 0.7;
  _wsize = 7;
  _minImageNum = 3;
  _CPU = 4;
  _setEdge = 0.0f;
  _useBound = 0;
  _useVisData = 0;
  _sequence = -1;
  _tflag = -10;
  _oflag = -10;

  // Max angle must be at least this big
  _maxAngleThreshold = 10.0f * M_PI / 180.0f;
  // The smaller the tighter
  _quadThreshold = 2.5f;
}

void PMVS3::SOption::init(const std::string prefix, const std::string option) {
  _prefix = prefix;
  _option = option;
  std::ifstream ifstr;
  std::string optionfile = prefix + option;
  ifstr.open(optionfile.c_str());
  while (1) {
    std::string name;
    ifstr >> name;
    if (ifstr.eof()) break;

    if (name[0] == '#') {
      char buffer[1024];
      ifstr.putback('#');
      ifstr.getline(buffer, 1024);
      continue;
    }
    if (name == "level") {
      ifstr >> _level;
    } else if (name == "csize") {
      ifstr >> _csize;
    } else if (name == "threshold") {
      ifstr >> _threshold;
    } else if (name == "wsize") {
      ifstr >> _wsize;
    } else if (name == "minImageNum") {
      ifstr >> _minImageNum;
    } else if (name == "CPU") {
      ifstr >> _CPU;
    } else if (name == "setEdge") {
      ifstr >> _setEdge;
    } else if (name == "useBound") {
      ifstr >> _useBound;
    } else if (name == "useVisData") {
      ifstr >> _useVisData;
    } else if (name == "sequence") {
      ifstr >> _sequence;
    } else if (name == "timages") {
      ifstr >> _tflag;
      if (_tflag == -1) {
        int firstimage, lastimage;
        ifstr >> firstimage >> lastimage;
        for (int i = firstimage; i < lastimage; ++i) {
          _timages.push_back(i);
        }
      }
      else if (0 < _tflag) {
        for (int i = 0; i < _tflag; ++i) {
          int itmp;
          ifstr >> itmp;
          _timages.push_back(itmp);
        }
      } else {
        std::cerr << "tflag is not valid: " << _tflag << std::endl;   exit (1);
      }
    } else if (name == "oimages") {
      ifstr >> _oflag;
      if (_oflag == -1) {
        int firstimage, lastimage;
        ifstr >> firstimage >> lastimage;
        for (int i = firstimage; i < lastimage; ++i) {
          _oimages.push_back(i);
        }
      } else if (0 <= _oflag) {
        for (int i = 0; i < _oflag; ++i) {
          int itmp;
          ifstr >> itmp;
          _oimages.push_back(itmp);
        }
      } else if (_oflag != -2 && _oflag != -3) {
        std::cerr << "oflag is not valid: " << _oflag << std::endl;   exit (1);
      }
    } else if (name == "quad") {
      ifstr >> _quadThreshold;
    } else if (name == "maxAngle") {
      ifstr >> _maxAngleThreshold;
      _maxAngleThreshold *= M_PI / 180.0f;
    } else {
      std::cerr << "Unrecognizable option: " << name << std::endl;   exit (1);
    }
  }
  ifstr.close();

  if (_tflag == -10 || _oflag == -10) {
    std::cerr << "_tflag and _oflag not specified: " << _tflag << ' ' << _oflag << std::endl;
    exit (1);
  }

  //----------------------------------------------------------------------
  std::string sbimages = prefix + std::string("bimages.dat");

  for (int i = 0; i < (int)_timages.size(); ++i) {
    _dict[_timages[i]] = i;
  }

  initOimages();
  initVisdata();

  if (_useBound) {
    initBindexes(sbimages);
  }

  std::cerr << "--------------------------------------------------" << std::endl
    << "--- Summary of specified options ---" << std::endl;
  std::cerr << "# of timages: " << (int)_timages.size();

  if (_tflag == -1) {
    std::cerr << " (range specification)" << std::endl;
  } else {
    std::cerr << " (enumeration)" << std::endl;
  }
  std::cerr << "# of oimages: " << (int)_oimages.size();

  if (_oflag == -1) {
    std::cerr << " (range specification)" << std::endl;
  } else if (0 <= _oflag) {
    std::cerr << " (enumeration)" << std::endl;
  } else if (_oflag == -2) {
    std::cerr << " (vis.dat is used)" << std::endl;
  } else if (_oflag == -3) {
    std::cerr << " (not used)" << std::endl;
  }

  std::cerr << "level: " << _level << "  csize: " << _csize << std::endl
    << "threshold: " << _threshold << "  wsize: " << _wsize << std::endl
    << "minImageNum: " << _minImageNum << "  CPU: " << _CPU << std::endl
    << "useVisData: " << _useVisData << "  sequence: " << _sequence << std::endl;
  std::cerr << "--------------------------------------------------" << std::endl;
}

void PMVS3::SOption::initOimages(void) {
  if (_oflag != -2) return;

  std::string svisdata = _prefix + std::string("vis.dat");
  std::ifstream ifstr;
  ifstr.open(svisdata.c_str());
  if (!ifstr.is_open()) {
    std::cerr << "No vis.dat although specified to initOimages: " << std::endl
      << svisdata << std::endl;
    exit (1);
  }

  std::string header;  int num2;
  ifstr >> header >> num2;

  _oimages.clear();
  for (int c = 0; c < num2; ++c) {
    int index0;
    auto ite0 = _dict.find(c);
    if (ite0 == _dict.end()) {
      index0 = -1;
    } else {
      index0 = ite0->second;
    }

    int itmp;
    ifstr >> itmp >> itmp;
    for (int i = 0; i < itmp; ++i) {
      int itmp2;
      ifstr >> itmp2;
      if (index0 != -1 && _dict.find(itmp2) == _dict.end()) {
        _oimages.push_back(itmp2);
      }
    }
  }
  ifstr.close();

  std::sort(_oimages.begin(), _oimages.end());
  _oimages.erase(unique(_oimages.begin(), _oimages.end()), _oimages.end());
}

// When do not use vis.dat
void PMVS3::SOption::initVisdata(void) {
  // Case classifications. Set _visdata by using vis.dat or not.
  if (_useVisData == 0) {
    const int tnum = (int)_timages.size();
    const int onum = (int)_oimages.size();
    const int num = tnum + onum;
    _visdata.resize(num);
    _visdata2.resize(num);
    for (int y = 0; y < num; ++y) {
      _visdata[y].resize(num);
      for (int x = 0; x < num; ++x) {
        if (x == y) {
          _visdata[y][x] = 0;
        } else {
          _visdata[y][x] = 1;
          _visdata2[y].push_back(x);
        }
      }
    }
  } else {
    initVisdata2();
  }
}

// Given _timages and _oimages, set _visdata, _visdata2
void PMVS3::SOption::initVisdata2(void) {
  std::string svisdata = _prefix + std::string("vis.dat");

  std::vector<int> images;
  images.insert(images.end(), _timages.begin(), _timages.end());
  images.insert(images.end(), _oimages.begin(), _oimages.end());
  std::map<int, int> dict2;
  for (int i = 0; i < (int)images.size(); ++i) {
    dict2[images[i]] = i;
  }

  std::ifstream ifstr;
  ifstr.open(svisdata.c_str());
  if (!ifstr.is_open()) {
    std::cerr << "No vis.dat although specified to initVisdata2: " << std::endl
         << svisdata << std::endl;
    exit (1);
  }

  std::string header;
  int num2;
  ifstr >> header >> num2;

  _visdata2.resize((int)images.size());
  for (int c = 0; c < num2; ++c) {
    int index0;
    auto ite0 = dict2.find(c);
    if (ite0 == dict2.end()) {
      index0 = -1;
    } else {
      index0 = ite0->second;
    }

    int itmp;
    ifstr >> itmp >> itmp;
    for (int i = 0; i < itmp; ++i) {
      int itmp2;
      ifstr >> itmp2;
      int index1;
      auto ite1 = dict2.find(itmp2);
      if (ite1 == dict2.end()) {
        index1 = -1;
      } else {
        index1 = ite1->second;
      }

      if (index0 != -1 && index1 != -1) {
        _visdata2[index0].push_back(index1);
      }
    }
  }
  ifstr.close();

  const int num = (int)images.size();
  _visdata.clear();
  _visdata.resize(num);
  for (int y = 0; y < num; ++y) {
    _visdata[y].resize(num);
    fill(_visdata[y].begin(), _visdata[y].end(), 0);
    for (int x = 0; x < (int)_visdata2[y].size(); ++x) {
      _visdata[y][_visdata2[y][x]] = 1;
    }
  }

  // check symmetry
  for (int i = 0; i < (int)_visdata.size(); ++i) {
    for (int j = i+1; j < (int)_visdata.size(); ++j) {
      if (_visdata[i][j] != _visdata[j][i]) {
        _visdata[i][j] = _visdata[j][i] = 1;
      }
    }
  }
}

void PMVS3::SOption::initBindexes(const std::string sbimages) {
  if (sbimages.empty()) return;

  _bindexes.clear();
  std::ifstream ifstr;
  ifstr.open(sbimages.c_str());
  if (!ifstr.is_open()) {
    std::cerr << "File not found: " << sbimages << std::endl;
    exit (1);
  }

  std::cerr << "Reading bimages" << std::endl;
  int itmp;
  ifstr >> itmp;
  for (int i = 0; i < itmp; ++i) {
    int itmp0;
    ifstr >> itmp0;

    if (_dict.find(itmp0) != _dict.end()) {
      _bindexes.push_back(_dict[itmp0]);
    }
  }
  ifstr.close();
}
