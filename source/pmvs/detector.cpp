#include <algorithm>
#include <set>

#include "pmvs/detector.hpp"
#include "pmvs/point.hpp"

void PMVS3::CDetector::setGaussD(const float sigmaD, std::vector<float>& gaussD) {
  const int marginD = (int)ceil(2 * sigmaD);
  const int sizeD = 2 * marginD + 1;

  gaussD.resize(sizeD);

  // set _gaussD
  float denom = 0.0;
  for (int x = 0; x < sizeD; ++x) {
    int xtmp = x - marginD;
    const float dtmp = xtmp * exp(- (xtmp * xtmp) / (2 * sigmaD * sigmaD));
    gaussD[x] = dtmp;
    if (0.0 < dtmp) {
      denom += dtmp;
    }
  }

  for (int x = 0; x < sizeD; ++x) {
    gaussD[x] /= denom;
  }
}

void PMVS3::CDetector::setGaussI(const float sigmaI, std::vector<float>& gaussI) {
  const int marginI = (int)ceil(2 * sigmaI);
  const int sizeI = 2 * marginI + 1;

  gaussI.resize(sizeI);

  // set _gaussI
  float denom = 0.0;
  for (int x = 0; x < sizeI; ++x) {
    int xtmp = x - marginI;
    const float dtmp = exp(- (xtmp * xtmp) / (2 * sigmaI * sigmaI));
    gaussI[x] = dtmp;
    denom += dtmp;
  }

  for (int x = 0; x < sizeI; ++x) {
    gaussI[x] /= denom;
  }
}


float PMVS3::CDetector::setThreshold(std::multiset<CPoint>& grid) {
  auto begin = grid.begin();
  auto end = grid.end();

  float ave = 0.0;
  float ave2 = 0.0;
  int count = 0;

  while (begin != end) {
    count++;
    ave += begin->_response;
    ave2 += begin->_response * begin->_response;
    begin++;
  }

  if (count == 0) {
    count = 1;
  }

  ave /= count;
  ave2 /= count;
  ave2 = sqrt(std::max(0.0f, ave2 - ave * ave));

  const float threshold = ave + ave2;
  return threshold;
}

int PMVS3::CDetector::isCloseBoundary(const int x, const int y, const int margin) const {
  if (_mask.empty()) return 0;

  if (x - margin < 0 || _width <= x + margin || y - margin < 0 || _height <= y + margin) return 1;

  for (int j = -margin; j <= margin; ++j) {
    const int ytmp = y + j;
    for (int i = -margin; i <= margin; ++i) {
      const int xtmp = x + i;
      if (_mask[ytmp][xtmp] == 0) return 1;
    }
  }

  return 0;
}
