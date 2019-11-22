#pragma once

#include <set>
#include <vector>

#include "detector.hpp"
#include "numeric/vec3.hpp"
#include "point.hpp"

namespace PMVS3 {
  class CHarris: public CDetector {
  public:
    void run(
      const std::vector<unsigned char>& image,
      const std::vector<unsigned char>& mask,
      const std::vector<unsigned char>& edge,
      const int width, const int height,
      const int gspeedup, const float sigma,
      std::multiset<CPoint> & result
    );

    virtual ~CHarris() {}

  protected:
    float _sigmaD;
    float _sigmaI;

    std::vector<float> _gaussD;
    std::vector<float> _gaussI;

    std::vector<std::vector<Vec3f> > _dIdx;
    std::vector<std::vector<Vec3f> > _dIdy;

    std::vector<std::vector<float> > _dIdxdIdx;
    std::vector<std::vector<float> > _dIdydIdy;
    std::vector<std::vector<float> > _dIdxdIdy;

    std::vector<std::vector<float> > _response;

    void init(
      const std::vector<unsigned char>& image,
      const std::vector<unsigned char>& mask,
      const std::vector<unsigned char>& edge
    );

    void setDerivatives(void);
    void preprocess(void);
    void preprocess2(void);
    void setResponse(void);
  };
}
