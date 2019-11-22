#pragma once

#include <set>
#include <vector>

#include "detector.hpp"
#include "numeric/vec3.hpp"
#include "point.hpp"

namespace PMVS3 {
  class CDifferenceOfGaussians: public CDetector {
  public:
    void run (
      const std::vector<unsigned char>& image,
      const std::vector<unsigned char>& mask,
      const std::vector<unsigned char>& edge,
      const int width,	    const int height,
      const int gspeedup,
      const float firstScale,   // 1.4f
      const float lastScale,    // 4.0f
      std::multiset<CPoint> & result
    );

    virtual ~CDifferenceOfGaussians() {}

  protected:
    float _firstScale;
    float _lastScale;

    void init(
      const std::vector<unsigned char>& image,
      const std::vector<unsigned char>& mask,
      const std::vector<unsigned char>& edge
    );

    void setRes(const float sigma, std::vector<std::vector<float> >& res);

    static int isLocalMax(
      const std::vector<std::vector<float> >& pdog,
      const std::vector<std::vector<float> >& cdog,
      const std::vector<std::vector<float> >& ndog,
      const int x, const int y
    );

    static int isLocalMax(const std::vector<std::vector<float> >& dog, const int x, const int y);

    static int notOnEdge(const std::vector<std::vector<float> >& dog, int x, int y);

    static float getResponse(
      const std::vector<std::vector<float>>& pdog,
      const std::vector<std::vector<float>>& cdog,
      const std::vector<std::vector<float>>& ndog,
      const int x, const int y
    );

    static float getResponse(const std::vector<std::vector<float> >& dog, const int x, const int y);

    static void setDOG(
      const std::vector<std::vector<float> >& cres,
      const std::vector<std::vector<float> >& nres,
      std::vector<std::vector<float> >& dog
    );
  };
}
