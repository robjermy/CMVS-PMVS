#pragma once

#include "camera.hpp"
#include "image.hpp"
#include "numeric/vec4.hpp"

namespace Image {
  // CPhoto is an image with camera parameters
  class CPhoto : public CImage, public CCamera {
  public:
    CPhoto();
    virtual ~CPhoto();

    virtual void init(const std::string name, const std::string mname, const std::string cname, const int maxLevel = 1);
    virtual void init(const std::string name, const std::string mname, const std::string ename, const std::string cname, const int maxLevel = 1);

    // grabTex given 2D sampling information
    void grabTex(const int level, const Vec2f& icoord, const Vec2f& xaxis, const Vec2f& yaxis, const int size, std::vector<Vec3f>& tex, const int normalizef = 1) const;
    // grabTex given 3D sampling information
    void grabTex(const int level, const Vec4f& coord, const Vec4f& pxaxis, const Vec4f& pyaxis, const Vec4f& pzaxis, const int size, std::vector<Vec3f>& tex, float& weight, const int normalizef = 1) const;

    inline Vec3f getColor(const float fx, const float fy, const int level) const;
    inline Vec3f getColor(const Vec4f& coord, const int level) const;
    inline int getMask(const Vec4f& coord, const int level) const;
    inline int getEdge(const Vec4f& coord, const int level) const;

    static float idot(const std::vector<Vec3f>& tex0, const std::vector<Vec3f>& tex1);
    static void idotC(const std::vector<Vec3f>& tex0, const std::vector<Vec3f>& tex1, double* idc);

    static void normalize(std::vector<Vec3f>& tex);

    static float ssd(const std::vector<Vec3f>& tex0, const std::vector<Vec3f>& tex1);
  };

  Vec3f CPhoto::getColor(const float fx, const float fy, const int level) const {
    return CImage::getColor(fx, fy, level);
  };

  Vec3f CPhoto::getColor(const Vec4f& coord, const int level) const {
    const Vec3f icoord = project(coord, level);
    return CImage::getColor(icoord[0], icoord[1], level);
  };

  int CPhoto::getMask(const Vec4f& coord, const int level) const {
    if (_masks[level].empty()) return 1;

    const Vec3f icoord = project(coord, level);
    return CImage::getMask(icoord[0], icoord[1], level);
  };

  int CPhoto::getEdge(const Vec4f& coord, const int level) const {
    if (_edges[level].empty()) return 1;

    const Vec3f icoord = project(coord, level);

    if (icoord[0] < 0 || _widths[level] - 1 <= icoord[0] || icoord[1] < 0 || _heights[level] - 1 <= icoord[1]) return 0;

    return CImage::getEdge(icoord[0], icoord[1], level);
  };
}
