#pragma once

#include <map>

#include "photo.hpp"

namespace Image {
  class CPhotoSetS {
  public:
    CPhotoSetS();
    virtual ~CPhotoSetS();

    void init(const std::vector<int>& images, const std::string prefix, const int maxLevel, const int size, const int alloc);

    // grabTex given 2D sampling information
    void grabTex(const int index, const int level, const Vec2f& icoord, const Vec2f& xaxis, const Vec2f& yaxis, std::vector<Vec3f>& tex, const int normalizef = 1) const;

    // grabTex given 3D sampling information
    void grabTex(const int index, const int level, const Vec4f& coord, const Vec4f& pxaxis, const Vec4f& pyaxis, const Vec4f& pzaxis, std::vector<Vec3f>& tex, float& weight, const int normalizef = 1) const;

    void write(const std::string outdir);
    void free(void);
    void free(const int level);

    void setEdge(const float threshold);

    inline Vec3f project(const int index, const Vec4f& coord, const int level) const;
    inline Vec3f mult(const int index, const Vec4f& coord, const int level) const;

    inline int getWidth(const int index, const int level) const;
    inline int getHeight(const int index, const int level) const;

    inline Vec3f getColor(const Vec4f& coord, const int index, const int level) const;
    inline Vec3f getColor(const int index, const float fx, const float fy, const int level) const;
    inline Vec3f getColor(const int index, const int ix, const int iy, const int level) const;

    inline int getMask(const Vec4f& coord, const int level) const;
    inline int getMask(const Vec4f& coord, const int index, const int level) const;
    inline int getMask(const int index, const float fx, const float fy, const int level) const;
    inline int getMask(const int index, const int ix, const int iy, const int level) const;

    inline int getEdge(const Vec4f& coord, const int index, const int level) const;
    inline int getEdge(const int index, const float fx, const float fy, const int level) const;
    inline int getEdge(const int index, const int ix, const int iy, const int level) const;

    static float incc(const std::vector<std::vector<Vec3f> >& texs, const std::vector<float>& weights);

    int checkAngles(const Vec4f& coord, const std::vector<int>& indexes, const float minAngle, const float maxAngle, const int tau) const;

    void getMinMaxAngles(const Vec4f& coord, const std::vector<int>& indexes, float& minAngle, float& maxAngle) const;

    float computeDepth(const int index, const Vec4f& coord) const;

    // Take care of indexes
    std::vector<int> _images;
    std::vector<CPhoto> _photos;

    const CPhoto& Photo(const int index) const { return _photos[index]; }
    // CPhoto& Photo(const int index) { return _photos[index]; }

    int image2index(const int image) const;
    std::map<int, int> _dict;

    // Number of cameras.
    int _num;
    // Root directory
    std::string _prefix;
    // maximum level
    int _maxLevel;
    // Window size used to refine location
    int _size;

    // getPAxes
    void getPAxes(const int index, const Vec4f& coord, const Vec4f& normal, Vec4f& pxaxis, Vec4f& pyaxis) const;

    // pairwise distance based on optical center and viewing direction
    void setDistances(void);
    std::vector<std::vector<float> > _distances;
  };

  Vec3f CPhotoSetS::project(const int index, const Vec4f& coord, const int level) const{
    return _photos[index].project(coord, level);
  };

  Vec3f CPhotoSetS::mult(const int index, const Vec4f& coord, const int level) const{
    return _photos[index].mult(coord, level);
  };

  int CPhotoSetS::getWidth(const int index, const int level) const {
    return _photos[index].getWidth(level);
  };

  int CPhotoSetS::getHeight(const int index, const int level) const {
    return _photos[index].getHeight(level);
  };

  Vec3f CPhotoSetS::getColor(const Vec4f& coord, const int index, const int level) const {
    return _photos[index].getColor(coord, level);
  };

  Vec3f CPhotoSetS::getColor(const int index, const float fx, const float fy, const int level) const {
    return _photos[index].Image::CImage::getColor(fx, fy, level);
  };

  Vec3f CPhotoSetS::getColor(const int index, const int ix, const int iy, const int level) const {
    return _photos[index].Image::CImage::getColor(ix, iy, level);
  };

  int CPhotoSetS::getMask(const Vec4f& coord, const int level) const {
    for (int index = 0; index < _num; ++index) {
      if (getMask(coord, index, level) == 0) {
        return 0;
      }
    }
    return 1;
  };

  int CPhotoSetS::getMask(const Vec4f& coord, const int index, const int level) const {
    return _photos[index].getMask(coord, level);
  };

  int CPhotoSetS::getMask(const int index, const float fx, const float fy, const int level) const {
    return _photos[index].Image::CImage::getMask(fx, fy, level);
  };

  int CPhotoSetS::getMask(const int index, const int ix, const int iy, const int level) const {
    return _photos[index].Image::CImage::getMask(ix, iy, level);
  };

  int CPhotoSetS::getEdge(const Vec4f& coord, const int index, const int level) const {
    return _photos[index].getEdge(coord, level);
  };

  int CPhotoSetS::getEdge(const int index, const float fx, const float fy, const int level) const {
    return _photos[index].Image::CImage::getEdge(fx, fy, level);
  };

  int CPhotoSetS::getEdge(const int index, const int ix, const int iy, const int level) const {
    return _photos[index].Image::CImage::getEdge(ix, iy, level);
  };
}
