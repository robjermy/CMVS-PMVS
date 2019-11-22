#include <fstream>
#include <algorithm>
#include "photoSetS.hpp"
#define _USE_MATH_DEFINES
#include <math.h>

Image::CPhotoSetS::CPhotoSetS() {}
Image::CPhotoSetS::~CPhotoSetS() {}

void Image::CPhotoSetS::init(const std::vector<int>& images, const std::string prefix, const int maxLevel, const int size, const int alloc) {
  _images = images;
  _num = (int)images.size();

  for (int i = 0; i < (int)images.size(); ++i) {
    _dict[images[i]] = i;
  }

  _prefix = prefix;
  _maxLevel = std::max(1, maxLevel);
  _photos.resize(_num);
  std::cerr << "Reading images: " << std::flush;
  for (int index = 0; index < _num; ++index) {
    const int image = _images[index];

    char test0[1024], test1[1024];
    char test2[1024], test3[1024];
    sprintf(test0, "%svisualize/%08d.ppm", prefix.c_str(), image);
    sprintf(test1, "%svisualize/%08d.jpg", prefix.c_str(), image);
    sprintf(test2, "%svisualize/%08d.png", prefix.c_str(), image);
    sprintf(test3, "%svisualize/%08d.tiff", prefix.c_str(), image);
    if (std::ifstream(test0) || std::ifstream(test1)
#if defined(PMVS_HAVE_PNG)
      || std::ifstream(test2)
#endif
#if defined(PMVS_HAVE_TIFF)
      || std::ifstream(test3)
#endif
    ) {
      char name[1024], mname[1024], ename[1024], cname[1024];

      // Set name
      sprintf(name, "%svisualize/%08d", prefix.c_str(), image);
      sprintf(mname, "%smasks/%08d", prefix.c_str(), image);
      sprintf(ename, "%sedges/%08d", prefix.c_str(), image);
      sprintf(cname, "%stxt/%08d.txt", prefix.c_str(), image);

      _photos[index].init(name, mname, ename, cname, _maxLevel);
      if (alloc) {
        _photos[index].alloc();
      } else {
        _photos[index].alloc(1);
      }
      std::cerr << '*' << std::flush;
    } else { // try 4 digits
      char name[1024], mname[1024], ename[1024], cname[1024];

      // Set name
      sprintf(name, "%svisualize/%04d", prefix.c_str(), image);
      sprintf(mname, "%smasks/%04d", prefix.c_str(), image);
      sprintf(ename, "%sedges/%04d", prefix.c_str(), image);
      sprintf(cname, "%stxt/%04d.txt", prefix.c_str(), image);

      _photos[index].init(name, mname, ename, cname, _maxLevel);
      if (alloc) {
        _photos[index].alloc();
      } else {
        _photos[index].alloc(1);
      }
      std::cerr << '*' << std::flush;
    }
  }
  std::cerr << std::endl;
  const int margin = size / 2;
  _size = 2 * margin + 1;
}

void Image::CPhotoSetS::free(void) {
  for (int index = 0; index < (int)_photos.size(); ++index) {
    _photos[index].free();
  }
}

void Image::CPhotoSetS::free(const int level) {
  for (int index = 0; index < (int)_photos.size(); ++index) {
    _photos[index].free(level);
  }
}

void Image::CPhotoSetS::setEdge(const float threshold) {
  for (int index = 0; index < _num; ++index) {
    _photos[index].setEdge(threshold);
  }
}

void Image::CPhotoSetS::write(const std::string outdir) {
  for (int index = 0; index < _num; ++index) {
    const int image = _images[index];
    char buffer[1024];
    sprintf(buffer, "%s%08d.txt", outdir.c_str(), image);

    _photos[index].write(buffer);
  }
}

// get x and y axis to collect textures given reference index and normal
void Image::CPhotoSetS::getPAxes(const int index, const Vec4f& coord, const Vec4f& normal, Vec4f& pxaxis, Vec4f& pyaxis) const{
  _photos[index].getPAxes(coord, normal, pxaxis, pyaxis);
}

void Image::CPhotoSetS::grabTex(const int index, const int level, const Vec2f& icoord, const Vec2f& xaxis, const Vec2f& yaxis, std::vector<Vec3f>& tex, const int normalizef) const{
  _photos[index].grabTex(level, icoord, xaxis, yaxis, _size, tex, normalizef);
}

// grabTex given 3D sampling information
void Image::CPhotoSetS::grabTex(const int index, const int level, const Vec4f& coord, const Vec4f& pxaxis, const Vec4f& pyaxis, const Vec4f& pzaxis, std::vector<Vec3f>& tex, float& weight, const int normalizef) const {
  _photos[index].grabTex(level, coord, pxaxis, pyaxis, pzaxis, _size, tex, weight, normalizef);
}

float Image::CPhotoSetS::incc(const std::vector<std::vector<Vec3f> >& texs, const std::vector<float>& weights) {
  float incctmp = 0.0;
  float denom = 0.0;
  for (int i = 0; i < (int)weights.size(); ++i) {
    if (texs[i].empty()) continue;

    for (int j = i+1; j < (int)weights.size(); ++j) {
      if (texs[j].empty()) continue;

      const float weight = weights[i] * weights[j];
      const float ftmp = CPhoto::idot(texs[i], texs[j]);
      incctmp += ftmp * weight;
      denom += weight;
    }
  }

  if (denom == 0.0) {
    return 2.0f;
  }

  return incctmp / denom;
}

void Image::CPhotoSetS::getMinMaxAngles(const Vec4f& coord, const std::vector<int>& indexes, float& minAngle, float& maxAngle) const {
  minAngle = M_PI;
  maxAngle = 0.0f;
  std::vector<Vec4f> rays;  rays.resize((int)indexes.size());
  for (int i = 0; i < (int)indexes.size(); ++i) {
    const int index = indexes[i];
    rays[i] = _photos[index].OpticalCenter() - coord;
    unitize(rays[i]);
  }

  for (int i = 0; i < (int)indexes.size(); ++i) {
    for (int j = i+1; j < (int)indexes.size(); ++j) {
      const float dot = std::max(-1.0f, std::min(1.0f, rays[i] * rays[j]));
      const float angle = acos(dot);
      minAngle = std::min(angle, minAngle);
      maxAngle = std::max(angle, maxAngle);
    }
  }
}

int Image::CPhotoSetS::checkAngles(const Vec4f& coord, const std::vector<int>& indexes, const float minAngle, const float maxAngle, const int num) const {
  int count = 0;

  std::vector<Vec4f> rays;  rays.resize((int)indexes.size());
  for (int i = 0; i < (int)indexes.size(); ++i) {
    const int index = indexes[i];
    rays[i] = _photos[index].OpticalCenter() - coord;
    unitize(rays[i]);
  }

  for (int i = 0; i < (int)indexes.size(); ++i) {
    for (int j = i+1; j < (int)indexes.size(); ++j) {
      const float dot = std::max(-1.0f, std::min(1.0f, rays[i] * rays[j]));
      const float angle = acos(dot);
      if (minAngle < angle && angle < maxAngle) {
        ++count;
      }
    }
  }

  if (count < 1) {
    return 1;
  }

  return 0;
}

float Image::CPhotoSetS::computeDepth(const int index, const Vec4f& coord) const {
  return _photos[index].computeDepth(coord);
}

void Image::CPhotoSetS::setDistances(void) {
  _distances.resize(_num);
  float avedis = 0.0f;
  int denom = 0;
  for (int i = 0; i < _num; ++i) {
    _distances[i].resize(_num);
    for (int j = 0; j < _num; ++j) {
      if (i == j) {
        _distances[i][j] = 0.0f;
      } else {
        const float ftmp = norm(_photos[i].OpticalCenter() - _photos[j].OpticalCenter());
        _distances[i][j] = ftmp;
        avedis += ftmp;
        denom++;
      }
    }
  }

  if (denom == 0) return;

  avedis /= denom;
  if (avedis == 0.0f) {
    std::cerr << "All the optical centers are identical..?" << std::endl;
    exit (1);
  }

  // plus angle difference
  for (int i = 0; i < _num; ++i) {
    Vec4f ray0 = _photos[i].OpticalAxis();
    ray0[3] = 0.0f;
    for (int j = 0; j < _num; ++j) {
      Vec4f ray1 = _photos[j].OpticalAxis();
      ray1[3] = 0.0f;

      _distances[i][j] /= avedis;
      const float margin = cos(10.0f * M_PI / 180.0f);
      const float dis = std::max(0.0f, 1.0f - ray0 * ray1 - margin);
      _distances[i][j] += dis;
    }
  }
}

int Image::CPhotoSetS::image2index(const int image) const {
  std::map<int, int>::const_iterator pos = _dict.find(image);
  if (pos == _dict.end()) {
    return -1;
  }

  return pos->second;
}
