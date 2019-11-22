#pragma once

#include <algorithm>
#include <climits>
#include <string>
#include <vector>

#include "numeric/mat3.hpp"
#include "numeric/mat4.hpp"
#include "numeric/vec4.hpp"

namespace Image {
  class CCamera {
  public:
    CCamera();
    virtual ~CCamera();

    // Update projection matrices from intrinsics and extrinsics
    void updateProjection(void);

    // Update all the camera related parameters
    void updateCamera(void);

    virtual void init(const std::string& cname, const int maxLevel);
    void write(const std::string file);

    inline Vec3f project(const Vec4f& coord, const int level) const;
    inline Vec3f mult(const Vec4f& coord, const int level) const;

    static void setProjection(const std::vector<float>& intrinsics, const std::vector<float>& extrinsics, std::vector<Vec4f>& projection, const int txtType);

    float getScale(const Vec4f& coord, const int level) const;
    void getPAxes(const Vec4f& coord, const Vec4f& normal, Vec4f& pxaxis, Vec4f& pyaxis, const int level = 0) const;

    void setAxesScale(const float axesScale);

    static void proj2q(Mat4& mat, double q[6]);
    static void q2proj(const double q[6], Mat4& mat);
    static void setProjectionSub(double params[], std::vector<Vec4f>& projection, const int level);

    float computeDistance(const Vec4f& point) const;
    float computeDepth(const Vec4f& point) const;
    float computeDepthDif(const Vec4f& rhs, const Vec4f& lhs) const;

    // Compute where the viewing ray passing through coord intersects
    // with the plane abcd.
    Vec4f intersect(const Vec4f& coord, const Vec4f& abcd) const;
    void intersect(const Vec4f& coord, const Vec4f& abcd, Vec4f& cross, float& distance) const;
    // Computer a 3D coordinate that projects to a given image
    // coordinate. You can specify a different depth by the third
    // component of icoord.
    Vec4f unproject(const Vec3f& icoord, const int _level) const;

    void setK(Mat3f& K) const;
    void setRT(Mat4f& RT) const;

    void getR(Mat3f& R) const;

    // Getters
    const Vec4f OpticalCenter() const { return _center; }
    const Vec4f OpticalAxis() const { return _oaxis; }
    const std::vector<std::vector<Vec4f>> ProjectionMatrix() const { return _projection; }

  protected: // variables
    std::string _cname; // txt file name

    Vec4f _center; // Optical center
    Vec4f _oaxis; // Optical axis
    Vec3f _xaxis;
    Vec3f _yaxis;
    Vec3f _zaxis;

    std::vector<std::vector<Vec4f>> _projection; // 3x4 projection matrix
    float _ipscale;

    // intrinsic and extrinsic camera parameters. Compact form.
    std::vector<float> _intrinsics;
    std::vector<float> _extrinsics;

    int _txtType; // camera parameter type
    int _maxLevel;

    float _axesScale;

  protected: // methods
    Vec4f getOpticalCenter(void) const;
  };

  inline Vec3f CCamera::project(const Vec4f& coord, const int level) const {
    Vec3f vtmp;
    for (int i = 0; i < 3; ++i) {
      vtmp[i] = _projection[level][i] * coord;
    }

    if (vtmp[2] <= 0.0) {
      vtmp[0] = -0xffff;
      vtmp[1] = -0xffff;
      vtmp[2] = -1.0f;
      return vtmp;
    } else {
      vtmp /= vtmp[2];
    }

    vtmp[0] = std::max((float)(INT_MIN + 3.0f), std::min((float)(INT_MAX - 3.0f), vtmp[0]));
    vtmp[1] = std::max((float)(INT_MIN + 3.0f), std::min((float)(INT_MAX - 3.0f), vtmp[1]));

    return vtmp;
  };

  inline Vec3f CCamera::mult(const Vec4f& coord, const int level) const {
    Vec3f vtmp;
    for (int i = 0; i < 3; ++i) {
      vtmp[i] = _projection[level][i] * coord;
    }

    return vtmp;
  };

  template<class T>
  float computeEPD(const TMat3<T>& F, const TVec3<T>& p0, const TVec3<T>& p1) {
    TVec3<T> line = F * p1;
    const T ftmp = sqrt(line[0] * line[0] + line[1] * line[1]);
    if (ftmp == 0.0) return 0.0;

    line /= ftmp;
    return fabs(line * p0);
  };

  template<class T>
  void setF(const Image::CCamera& lhs, const Image::CCamera& rhs,
      TMat3<T>& F, const int level = 0) {
    const TVec4<T>& p00 = lhs.ProjectionMatrix()[level][0];
    const TVec4<T>& p01 = lhs.ProjectionMatrix()[level][1];
    const TVec4<T>& p02 = lhs.ProjectionMatrix()[level][2];

    const TVec4<T>& p10 = rhs.ProjectionMatrix()[level][0];
    const TVec4<T>& p11 = rhs.ProjectionMatrix()[level][1];
    const TVec4<T>& p12 = rhs.ProjectionMatrix()[level][2];

    F[0][0] = det(TMat4<T>(p01, p02, p11, p12));
    F[0][1] = det(TMat4<T>(p01, p02, p12, p10));
    F[0][2] = det(TMat4<T>(p01, p02, p10, p11));

    F[1][0] = det(TMat4<T>(p02, p00, p11, p12));
    F[1][1] = det(TMat4<T>(p02, p00, p12, p10));
    F[1][2] = det(TMat4<T>(p02, p00, p10, p11));

    F[2][0] = det(TMat4<T>(p00, p01, p11, p12));
    F[2][1] = det(TMat4<T>(p00, p01, p12, p10));
    F[2][2] = det(TMat4<T>(p00, p01, p10, p11));
  };
};
