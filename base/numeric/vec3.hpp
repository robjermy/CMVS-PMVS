#pragma once

#include <cmath>
#include "vec2.hpp"

template<typename T>
class TVec3 {
private:
  T _elt[3];

 public:
  // Standard constructors
  //
  TVec3(T s = T()) { *this = s; }
  TVec3(T x, T y, T z) {
    _elt[0] = x;
    _elt[1] = y;
    _elt[2] = z;
  }

  // Copy constructors & assignment operators
  template<typename U> TVec3(const TVec3<U>& v) { *this = v; }
  template<typename U> TVec3(const TVec2<U>& v,T w) {
    _elt[0] = v[0];
    _elt[1] = v[1];
    _elt[2] = w;
  }

  template<typename U> TVec3(const U v[3]) {
    _elt[0]=v[0];
    _elt[1]=v[1];
    _elt[2]=v[2];
  }

  template<typename U> TVec3& operator=(const TVec3<U>& v) {
    _elt[0]=v[0];
    _elt[1]=v[1];
    _elt[2]=v[2];
    return *this;
  }

  TVec3& operator=(T s) {
    _elt[0] = _elt[1] = _elt[2] = s;
    return *this;
  }

  // Descriptive interface
  typedef T value_type;
  static int size() { return 3; }

  // Access methods
  operator T*() {
    return _elt;
  }
  operator const T*() const {
    return _elt;
  }

  T& operator[](int i) {
    return _elt[i];
  }
  T operator[](int i) const {
    return _elt[i];
  }
  operator const T*() {
    return _elt;
  }

  // Assignment and in-place arithmetic methods
  inline TVec3& operator+=(const TVec3& v);
  inline TVec3& operator-=(const TVec3& v);
  inline TVec3& operator*=(T s);
  inline TVec3& operator/=(T s);

  inline bool operator==(const TVec3& v) const;
  inline bool operator!=(const TVec3& v) const;

  inline T norm2(void) const {
    return (*this) * (*this);
  }

  inline T norm(void) const {
    return sqrt(norm2());
  }

  inline void unitize(void) {
    const T denom2 = norm2();

    if(denom2 != 1.0 && denom2 != 0.0) {
      const T denom = sqrt(denom2);
      _elt[0] /= denom;
      _elt[1] /= denom;
      _elt[2] /= denom;
    }
  }

  inline T sum(void) const {
    return _elt[0] + _elt[1] + _elt[2];
  }
};

////////////////////////////////////////////////////////////////////////
//
// Method definitions
//

template<typename T> inline TVec3<T>& TVec3<T>::operator+=(const TVec3<T>& v) {
  _elt[0] += v[0];
  _elt[1] += v[1];
  _elt[2] += v[2];
  return *this;
}

template<typename T> inline TVec3<T>& TVec3<T>::operator-=(const TVec3<T>& v) {
  _elt[0] -= v[0];
  _elt[1] -= v[1];
  _elt[2] -= v[2];
  return *this;
}

template<typename T> inline TVec3<T>& TVec3<T>::operator*=(T s) {
  _elt[0] *= s;
  _elt[1] *= s;
  _elt[2] *= s;
  return *this;
}

template<typename T> inline TVec3<T>& TVec3<T>::operator/=(T s) {
  _elt[0] /= s;
  _elt[1] /= s;
  _elt[2] /= s;
  return *this;
}

template<typename T> inline bool TVec3<T>::operator==(const TVec3<T>& v) const{
  return _elt[0] == v._elt[0] && _elt[1] == v._elt[1] && _elt[2] == v._elt[2];
}

template<typename T> inline bool TVec3<T>::operator!=(const TVec3<T>& v) const{
  return !(*this == v);
}

////////////////////////////////////////////////////////////////////////
// Operator definitions
template<typename T>
inline TVec3<T> operator+(const TVec3<T> &u, const TVec3<T>& v) {
  return TVec3<T>(
    u[0]+v[0],
    u[1]+v[1],
    u[2]+v[2]
  );
}

template<typename T>
inline TVec3<T> operator-(const TVec3<T> &u, const TVec3<T>& v) {
  return TVec3<T>(
    u[0]-v[0],
    u[1]-v[1],
    u[2]-v[2]
  );
}

template<typename T> inline TVec3<T> operator-(const TVec3<T> &v) {
  return TVec3<T>(
    -v[0],
    -v[1],
    -v[2]
  );
}

template<typename T, typename N> inline TVec3<T> operator*(N s, const TVec3<T> &v) {
  return TVec3<T>(
    v[0]*s,
    v[1]*s,
    v[2]*s
  );
}
template<typename T, typename N> inline TVec3<T> operator*(const TVec3<T> &v, N s) {
  return s*v;
}

template<typename T, typename N> inline TVec3<T> operator/(const TVec3<T> &v, N s) {
  return TVec3<T>(
    v[0]/s,
    v[1]/s,
    v[2]/s
  );
}

template<typename T> inline T operator*(const TVec3<T> &u, const TVec3<T>& v) {
  return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

template<typename T> inline TVec3<T> cross(const TVec3<T>& u, const TVec3<T>& v)
{
    return TVec3<T>(
      u[1]*v[2] - v[1]*u[2],
      -u[0]*v[2] + v[0]*u[2],
      u[0]*v[1] - v[0]*u[1]
    );
}

template<typename T>
inline TVec3<T> operator^(const TVec3<T>& u, const TVec3<T>& v) {
  return cross(u, v);
}


template<typename T>
inline std::ostream &operator<<(std::ostream &out, const TVec3<T>& v) {
  return out << v[0] << " " << v[1] << " " << v[2];
}

template<typename T>
inline std::istream &operator>>(std::istream &in, TVec3<T>& v) {
  return in >> v[0] >> v[1] >> v[2];
}

////////////////////////////////////////////////////////////////////////
//
// Misc. function definitions
//

template<typename T> inline T norm2(const TVec3<T>& v)  {
  return v*v;
}

template<typename T> inline T norm(const TVec3<T>& v)   {
  return sqrt(norm2(v));
}

template<typename T> inline void unitize(TVec3<T>& v) {
  T l = norm2(v);
  if(l != 1.0 && l != 0.0) {
    v /= sqrt(l);
  }
}

template<typename T> inline TVec2<T> proj(const TVec3<T>& v) {
  TVec2<T> u(v[0], v[1]);
    if(v[2] != 1.0 && v[2] != 0.0) {
	    u /= v[2];
    }
    return u;
}

template<typename T> inline void ortho(const TVec3<T>& z, TVec3<T>& x, TVec3<T>& y) {
  if (fabs(z[0]) > 0.5) {
    x[0] = z[1];
    x[1] = -z[0];
    x[2] = 0;
  } else if (fabs(z[1]) > 0.5) {
    x[1] = z[2];
    x[2] = -z[1];
    x[0] = 0;
  } else {
    x[2] = z[0];
    x[0] = -z[2];
    x[1] = 0;
  }
  unitize(x);
  y = cross(z, x);
}

template<typename T>
bool predVec30(const TVec3<T>& lhs, const TVec3<T>& rhs) {
  return lhs[0] < rhs[0];
}

template<typename T>
bool predVec31(const TVec3<T>& lhs, const TVec3<T>& rhs) {
  return lhs[1] < rhs[1];
}

template<typename T>
bool predVec32(const TVec3<T>& lhs, const TVec3<T>& rhs) {
  return lhs[2] < rhs[2];
}

typedef TVec3<double> Vec3;
typedef TVec3<float>  Vec3f;
typedef TVec3<int>    Vec3i;

template<typename T>
struct Svec3cmp {
  bool operator()(const TVec3<T>& lhs, const TVec3<T>& rhs) const {
    return lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]) || (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] < rhs[2]);
  }
};
