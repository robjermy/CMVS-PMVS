#pragma once

#include "vec2.hpp"

template<typename T>
class TMat2
{
 private:
  TVec2<T> _row[2];

 public:
  // Standard constructors
  //
  TMat2() { *this = 0.0; }
  TMat2(T a, T b, T c, T d) {
    _row[0][0]=a;
    _row[0][1]=b;
    _row[1][0]=c;
    _row[1][1]=d;
  }
  TMat2(const TVec2<T> &r0,const TVec2<T> &r1) {
    _row[0] = r0;
    _row[1] = r1;
  }
  TMat2(const TMat2 &m) { *this = m; }

  // Descriptive interface
  typedef T value_type;
  typedef TVec2<T> vector_type;
  typedef TMat2 inverse_type;
  static int dim() { return 2; }

  // Access methods - note: A(i, j) == row i, col j
  T& operator()(int i, int j) {
    return _row[i][j];
  }
  T operator()(int i, int j) const {
    return _row[i][j];
  }
  TVec2<T>& operator[](int i) {
    return _row[i];
  }
  const TVec2<T>& operator[](int i) const {
    return _row[i];
  }
  inline TVec2<T> col(int i) const {
    return TVec2<T>(
      _row[0][i],
      _row[1][i]
    );
  }

  operator T*() {
    return _row[0];
  }
  operator const T*() {
    return _row[0];
  }
  operator const T*() const {
    return _row[0];
  }

  // Assignment methods
  //
  inline TMat2& operator=(const TMat2& m);
  inline TMat2& operator=(T s);

  inline TMat2& operator+=(const TMat2& m);
  inline TMat2& operator-=(const TMat2& m);
  inline TMat2& operator*=(T s);
  inline TMat2& operator/=(T s);

  inline bool operator==(const TMat2& m) const;
  inline bool operator!=(const TMat2& m) const;

  // Construction of standard matrices
  static TMat2 I();
  static TMat2 outer_product(const TVec2<T> &u, const TVec2<T> &v) {
    return TMat2(
      u[0]*v[0],
      u[0]*v[1],
      u[1]*v[0],
      u[1]*v[1]
    );
  }
  static TMat2 outer_product(const TVec2<T> &u) {
    return outer_product(u,u);
  }

  TMat2 &diag(T d);
  TMat2 &ident() { return diag(1.0); }
};

////////////////////////////////////////////////////////////////////////
//
// Method definitions
//

template<typename T>
inline TMat2<T>& TMat2<T>::operator=(const TMat2<T>& m) {
  _row[0] = m[0];
  _row[1] = m[1];
  return *this;
}

template<typename T>
inline TMat2<T>& TMat2<T>::operator=(T s) {
  _row[0] = s;
  _row[1] = s;
  return *this;
}

template<typename T>
inline TMat2<T>& TMat2<T>::operator+=(const TMat2<T>& m) {
  _row[0] += m._row[0];
  _row[1] += m._row[1];
  return *this;
}

template<typename T>
inline TMat2<T>& TMat2<T>::operator-=(const TMat2<T>& m) {
  _row[0] -= m._row[0];
  _row[1] -= m._row[1];
  return *this;
}

template<typename T>
inline TMat2<T>& TMat2<T>::operator*=(T s) {
  _row[0] *= s;
  _row[1] *= s;
  return *this;
}

template<typename T>
inline TMat2<T>& TMat2<T>::operator/=(T s) {
  _row[0] /= s;
  _row[1] /= s;
  return *this;
}

template<typename T>
inline bool TMat2<T>::operator==(const TMat2<T>& m) const {
  return _row[0] == m._row[0] && _row[1] == m._row[1];
}

template<typename T>
inline bool TMat2<T>::operator!=(const TMat2<T>& m) const {
  return !(*this == m);
}

////////////////////////////////////////////////////////////////////////
//
// Operator definitions
//

template<typename T>
inline TMat2<T> operator+(const TMat2<T> &n, const TMat2<T> &m) {
  return TMat2<T>(
    n[0] + m[0],
    n[1] + m[1]
  );
}

template<typename T>
inline TMat2<T> operator-(const TMat2<T> &n, const TMat2<T> &m) {
  return TMat2<T>(
    n[0] - m[0],
    n[1] - m[1]
  );
}

template<typename T>
inline TMat2<T> operator-(const TMat2<T> &m) {
  return TMat2<T>(
    -m[0],
    -m[1]
  );
}

template<typename T>
inline TMat2<T> operator*(T s, const TMat2<T> &m) {
  return TMat2<T>(
    m[0]*s,
    m[1]*s
  );
}

template<typename T>
inline TMat2<T> operator*(const TMat2<T> &m, T s) {
  return s*m;
}

template<typename T>
inline TMat2<T> operator/(const TMat2<T> &m, T s) {
  return TMat2<T>(
    m[0]/s,
    m[1]/s
  );
}

template<typename T>
inline TVec2<T> operator*(const TMat2<T> &m, const TVec2<T> &v) {
  return TVec2<T>(
    m[0]*v,
    m[1]*v
  );
}

template<typename T>
extern TMat2<T> operator*(const TMat2<T> &n, const TMat2<T> &m);

template<typename T>
inline std::ostream &operator<<(std::ostream &out, const TMat2<T>& M) {
  return out << M[0] << std::endl  << M[1];
}

template<typename T>
inline std::istream &operator>>(std::istream &in, TMat2<T>& M) {
  return in >> M[0] >> M[1];
}

////////////////////////////////////////////////////////////////////////
// Misc. function definitions

template<typename T>
inline T det(const TMat2<T> &m) {
  return m(0,0)*m(1,1) - m(0,1)*m(1,0);
}

template<typename T>
inline T trace(const TMat2<T> &m) {
  return m(0,0) + m(1,1);
}

template<typename T>
inline TMat2<T> transpose(const TMat2<T> &m) {
  return TMat2<T>(m.col(0), m.col(1));
}

template<typename T>
inline TMat2<T> adjoint(const TMat2<T> &m) {
  return TMat2<T>(perp(m[1]), -perp(m[0]));
}

template<typename T>
TMat2<T> TMat2<T>::I() {
  return TMat2<T>(
    1,0,
    0,1
  );
}

template<typename T>
TMat2<T> &TMat2<T>::diag(T d) {
    _row[0][0] = d; _row[0][1] = 0;
    _row[1][0] = 0; _row[1][1] = d;

    return *this;
}

template<typename T>
TMat2<T> operator*(const TMat2<T> &n, const TMat2<T>& m) {
    TMat2<T> A;

    for(int i = 0; i < 2; ++i) {
      for(int j = 0; j < 2; ++j) {
        A(i,j) = n[i]*m.col(j);
      }
    }
    return A;
}

template<typename T>
T invert(TMat2<T> &inv, const TMat2<T> &m) {
  T d = det(m);

  if( d==0.0 )
    return 0.0;

  inv(0, 0) =  m(1,1)/d;
  inv(0, 1) = -m(0,1)/d;
  inv(1, 0) = -m(1,0)/d;
  inv(1, 1) =  m(0,0)/d;

  return d;
}

template<typename T>
bool eigenvalues(const TMat2<T>& M, TVec2<T>& evals) {
  T B = -M(0,0)-M(1,1);
  T C = det(M);

  T dis = B*B - 4.0*C;
  if( dis< 1e-6 ) {
	  return false;
  }

  T s = sqrt(dis);

  evals[0] = 0.5*(-B + s);
  evals[1] = 0.5*(-B - s);
  return true;
}

template<typename T>
bool eigenvectors(const TMat2<T>& M, const TVec2<T>& evals, TVec2<T> evecs[2]) {
  evecs[0] = Vec2(-M(0,1), M(0,0)-evals[0]);
  evecs[1] = Vec2(-M(0,1), M(0,0)-evals[1]);

  unitize(evecs[0]);
  unitize(evecs[1]);

  return true;
}

template<typename T>
bool eigen(const TMat2<T>& M, TVec2<T>& evals, TVec2<T> evecs[2]) {
  bool result = eigenvalues(M, evals);
  if (result) {
	  eigenvectors(M, evals, evecs);
  }
  return result;
}

typedef TMat2<double> Mat2;
typedef TMat2<float> Mat2f;
