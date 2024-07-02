
#include <cmath>

#include "vec.h"

namespace vortex {

template <typename R, typename S>
vecd<typename result_of<R, S>::type> operator+(const vecd<R>& x,
                                               const vecd<S>& y) {
  typedef typename result_of<R, S>::type T;
  ASSERT(x.m() == x.m());
  vecd<T> z(x.m());
  for (int i = 0; i < z.m(); i++) z(i) = x(i) + y(i);
  return z;
}

template <typename R, typename S>
vecd<typename result_of<R, S>::type> operator-(const vecd<R>& x,
                                               const vecd<S>& y) {
  typedef typename result_of<R, S>::type T;
  ASSERT(x.m() == x.m());
  vecd<T> z(x.m());
  for (int i = 0; i < z.m(); i++) z(i) = x(i) - y(i);
  return z;
}

template <typename R, typename S>
vecd<typename result_of<R, S>::type> operator*(const R& x, const vecd<S>& y) {
  typedef typename result_of<R, S>::type T;
  vecd<T> z(y.m());
  for (int i = 0; i < z.m(); i++) z(i) = x * y(i);
  return z;
}

template <typename R, typename S>
vecd<typename result_of<R, S>::type> operator*(const vecd<R>& x, const S& y) {
  typedef typename result_of<R, S>::type T;
  vecd<T> z(x.m());
  for (int i = 0; i < z.m(); i++) z(i) = x(i) * y;
  return z;
}

/**
 * \brief Computes the (static) vector addition x + y
 */
#define INSTANTIATE_VECADD(R, S, T)                                \
  template <int M>                                                 \
  vecs<M, T> operator+(const vecs<M, R>& u, const vecs<M, S>& v) { \
    vecs<M, T> w;                                                  \
    for (int i = 0; i < M; i++) w(i) = u(i) + v(i);                \
    return w;                                                      \
  }

/**
 * \brief Computes the (static) vector subtraction x - y
 */
#define INSTANTIATE_VECSUB(R, S, T)                                \
  template <int M>                                                 \
  vecs<M, T> operator-(const vecs<M, R>& u, const vecs<M, S>& v) { \
    vecs<M, T> w;                                                  \
    for (int i = 0; i < M; i++) w(i) = u(i) - v(i);                \
    return w;                                                      \
  }

/**
 * \brief Unary - operator, to write vector expression such as -u
 */
#define INSTANTIATE_VECMINUS(R)               \
  template <int M>                            \
  vecs<M, R> operator-(const vecs<M, R>& u) { \
    vecs<M, R> w;                             \
    for (int i = 0; i < M; i++) w(i) = -u(i); \
    return w;                                 \
  }

/**
 * \brief Unary += operator, to write vector expression such as u += v
 */
#define INSTANTIATE_VECINC(S, T)                               \
  template <int M>                                             \
  vecs<M, T>& operator+=(vecs<M, S>& u, const vecs<M, T>& v) { \
    for (int i = 0; i < M; i++) u(i) += v(i);                  \
    return u;                                                  \
  }

/**
 * \brief Unary += operator, to write vector expression such as u -= v
 */
#define INSTANTIATE_VECDEC(S, T)                               \
  template <int M>                                             \
  vecs<M, S>& operator-=(vecs<M, S>& u, const vecs<M, T>& v) { \
    for (int i = 0; i < M; i++) u(i) -= v(i);                  \
    return u;                                                  \
  }

/**
 * \brief Performs component-wise vector multiplication.
 *        This is often needed in graphics when evaluating diffuse or specular
 * light terms but will not be used often for manipulating meshes.
 */
#define INSTANTIATE_VECVECMUL(R, S, T)                             \
  template <int M>                                                 \
  vecs<M, T> operator*(const vecs<M, R>& u, const vecs<M, S>& v) { \
    vecs<M, T> w;                                                  \
    for (int i = 0; i < M; i++) w(i) = u(i) * v(i);                \
    return w;                                                      \
  }

/**
 * \brief Computes the vector-scalar multiplication x * a
 */
#define INSTANTIATE_VECSCAMUL_R(R, S, T)                  \
  template <int M>                                        \
  vecs<M, T> operator*(const vecs<M, R>& u, const S& a) { \
    vecs<M, T> v;                                         \
    for (int i = 0; i < M; i++) v(i) = a * u(i);          \
    return v;                                             \
  }

/**
 * \brief Computes the scalar-vector multiplication a * x
 */
#define INSTANTIATE_VECSCAMUL_L(R, S, T)                  \
  template <int M>                                        \
  vecs<M, T> operator*(const R& a, const vecs<M, S>& u) { \
    vecs<M, T> v;                                         \
    for (int i = 0; i < M; i++) v(i) = a * u(i);          \
    return v;                                             \
  }

/**
 * \brief Computes the vector-scalar multiplication x * a
 */
#define INSTANTIATE_VECSCADIV(R, S, T)                    \
  template <int M>                                        \
  vecs<M, T> operator/(const vecs<M, R>& u, const S& a) { \
    vecs<M, T> v;                                         \
    ASSERT(a != 0) << "divide by zero";                   \
    for (int i = 0; i < M; i++) v(i) = u(i) / a;          \
    return v;                                             \
  }

/**
 * \brief Computes the dot product u.v between two static vectors.
 *
 * \param[in] u - vector (vecs)
 * \param[in] v - vector (vecs)
 *
 * \return dot product u.v
 */
#define INSTANTIATE_DOT(R, S, T)                       \
  template <int M>                                     \
  T dot(const vecs<M, R>& u, const vecs<M, S>& v) {    \
    T result = 0;                                      \
    for (int i = 0; i < M; i++) result += u(i) * v(i); \
    return result;                                     \
  }

/**
 * \brief Computes the outter product uv^T between two static vectors.
 * \param[in] u - vector (vecs)
 * \param[in] v - vector (vecs)
 *
 * \return outter product uv^T (an matrix)
 */
#define INSTANTIATE_OUTER(R, S, T)                                \
  template <int M, int N>                                         \
  mats<M, N, T> outer(const vecs<M, R>& u, const vecs<N, S>& v) { \
    mats<M, N, T> O;                                              \
    for (int i = 0; i < M; i++) {                                 \
      for (int j = 0; j < N; j++) {                               \
        O(i, j) = u(i) * v(j);                                    \
      }                                                           \
    }                                                             \
    return O;                                                     \
  }

/**
 * \brief Computes the length (length) of a static vector
 *
 * \param[in] u - vector to compute the length of
 *
 * \return length of u: || u || = sqrt( u^T u )
 */
template <int M, typename T>
T length(const vecs<M, T>& u) {
  return sqrt(dot(u, u));
}

/**
 * \brief Normalizes a vector so it has a unit length.
 *
 * \param[in,out] u - vector (static) to normalize
 */
#define INSTANTIATE_NORMALIZE(T)                 \
  template <int M>                               \
  vecs<M, T> normalize(const vecs<M, T>& u) {    \
    vecs<M, T> v;                                \
    T n = length(u);                             \
    if (n == 0.0) return u;                      \
    for (int i = 0; i < M; i++) v[i] = u[i] / n; \
    return v;                                    \
  }

#define COMMA ,

INSTANTIATE_VECADD(double, double, double)
INSTANTIATE_VECADD(float, float, float)
INSTANTIATE_VECSUB(double, double, double)
INSTANTIATE_VECSUB(float, float, float)
INSTANTIATE_VECMINUS(double)

INSTANTIATE_VECADD(vecs<2 COMMA double>, vecs<2 COMMA double>,
                   vecs<2 COMMA double>)
INSTANTIATE_VECADD(vecs<3 COMMA double>, vecs<3 COMMA double>,
                   vecs<3 COMMA double>)

INSTANTIATE_VECINC(double, double)
INSTANTIATE_VECDEC(double, double)

INSTANTIATE_VECINC(vecs<2 COMMA double>, vecs<2 COMMA double>)
INSTANTIATE_VECINC(vecs<3 COMMA double>, vecs<3 COMMA double>)

INSTANTIATE_DOT(double, double, double)
INSTANTIATE_DOT(float, float, float)
INSTANTIATE_DOT(int, int, int)
INSTANTIATE_OUTER(double, double, double)

INSTANTIATE_VECSCAMUL_R(double, double, double)
INSTANTIATE_VECSCAMUL_R(float, float, float)
INSTANTIATE_VECSCAMUL_R(vecs<2 COMMA double>, double, vecs<2 COMMA double>)
INSTANTIATE_VECSCAMUL_R(vecs<3 COMMA double>, double, vecs<3 COMMA double>)

INSTANTIATE_VECSCAMUL_L(double, double, double)
INSTANTIATE_VECSCAMUL_L(float, float, float)

INSTANTIATE_VECSCADIV(double, double, double)
INSTANTIATE_VECSCADIV(float, float, float)

INSTANTIATE_NORMALIZE(double)
INSTANTIATE_NORMALIZE(float)

INSTANTIATE_VECVECMUL(double, double, double)

#undef COMMA

}  // namespace vortex
