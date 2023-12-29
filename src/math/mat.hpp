#include "linalg.h"
#include "log.h"
#include "mat.h"
#include "sym.h"
#include "vec.h"

namespace vortex {

template <typename T>
template <typename S>
matd<T>::matd(const symd<S>& M) : matd<T>(M.n(), M.n()) {
  for (int i = 0; i < m_; i++)
    for (int j = 0; j < n_; j++) (*this)(i, j) = M(i, j);
}

template <int M, int N, typename T>
template <typename S>
mats<M, N, T>::mats(const syms<M, S>& A) : mats<M, N, T>() {
  static_assert(M == N, "matrix shoould be square");
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) (*this)(i, j) = A(i, j);
}

template <typename T>
void matd<T>::set_row(int i, const vecd<T>& row) {
  ASSERT(row.m() == n_);
  for (int j = 0; j < row.m(); j++) (*this)(i, j) = row(j);
}

template <typename T>
void matd<T>::get_row(int i, vecd<T>& row) const {
  ASSERT(row.m() == n_);
  for (int j = 0; j < row.m(); j++) row(j) = (*this)(i, j);
}

template <typename R, typename S>
matd<typename result_of<R, S>::type> operator*(const matd<R>& A,
                                               const matd<S>& B) {
  typedef typename result_of<R, S>::type T;
  ASSERT(A.n() == B.m()) << "bad matrix sizes";
  matd<T> C(A.m(), B.n());
  for (int i = 0; i < A.m(); i++) {
    for (int j = 0; j < B.n(); j++) {
      T sum = 0;
      for (int k = 0; k < A.n(); k++) sum += A(i, k) * B(k, j);
      C(i, j) = sum;
    }
  }
  return C;
}

template <typename R, typename S>
matd<typename result_of<R, S>::type> operator+(const matd<R>& A,
                                               const matd<S>& B) {
  typedef typename result_of<R, S>::type T;
  ASSERT(A.m() == B.m()) << "bad matrix sizes";
  ASSERT(A.n() == B.n()) << "bad matrix sizes";
  matd<T> C(A.m(), A.n());
  for (int i = 0; i < C.m(); i++) {
    for (int j = 0; j < C.n(); j++) {
      C(i, j) = A(i, j) + B(i, j);
    }
  }
  return C;
}

template <typename R, typename S>
matd<typename result_of<R, S>::type> operator-(const matd<R>& A,
                                               const matd<S>& B) {
  typedef typename result_of<R, S>::type T;
  ASSERT(A.m() == B.m()) << "bad matrix sizes";
  ASSERT(A.n() == B.n()) << "bad matrix sizes";
  matd<T> C(A.m(), A.n());
  for (int i = 0; i < C.m(); i++) {
    for (int j = 0; j < C.n(); j++) {
      C(i, j) = A(i, j) - B(i, j);
    }
  }
  return C;
}

/**
 * \brief Computes the matrix addition A + B.
 */
#define INSTANTIATE_MATADD(R, S, T)                                         \
  template <int M, int N>                                                   \
  mats<M, N, T> operator+(const mats<M, N, R>& A, const mats<M, N, S>& B) { \
    mats<M, N, T> C;                                                        \
    for (int i = 0; i < M; i++)                                             \
      for (int j = 0; j < N; j++) C(i, j) = A(i, j) + B(i, j);              \
    return C;                                                               \
  }

/**
 * \brief Computes the matrix subtraction A - B.
 */
#define INSTANTIATE_MATSUB(R, S, T)                                         \
  template <int M, int N>                                                   \
  mats<M, N, T> operator-(const mats<M, N, R>& A, const mats<M, N, S>& B) { \
    mats<M, N, T> C;                                                        \
    for (int i = 0; i < M; i++)                                             \
      for (int j = 0; j < N; j++) C(i, j) = A(i, j) - B(i, j);              \
    return C;                                                               \
  }

/**
 * \brief Computes the matrix-matrix multiplication A * B.
 */
#define INSTANTIATE_MATMUL(R, S, T)                            \
  template <int MA, int NA, int MB, int NB>                    \
  mats<MA, NB, T> operator*(const mats<MA, NA, R>& A,          \
                            const mats<MB, NB, S>& B) {        \
    static_assert(NA == MB, "bad matrix sizes");               \
    mats<MA, NB, T> C;                                         \
    for (int i = 0; i < MA; i++) {                             \
      for (int j = 0; j < NB; j++) {                           \
        T sum = 0;                                             \
        for (int k = 0; k < NA; k++) sum += A(i, k) * B(k, j); \
        C(i, j) = sum;                                         \
      }                                                        \
    }                                                          \
    return C;                                                  \
  }

/**
 * \brief Computes the matrix-scalar multiplication A * b.
 */
#define INSTANTIATE_MATSCAMUL_L(R, S, T)                        \
  template <int M, int N>                                       \
  mats<M, N, T> operator*(const mats<M, N, R>& A, const S& b) { \
    mats<M, N, T> C;                                            \
    for (int i = 0; i < M; i++)                                 \
      for (int j = 0; j < N; j++) C(i, j) = A(i, j) * b;        \
    return C;                                                   \
  }

/**
 * \brief Compute the scalar-matrix multiplication a * B.
 */
#define INSTANTIATE_MATSCAMUL_R(R, S, T)                        \
  template <int M, int N>                                       \
  mats<M, N, T> operator*(const R& b, const mats<M, N, S>& A) { \
    mats<M, N, T> C;                                            \
    for (int i = 0; i < M; i++)                                 \
      for (int j = 0; j < N; j++) C(i, j) = b * A(i, j);        \
    return C;                                                   \
  }

/**
 * \brief Computes the matrix-vector multiplication A * x.
 */
#define INSTANTIATE_MATVECMUL(R, S)                                         \
  template <int M, int N>                                                   \
  vecs<M, typename result_of<R, S>::type> operator*(const mats<M, N, R>& A, \
                                                    const vecs<N, S>& b) {  \
    typedef typename result_of<R, S>::type T;                               \
    vecs<M, T> c;                                                           \
    for (int i = 0; i < M; i++) {                                           \
      c(i) = 0;                                                             \
      for (int j = 0; j < N; j++) c(i) += A(i, j) * b(j);                   \
    }                                                                       \
    return c;                                                               \
  }

/**
 * \brief Unary + operator, for including A in expressions such as +A + B.
 */
#define INSTANTIATE_MATPLUS(R)                      \
  template <int M, int N>                           \
  mats<M, N, R> operator+(const mats<M, N, R>& A) { \
    return A;                                       \
  }

/**
 * \brief Unary - operator, for including A in expressions such as -A + B.
 */
#define INSTANTIATE_MATMINUS(R)                       \
  template <int M, int N>                             \
  mats<M, N, R> operator-(const mats<M, N, R>& A) {   \
    mats<M, N, R> B;                                  \
    for (int i = 0; i < M; i++)                       \
      for (int j = 0; j < N; j++) B(i, j) = -A(i, j); \
    return B;                                         \
  }

template <int M, int N, typename T>
mats<N, M, T> transpose(const mats<M, N, T>& A) {
  mats<N, M, T> At;
  for (int i = 0; i < M; i++)
    for (int j = 0; j < N; j++) At(j, i) = A(i, j);
  return At;
}

template <typename T>
T det(const mats<2, 2, T>& A) {
  return A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
}

template <typename T>
T det(const mats<3, 3, T>& A) {
  return A(0, 0) * (A(2, 2) * A(1, 1) - A(2, 1) * A(1, 2)) -
         A(1, 0) * (A(2, 2) * A(0, 1) - A(2, 1) * A(0, 2)) +
         A(2, 0) * (A(1, 2) * A(0, 1) - A(1, 1) * A(0, 2));
}

template <typename T>
mats<2, 2, T> inverse(const mats<2, 2, T>& A) {
  mats<2, 2, T> Ainv;
  const T id = 1.0 / det(A);
  Ainv(0, 0) = A(1, 1) * id;
  Ainv(0, 1) = -A(0, 1) * id;
  Ainv(1, 0) = -A(1, 0) * id;
  Ainv(1, 1) = A(0, 0) * id;
  return Ainv;
}

template <typename T>
mats<3, 3, T> inverse(const mats<3, 3, T>& M) {
  mats<3, 3, T> Minv;
  const T idetM = 1.0 / det(M);
  const T a1_1 = M(0, 0);
  const T a1_2 = M(0, 1);
  const T a1_3 = M(0, 2);
  const T a2_1 = M(1, 0);
  const T a2_2 = M(1, 1);
  const T a2_3 = M(1, 2);
  const T a3_1 = M(2, 0);
  const T a3_2 = M(2, 1);
  const T a3_3 = M(2, 2);
  Minv(0, 0) = (a2_2 * a3_3 - a2_3 * a3_2) * idetM;
  Minv(0, 1) = (a1_3 * a3_2 - a1_2 * a3_3) * idetM;
  Minv(0, 2) = (a1_2 * a2_3 - a1_3 * a2_2) * idetM;
  Minv(1, 0) = (a2_3 * a3_1 - a2_1 * a3_3) * idetM;
  Minv(1, 1) = (a1_1 * a3_3 - a1_3 * a3_1) * idetM;
  Minv(1, 2) = (a1_3 * a2_1 - a1_1 * a2_3) * idetM;
  Minv(2, 0) = (a2_1 * a3_2 - a2_2 * a3_1) * idetM;
  Minv(2, 1) = (a1_2 * a3_1 - a1_1 * a3_2) * idetM;
  Minv(2, 2) = (a1_1 * a2_2 - a1_2 * a2_1) * idetM;
  return Minv;
}

template <typename T>
mats<4, 4, T> inverse(const mats<4, 4, T>& M) {
  mats<4, 4, T> Minv;
  const T idetM = 1.0 / det(M);
  const T a1_1 = M(0, 0);
  const T a1_2 = M(0, 1);
  const T a1_3 = M(0, 2);
  const T a1_4 = M(0, 3);
  const T a2_1 = M(1, 0);
  const T a2_2 = M(1, 1);
  const T a2_3 = M(1, 2);
  const T a2_4 = M(1, 3);
  const T a3_1 = M(2, 0);
  const T a3_2 = M(2, 1);
  const T a3_3 = M(2, 2);
  const T a3_4 = M(2, 3);
  const T a4_1 = M(3, 0);
  const T a4_2 = M(3, 1);
  const T a4_3 = M(3, 2);
  const T a4_4 = M(3, 3);

  Minv(0, 0) = (a2_2 * a3_3 * a4_4 - a2_2 * a3_4 * a4_3 - a2_3 * a3_2 * a4_4 +
                a2_3 * a3_4 * a4_2 + a2_4 * a3_2 * a4_3 - a2_4 * a3_3 * a4_2) *
               idetM;
  Minv(0, 1) = (-a1_2 * a3_3 * a4_4 + a1_2 * a3_4 * a4_3 + a1_3 * a3_2 * a4_4 -
                a1_3 * a3_4 * a4_2 - a1_4 * a3_2 * a4_3 + a1_4 * a3_3 * a4_2) *
               idetM;
  Minv(0, 2) = (a1_2 * a2_3 * a4_4 - a1_2 * a2_4 * a4_3 - a1_3 * a2_2 * a4_4 +
                a1_3 * a2_4 * a4_2 + a1_4 * a2_2 * a4_3 - a1_4 * a2_3 * a4_2) *
               idetM;
  Minv(0, 3) = (-a1_2 * a2_3 * a3_4 + a1_2 * a2_4 * a3_3 + a1_3 * a2_2 * a3_4 -
                a1_3 * a2_4 * a3_2 - a1_4 * a2_2 * a3_3 + a1_4 * a2_3 * a3_2) *
               idetM;
  Minv(1, 0) = (-a2_1 * a3_3 * a4_4 + a2_1 * a3_4 * a4_3 + a2_3 * a3_1 * a4_4 -
                a2_3 * a3_4 * a4_1 - a2_4 * a3_1 * a4_3 + a2_4 * a3_3 * a4_1) *
               idetM;
  Minv(1, 1) = (a1_1 * a3_3 * a4_4 - a1_1 * a3_4 * a4_3 - a1_3 * a3_1 * a4_4 +
                a1_3 * a3_4 * a4_1 + a1_4 * a3_1 * a4_3 - a1_4 * a3_3 * a4_1) *
               idetM;
  Minv(1, 2) = (-a1_1 * a2_3 * a4_4 + a1_1 * a2_4 * a4_3 + a1_3 * a2_1 * a4_4 -
                a1_3 * a2_4 * a4_1 - a1_4 * a2_1 * a4_3 + a1_4 * a2_3 * a4_1) *
               idetM;
  Minv(1, 3) = (a1_1 * a2_3 * a3_4 - a1_1 * a2_4 * a3_3 - a1_3 * a2_1 * a3_4 +
                a1_3 * a2_4 * a3_1 + a1_4 * a2_1 * a3_3 - a1_4 * a2_3 * a3_1) *
               idetM;
  Minv(2, 0) = (a2_1 * a3_2 * a4_4 - a2_1 * a3_4 * a4_2 - a2_2 * a3_1 * a4_4 +
                a2_2 * a3_4 * a4_1 + a2_4 * a3_1 * a4_2 - a2_4 * a3_2 * a4_1) *
               idetM;
  Minv(2, 1) = (-a1_1 * a3_2 * a4_4 + a1_1 * a3_4 * a4_2 + a1_2 * a3_1 * a4_4 -
                a1_2 * a3_4 * a4_1 - a1_4 * a3_1 * a4_2 + a1_4 * a3_2 * a4_1) *
               idetM;
  Minv(2, 2) = (a1_1 * a2_2 * a4_4 - a1_1 * a2_4 * a4_2 - a1_2 * a2_1 * a4_4 +
                a1_2 * a2_4 * a4_1 + a1_4 * a2_1 * a4_2 - a1_4 * a2_2 * a4_1) *
               idetM;
  Minv(2, 3) = (-a1_1 * a2_2 * a3_4 + a1_1 * a2_4 * a3_2 + a1_2 * a2_1 * a3_4 -
                a1_2 * a2_4 * a3_1 - a1_4 * a2_1 * a3_2 + a1_4 * a2_2 * a3_1) *
               idetM;
  Minv(3, 0) = (-a2_1 * a3_2 * a4_3 + a2_1 * a3_3 * a4_2 + a2_2 * a3_1 * a4_3 -
                a2_2 * a3_3 * a4_1 - a2_3 * a3_1 * a4_2 + a2_3 * a3_2 * a4_1) *
               idetM;
  Minv(3, 1) = (a1_1 * a3_2 * a4_3 - a1_1 * a3_3 * a4_2 - a1_2 * a3_1 * a4_3 +
                a1_2 * a3_3 * a4_1 + a1_3 * a3_1 * a4_2 - a1_3 * a3_2 * a4_1) *
               idetM;
  Minv(3, 2) = (-a1_1 * a2_2 * a4_3 + a1_1 * a2_3 * a4_2 + a1_2 * a2_1 * a4_3 -
                a1_2 * a2_3 * a4_1 - a1_3 * a2_1 * a4_2 + a1_3 * a2_2 * a4_1) *
               idetM;
  Minv(3, 3) = (a1_1 * a2_2 * a3_3 - a1_1 * a2_3 * a3_2 - a1_2 * a2_1 * a3_3 +
                a1_2 * a2_3 * a3_1 + a1_3 * a2_1 * a3_2 - a1_3 * a2_2 * a3_1) *
               idetM;

  return Minv;
}

template <int N, typename T>
mats<N, N, T> inverse(const mats<N, N, T>& A) {
  matd<T> a(N, N);
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) a(i, j) = A(i, j);
  matd<T> ai = inverse(a);
  mats<N, N, T> Ai;
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) Ai(i, j) = ai(i, j);
  ASSERT(N <= 4);

  return Ai;
}

INSTANTIATE_MATADD(double, double, double)
INSTANTIATE_MATSUB(double, double, double)

INSTANTIATE_MATMUL(double, double, double)
INSTANTIATE_MATMUL(float, float, float)

INSTANTIATE_MATVECMUL(float, float)
INSTANTIATE_MATVECMUL(double, double)

INSTANTIATE_MATSCAMUL_L(double, double, double)
INSTANTIATE_MATSCAMUL_R(double, double, double)

INSTANTIATE_MATPLUS(double)
INSTANTIATE_MATMINUS(double)

}  // namespace vortex
