#include "linalg.h"
#include "log.h"
#include "mat.h"
#include "vec.h"

namespace vortex {

template <typename R, typename S>
matd<typename result_of<R, S>::type> operator*(const symd<R>& A,
                                               const symd<S>& B) {
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

template <int N, typename R, typename S>
mats<N, N, typename result_of<R, S>::type> operator*(const syms<N, R>& A,
                                                     const syms<N, S>& B) {
  typedef typename result_of<R, S>::type T;
  mats<N, N, T> C;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      T sum = 0;
      for (int k = 0; k < N; k++) sum += A(i, k) * B(k, j);
      C(i, j) = sum;
    }
  }
  return C;
}

template <typename R, typename S>
symd<typename result_of<R, S>::type> operator-(const symd<R>& A,
                                               const symd<S>& B) {
  typedef typename result_of<R, S>::type T;
  ASSERT(A.m() == B.m()) << "bad matrix sizes";
  ASSERT(A.n() == B.n()) << "bad matrix sizes";
  symd<T> C(A.m(), A.n());
  for (int i = 0; i < C.m(); i++) {
    for (int j = 0; j < C.n(); j++) {
      C(i, j) = A(i, j) - B(i, j);
    }
  }
  return C;
}

template <int N, typename R, typename S>
syms<N, typename result_of<R, S>::type> operator-(const syms<N, R>& A,
                                                  const syms<N, S>& B) {
  typedef typename result_of<R, S>::type T;
  syms<N, T> C;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      C(i, j) = A(i, j) - B(i, j);
    }
  }
  return C;
}

template <typename R, typename S>
matd<typename result_of<R, S>::type> operator-(const symd<R>& A,
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

template <int N, typename R, typename S>
mats<N, N, typename result_of<R, S>::type> operator-(const syms<N, R>& A,
                                                     const mats<N, N, S>& B) {
  typedef typename result_of<R, S>::type T;
  mats<N, N, T> C;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      C(i, j) = A(i, j) - B(i, j);
    }
  }
  return C;
}

template <typename R, typename S>
matd<typename result_of<R, S>::type> operator-(const matd<R>& A,
                                               const symd<S>& B) {
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

template <int N, typename R, typename S>
mats<N, N, typename result_of<R, S>::type> operator-(const mats<N, N, R>& A,
                                                     const syms<N, S>& B) {
  typedef typename result_of<R, S>::type T;
  mats<N, N, T> C;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      C(i, j) = A(i, j) - B(i, j);
    }
  }
  return C;
}

}  // namespace vortex
