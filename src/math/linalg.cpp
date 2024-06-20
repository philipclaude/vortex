#include "linalg.h"

#include <algorithm>
#include <vector>

#include "log.h"
#include "mat.h"
#include "mat.hpp"
#include "sym.h"
#include "vec.h"

namespace vortex {

template <typename T>
symd<T> interp(const std::vector<double>& alpha,
               const std::vector<symd<T>>& tensors) {
  ASSERT(alpha.size() == tensors.size());
  ASSERT(tensors.size() > 0);
  const int n = tensors[0].n();
  symd<T> m(n);
  m.zero();
  for (size_t k = 0; k < tensors.size(); k++)
    m = m + logm(tensors[k]) * alpha[k];
  return expm(m);
}

template <int N, typename T>
syms<N, T> interp(const std::vector<double>& alpha,
                  const std::vector<syms<N, T>>& tensors) {
  ASSERT(alpha.size() == tensors.size());
  ASSERT(tensors.size() > 0);
  syms<N, T> m;
  m.zero();
  for (size_t k = 0; k < tensors.size(); k++)
    m = m + logm(tensors[k]) * alpha[k];
  return expm(m);
}

template <typename T>
matd<T> transpose(const matd<T>& A) {
  matd<T> At(A.n(), A.m());
  for (int i = 0; i < A.m(); i++)
    for (int j = 0; j < A.n(); j++) At(j, i) = A(i, j);
  return At;
}

template <typename T>
matd<T> diag(const vecd<T>& d) {
  matd<T> A(d.m(), d.m());  // initializes to zero
  for (int i = 0; i < d.m(); i++) A(i, i) = d(i);
  return A;
}

template <int N, typename T>
mats<N, N, T> diag(const vecs<N, T>& d) {
  mats<N, N, T> A;
  for (int i = 0; i < N; i++) A(i, i) = d(i);
  return A;
}

template <typename T>
void decomposeLUP(const matd<T>& A, matd<T>& LU, std::vector<int>& P) {
  double tol = 1e-12;

  int m = A.m();
  int n = A.n();
  ASSERT(m == n) << fmt::format("matrix not square: {} x {}", m, n);

  // define unit permutation matrix
  P.resize(n + 1);
  for (int j = 0; j <= n; j++) P[j] = j;

  // copy A into B
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++) LU(i, j) = A(i, j);

  // initialize the rows that will be used to swap when pivoting
  vecd<T> row1(n);
  vecd<T> row2(n);
  for (int i = 0; i < m; i++) {
    // find the row with the maximum value
    T maxa = 0;
    int imax = i;

    for (int k = i; k < n; k++) {
      T absa = fabs(LU(k, i));
      if (absa > maxa) {
        maxa = absa;
        imax = k;
      }
    }
    ASSERT(maxa > tol) << "matrix is degenerate";

    if (imax != i) {
      // pivot P
      int j = P[i];
      P[i] = P[imax];
      P[imax] = j;

      // pivot rows of LU
      LU.get_row(i, row1);
      LU.get_row(imax, row2);
      LU.set_row(i, row2);
      LU.set_row(imax, row1);

      // increase the pivot counter for determinant
      P[n]++;
    }

    for (int j = i + 1; j < n; j++) {
      LU(j, i) /= LU(i, i);

      for (int k = i + 1; k < n; k++) LU(j, k) -= LU(j, i) * LU(i, k);
    }
  }
}

template <typename T>
void solveLUP(const matd<T>& LU, const std::vector<int>& P, const vecd<T>& b,
              vecd<T>& x) {
  int n = LU.n();
  ASSERT(x.m() == b.m());
  ASSERT(n == LU.m());

  // pivot the vector according to P
  for (int i = 0; i < n; i++) {
    x(i) = b(P[i]);
    for (int k = 0; k < i; k++) x(i) -= LU(i, k) * x(k);
  }

  // backwards substitution
  for (int i = n - 1; i >= 0; i--) {
    for (int k = i + 1; k < n; k++) x(i) -= LU(i, k) * x(k);
    x(i) /= LU(i, i);
  }
}

template <typename T>
void solveLUP(const matd<T>& A, const vecd<T>& b, vecd<T>& x) {
  matd<T> LU(A.m(), A.n());
  std::vector<int> P(A.n() + 1);
  decomposeLUP(A, LU, P);
  solveLUP(LU, P, b, x);
}

template <typename T>
void inverseLUP(const matd<T>& LU, const std::vector<int>& P,
                matd<double>& Ainv) {
  int n = LU.n();
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      Ainv(i, j) = (P[i] == j) ? 1.0 : 0.0;

      for (int k = 0; k < i; k++) Ainv(i, j) -= LU(i, k) * Ainv(k, j);
    }

    for (int i = n - 1; i >= 0; i--) {
      for (int k = i + 1; k < n; k++) Ainv(i, j) -= LU(i, k) * Ainv(k, j);
      Ainv(i, j) /= LU(i, i);
    }
  }
}

template <typename T>
void inverseLUP(const matd<T>& A, matd<T>& Ainv) {
  matd<T> LU(A.m(), A.n());
  std::vector<int> P(A.n() + 1);
  decomposeLUP(A, LU, P);
  inverseLUP(LU, P, Ainv);
}

template <typename T>
matd<T> inverse(const matd<T>& M) {
  ASSERT(M.m() == M.n());
  const T idetM = 1. / det(M);

  matd<T> Minv(M.m(), M.n());

  if (M.n() == 1) {
    Minv(0, 0) = idetM;
  } else if (M.n() == 2) {
    Minv(0, 0) = M(1, 1) * idetM;
    Minv(0, 1) = -M(0, 1) * idetM;
    Minv(1, 0) = -M(1, 0) * idetM;
    Minv(1, 1) = M(0, 0) * idetM;
  } else if (M.n() == 3) {
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
  } else if (M.n() == 4) {
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

    Minv(0, 0) =
        (a2_2 * a3_3 * a4_4 - a2_2 * a3_4 * a4_3 - a2_3 * a3_2 * a4_4 +
         a2_3 * a3_4 * a4_2 + a2_4 * a3_2 * a4_3 - a2_4 * a3_3 * a4_2) *
        idetM;
    Minv(0, 1) =
        (-a1_2 * a3_3 * a4_4 + a1_2 * a3_4 * a4_3 + a1_3 * a3_2 * a4_4 -
         a1_3 * a3_4 * a4_2 - a1_4 * a3_2 * a4_3 + a1_4 * a3_3 * a4_2) *
        idetM;
    Minv(0, 2) =
        (a1_2 * a2_3 * a4_4 - a1_2 * a2_4 * a4_3 - a1_3 * a2_2 * a4_4 +
         a1_3 * a2_4 * a4_2 + a1_4 * a2_2 * a4_3 - a1_4 * a2_3 * a4_2) *
        idetM;
    Minv(0, 3) =
        (-a1_2 * a2_3 * a3_4 + a1_2 * a2_4 * a3_3 + a1_3 * a2_2 * a3_4 -
         a1_3 * a2_4 * a3_2 - a1_4 * a2_2 * a3_3 + a1_4 * a2_3 * a3_2) *
        idetM;
    Minv(1, 0) =
        (-a2_1 * a3_3 * a4_4 + a2_1 * a3_4 * a4_3 + a2_3 * a3_1 * a4_4 -
         a2_3 * a3_4 * a4_1 - a2_4 * a3_1 * a4_3 + a2_4 * a3_3 * a4_1) *
        idetM;
    Minv(1, 1) =
        (a1_1 * a3_3 * a4_4 - a1_1 * a3_4 * a4_3 - a1_3 * a3_1 * a4_4 +
         a1_3 * a3_4 * a4_1 + a1_4 * a3_1 * a4_3 - a1_4 * a3_3 * a4_1) *
        idetM;
    Minv(1, 2) =
        (-a1_1 * a2_3 * a4_4 + a1_1 * a2_4 * a4_3 + a1_3 * a2_1 * a4_4 -
         a1_3 * a2_4 * a4_1 - a1_4 * a2_1 * a4_3 + a1_4 * a2_3 * a4_1) *
        idetM;
    Minv(1, 3) =
        (a1_1 * a2_3 * a3_4 - a1_1 * a2_4 * a3_3 - a1_3 * a2_1 * a3_4 +
         a1_3 * a2_4 * a3_1 + a1_4 * a2_1 * a3_3 - a1_4 * a2_3 * a3_1) *
        idetM;
    Minv(2, 0) =
        (a2_1 * a3_2 * a4_4 - a2_1 * a3_4 * a4_2 - a2_2 * a3_1 * a4_4 +
         a2_2 * a3_4 * a4_1 + a2_4 * a3_1 * a4_2 - a2_4 * a3_2 * a4_1) *
        idetM;
    Minv(2, 1) =
        (-a1_1 * a3_2 * a4_4 + a1_1 * a3_4 * a4_2 + a1_2 * a3_1 * a4_4 -
         a1_2 * a3_4 * a4_1 - a1_4 * a3_1 * a4_2 + a1_4 * a3_2 * a4_1) *
        idetM;
    Minv(2, 2) =
        (a1_1 * a2_2 * a4_4 - a1_1 * a2_4 * a4_2 - a1_2 * a2_1 * a4_4 +
         a1_2 * a2_4 * a4_1 + a1_4 * a2_1 * a4_2 - a1_4 * a2_2 * a4_1) *
        idetM;
    Minv(2, 3) =
        (-a1_1 * a2_2 * a3_4 + a1_1 * a2_4 * a3_2 + a1_2 * a2_1 * a3_4 -
         a1_2 * a2_4 * a3_1 - a1_4 * a2_1 * a3_2 + a1_4 * a2_2 * a3_1) *
        idetM;
    Minv(3, 0) =
        (-a2_1 * a3_2 * a4_3 + a2_1 * a3_3 * a4_2 + a2_2 * a3_1 * a4_3 -
         a2_2 * a3_3 * a4_1 - a2_3 * a3_1 * a4_2 + a2_3 * a3_2 * a4_1) *
        idetM;
    Minv(3, 1) =
        (a1_1 * a3_2 * a4_3 - a1_1 * a3_3 * a4_2 - a1_2 * a3_1 * a4_3 +
         a1_2 * a3_3 * a4_1 + a1_3 * a3_1 * a4_2 - a1_3 * a3_2 * a4_1) *
        idetM;
    Minv(3, 2) =
        (-a1_1 * a2_2 * a4_3 + a1_1 * a2_3 * a4_2 + a1_2 * a2_1 * a4_3 -
         a1_2 * a2_3 * a4_1 - a1_3 * a2_1 * a4_2 + a1_3 * a2_2 * a4_1) *
        idetM;
    Minv(3, 3) =
        (a1_1 * a2_2 * a3_3 - a1_1 * a2_3 * a3_2 - a1_2 * a2_1 * a3_3 +
         a1_2 * a2_3 * a3_1 + a1_3 * a2_1 * a3_2 - a1_3 * a2_2 * a3_1) *
        idetM;
  } else
    NOT_IMPLEMENTED;
  return Minv;
}

template <typename T>
symd<T> inverse(const symd<T>& M) {
  const T idetM = 1. / det(M);

  symd<T> Minv(M.m(), M.n());

  if (M.n() == 1) {
    Minv(0, 0) = idetM;
  } else if (M.n() == 2) {
    Minv(0, 0) = M(1, 1) * idetM;
    Minv(0, 1) = -M(0, 1) * idetM;
    Minv(1, 0) = -M(1, 0) * idetM;
    Minv(1, 1) = M(0, 0) * idetM;
  } else if (M.n() == 3) {
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
  } else if (M.n() == 4) {
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

    Minv(0, 0) =
        (a2_2 * a3_3 * a4_4 - a2_2 * a3_4 * a4_3 - a2_3 * a3_2 * a4_4 +
         a2_3 * a3_4 * a4_2 + a2_4 * a3_2 * a4_3 - a2_4 * a3_3 * a4_2) *
        idetM;
    Minv(0, 1) =
        (-a1_2 * a3_3 * a4_4 + a1_2 * a3_4 * a4_3 + a1_3 * a3_2 * a4_4 -
         a1_3 * a3_4 * a4_2 - a1_4 * a3_2 * a4_3 + a1_4 * a3_3 * a4_2) *
        idetM;
    Minv(0, 2) =
        (a1_2 * a2_3 * a4_4 - a1_2 * a2_4 * a4_3 - a1_3 * a2_2 * a4_4 +
         a1_3 * a2_4 * a4_2 + a1_4 * a2_2 * a4_3 - a1_4 * a2_3 * a4_2) *
        idetM;
    Minv(0, 3) =
        (-a1_2 * a2_3 * a3_4 + a1_2 * a2_4 * a3_3 + a1_3 * a2_2 * a3_4 -
         a1_3 * a2_4 * a3_2 - a1_4 * a2_2 * a3_3 + a1_4 * a2_3 * a3_2) *
        idetM;
    Minv(1, 0) =
        (-a2_1 * a3_3 * a4_4 + a2_1 * a3_4 * a4_3 + a2_3 * a3_1 * a4_4 -
         a2_3 * a3_4 * a4_1 - a2_4 * a3_1 * a4_3 + a2_4 * a3_3 * a4_1) *
        idetM;
    Minv(1, 1) =
        (a1_1 * a3_3 * a4_4 - a1_1 * a3_4 * a4_3 - a1_3 * a3_1 * a4_4 +
         a1_3 * a3_4 * a4_1 + a1_4 * a3_1 * a4_3 - a1_4 * a3_3 * a4_1) *
        idetM;
    Minv(1, 2) =
        (-a1_1 * a2_3 * a4_4 + a1_1 * a2_4 * a4_3 + a1_3 * a2_1 * a4_4 -
         a1_3 * a2_4 * a4_1 - a1_4 * a2_1 * a4_3 + a1_4 * a2_3 * a4_1) *
        idetM;
    Minv(1, 3) =
        (a1_1 * a2_3 * a3_4 - a1_1 * a2_4 * a3_3 - a1_3 * a2_1 * a3_4 +
         a1_3 * a2_4 * a3_1 + a1_4 * a2_1 * a3_3 - a1_4 * a2_3 * a3_1) *
        idetM;
    Minv(2, 0) =
        (a2_1 * a3_2 * a4_4 - a2_1 * a3_4 * a4_2 - a2_2 * a3_1 * a4_4 +
         a2_2 * a3_4 * a4_1 + a2_4 * a3_1 * a4_2 - a2_4 * a3_2 * a4_1) *
        idetM;
    Minv(2, 1) =
        (-a1_1 * a3_2 * a4_4 + a1_1 * a3_4 * a4_2 + a1_2 * a3_1 * a4_4 -
         a1_2 * a3_4 * a4_1 - a1_4 * a3_1 * a4_2 + a1_4 * a3_2 * a4_1) *
        idetM;
    Minv(2, 2) =
        (a1_1 * a2_2 * a4_4 - a1_1 * a2_4 * a4_2 - a1_2 * a2_1 * a4_4 +
         a1_2 * a2_4 * a4_1 + a1_4 * a2_1 * a4_2 - a1_4 * a2_2 * a4_1) *
        idetM;
    Minv(2, 3) =
        (-a1_1 * a2_2 * a3_4 + a1_1 * a2_4 * a3_2 + a1_2 * a2_1 * a3_4 -
         a1_2 * a2_4 * a3_1 - a1_4 * a2_1 * a3_2 + a1_4 * a2_2 * a3_1) *
        idetM;
    Minv(3, 0) =
        (-a2_1 * a3_2 * a4_3 + a2_1 * a3_3 * a4_2 + a2_2 * a3_1 * a4_3 -
         a2_2 * a3_3 * a4_1 - a2_3 * a3_1 * a4_2 + a2_3 * a3_2 * a4_1) *
        idetM;
    Minv(3, 1) =
        (a1_1 * a3_2 * a4_3 - a1_1 * a3_3 * a4_2 - a1_2 * a3_1 * a4_3 +
         a1_2 * a3_3 * a4_1 + a1_3 * a3_1 * a4_2 - a1_3 * a3_2 * a4_1) *
        idetM;
    Minv(3, 2) =
        (-a1_1 * a2_2 * a4_3 + a1_1 * a2_3 * a4_2 + a1_2 * a2_1 * a4_3 -
         a1_2 * a2_3 * a4_1 - a1_3 * a2_1 * a4_2 + a1_3 * a2_2 * a4_1) *
        idetM;
    Minv(3, 3) =
        (a1_1 * a2_2 * a3_3 - a1_1 * a2_3 * a3_2 - a1_2 * a2_1 * a3_3 +
         a1_2 * a2_3 * a3_1 + a1_3 * a2_1 * a3_2 - a1_3 * a2_2 * a3_1) *
        idetM;
  } else
    NOT_IMPLEMENTED;
  return Minv;
}

template <typename T>
std::pair<vecd<T>, matd<T>> eig(const symd<T>& m) {
  return m.eig();
}

template <typename T>
void eig(const symd<T>& m, vecd<T>& L, matd<T>& Q) {
  std::pair<vecd<T>, matd<T>> decomp = m.eig();
  L.set(decomp.first);
  Q.set(decomp.second);
}

template <int N, typename T>
std::pair<vecs<N, T>, mats<N, N, T>> eig(const syms<N, T>& m) {
  return m.eig();
}

template <int N, typename T>
void eig(const syms<N, T>& m, vecs<N, T>& L, mats<N, N, T>& Q) {
  std::pair<vecs<N, T>, mats<N, N, T>> decomp = m.eig();
  for (int i = 0; i < N; i++) {
    L[i] = decomp.first[i];
    for (int j = 0; j < N; j++) Q(i, j) = decomp.second(i, j);
  }
}

template <typename type>
type det(const matd<type>& X) {
  ASSERT(X.m() == X.n());
  if (X.m() == 1) return X(0, 0);
  if (X.m() == 2) return X(1, 1) * X(0, 0) - X(0, 1) * X(1, 0);
  if (X.m() == 3)
    return X(0, 0) * (X(2, 2) * X(1, 1) - X(2, 1) * X(1, 2)) -
           X(1, 0) * (X(2, 2) * X(0, 1) - X(2, 1) * X(0, 2)) +
           X(2, 0) * (X(1, 2) * X(0, 1) - X(1, 1) * X(0, 2));
  const type X1_1 = X(0, 0), X1_2 = X(0, 1), X1_3 = X(0, 2), X1_4 = X(0, 3);
  const type X2_1 = X(1, 0), X2_2 = X(1, 1), X2_3 = X(1, 2), X2_4 = X(1, 3);
  const type X3_1 = X(2, 0), X3_2 = X(2, 1), X3_3 = X(2, 2), X3_4 = X(2, 3);
  const type X4_1 = X(3, 0), X4_2 = X(3, 1), X4_3 = X(3, 2), X4_4 = X(3, 3);
  if (X.m() == 4) {
    return X1_1 * X2_2 * X3_3 * X4_4 - X1_1 * X2_2 * X3_4 * X4_3 -
           X1_1 * X2_3 * X3_2 * X4_4 + X1_1 * X2_3 * X3_4 * X4_2 +
           X1_1 * X2_4 * X3_2 * X4_3 - X1_1 * X2_4 * X3_3 * X4_2 -
           X1_2 * X2_1 * X3_3 * X4_4 + X1_2 * X2_1 * X3_4 * X4_3 +
           X1_2 * X2_3 * X3_1 * X4_4 - X1_2 * X2_3 * X3_4 * X4_1 -
           X1_2 * X2_4 * X3_1 * X4_3 + X1_2 * X2_4 * X3_3 * X4_1 +
           X1_3 * X2_1 * X3_2 * X4_4 - X1_3 * X2_1 * X3_4 * X4_2 -
           X1_3 * X2_2 * X3_1 * X4_4 + X1_3 * X2_2 * X3_4 * X4_1 +
           X1_3 * X2_4 * X3_1 * X4_2 - X1_3 * X2_4 * X3_2 * X4_1 -
           X1_4 * X2_1 * X3_2 * X4_3 + X1_4 * X2_1 * X3_3 * X4_2 +
           X1_4 * X2_2 * X3_1 * X4_3 - X1_4 * X2_2 * X3_3 * X4_1 -
           X1_4 * X2_3 * X3_1 * X4_2 + X1_4 * X2_3 * X3_2 * X4_1;
  }
  const type X1_5 = X(0, 4), X2_5 = X(1, 4), X3_5 = X(2, 4), X4_5 = X(3, 4);
  const type X5_1 = X(4, 0), X5_2 = X(4, 1), X5_3 = X(4, 2), X5_4 = X(4, 3),
             X5_5 = X(4, 4);
  if (X.m() == 5) {
    return X1_1 * X2_2 * X3_3 * X4_4 * X5_5 - X1_1 * X2_2 * X3_3 * X4_5 * X5_4 -
           X1_1 * X2_2 * X3_4 * X4_3 * X5_5 + X1_1 * X2_2 * X3_4 * X4_5 * X5_3 +
           X1_1 * X2_2 * X3_5 * X4_3 * X5_4 - X1_1 * X2_2 * X3_5 * X4_4 * X5_3 -
           X1_1 * X2_3 * X3_2 * X4_4 * X5_5 + X1_1 * X2_3 * X3_2 * X4_5 * X5_4 +
           X1_1 * X2_3 * X3_4 * X4_2 * X5_5 - X1_1 * X2_3 * X3_4 * X4_5 * X5_2 -
           X1_1 * X2_3 * X3_5 * X4_2 * X5_4 + X1_1 * X2_3 * X3_5 * X4_4 * X5_2 +
           X1_1 * X2_4 * X3_2 * X4_3 * X5_5 - X1_1 * X2_4 * X3_2 * X4_5 * X5_3 -
           X1_1 * X2_4 * X3_3 * X4_2 * X5_5 + X1_1 * X2_4 * X3_3 * X4_5 * X5_2 +
           X1_1 * X2_4 * X3_5 * X4_2 * X5_3 - X1_1 * X2_4 * X3_5 * X4_3 * X5_2 -
           X1_1 * X2_5 * X3_2 * X4_3 * X5_4 + X1_1 * X2_5 * X3_2 * X4_4 * X5_3 +
           X1_1 * X2_5 * X3_3 * X4_2 * X5_4 - X1_1 * X2_5 * X3_3 * X4_4 * X5_2 -
           X1_1 * X2_5 * X3_4 * X4_2 * X5_3 + X1_1 * X2_5 * X3_4 * X4_3 * X5_2 -
           X1_2 * X2_1 * X3_3 * X4_4 * X5_5 + X1_2 * X2_1 * X3_3 * X4_5 * X5_4 +
           X1_2 * X2_1 * X3_4 * X4_3 * X5_5 - X1_2 * X2_1 * X3_4 * X4_5 * X5_3 -
           X1_2 * X2_1 * X3_5 * X4_3 * X5_4 + X1_2 * X2_1 * X3_5 * X4_4 * X5_3 +
           X1_2 * X2_3 * X3_1 * X4_4 * X5_5 - X1_2 * X2_3 * X3_1 * X4_5 * X5_4 -
           X1_2 * X2_3 * X3_4 * X4_1 * X5_5 + X1_2 * X2_3 * X3_4 * X4_5 * X5_1 +
           X1_2 * X2_3 * X3_5 * X4_1 * X5_4 - X1_2 * X2_3 * X3_5 * X4_4 * X5_1 -
           X1_2 * X2_4 * X3_1 * X4_3 * X5_5 + X1_2 * X2_4 * X3_1 * X4_5 * X5_3 +
           X1_2 * X2_4 * X3_3 * X4_1 * X5_5 - X1_2 * X2_4 * X3_3 * X4_5 * X5_1 -
           X1_2 * X2_4 * X3_5 * X4_1 * X5_3 + X1_2 * X2_4 * X3_5 * X4_3 * X5_1 +
           X1_2 * X2_5 * X3_1 * X4_3 * X5_4 - X1_2 * X2_5 * X3_1 * X4_4 * X5_3 -
           X1_2 * X2_5 * X3_3 * X4_1 * X5_4 + X1_2 * X2_5 * X3_3 * X4_4 * X5_1 +
           X1_2 * X2_5 * X3_4 * X4_1 * X5_3 - X1_2 * X2_5 * X3_4 * X4_3 * X5_1 +
           X1_3 * X2_1 * X3_2 * X4_4 * X5_5 - X1_3 * X2_1 * X3_2 * X4_5 * X5_4 -
           X1_3 * X2_1 * X3_4 * X4_2 * X5_5 + X1_3 * X2_1 * X3_4 * X4_5 * X5_2 +
           X1_3 * X2_1 * X3_5 * X4_2 * X5_4 - X1_3 * X2_1 * X3_5 * X4_4 * X5_2 -
           X1_3 * X2_2 * X3_1 * X4_4 * X5_5 + X1_3 * X2_2 * X3_1 * X4_5 * X5_4 +
           X1_3 * X2_2 * X3_4 * X4_1 * X5_5 - X1_3 * X2_2 * X3_4 * X4_5 * X5_1 -
           X1_3 * X2_2 * X3_5 * X4_1 * X5_4 + X1_3 * X2_2 * X3_5 * X4_4 * X5_1 +
           X1_3 * X2_4 * X3_1 * X4_2 * X5_5 - X1_3 * X2_4 * X3_1 * X4_5 * X5_2 -
           X1_3 * X2_4 * X3_2 * X4_1 * X5_5 + X1_3 * X2_4 * X3_2 * X4_5 * X5_1 +
           X1_3 * X2_4 * X3_5 * X4_1 * X5_2 - X1_3 * X2_4 * X3_5 * X4_2 * X5_1 -
           X1_3 * X2_5 * X3_1 * X4_2 * X5_4 + X1_3 * X2_5 * X3_1 * X4_4 * X5_2 +
           X1_3 * X2_5 * X3_2 * X4_1 * X5_4 - X1_3 * X2_5 * X3_2 * X4_4 * X5_1 -
           X1_3 * X2_5 * X3_4 * X4_1 * X5_2 + X1_3 * X2_5 * X3_4 * X4_2 * X5_1 -
           X1_4 * X2_1 * X3_2 * X4_3 * X5_5 + X1_4 * X2_1 * X3_2 * X4_5 * X5_3 +
           X1_4 * X2_1 * X3_3 * X4_2 * X5_5 - X1_4 * X2_1 * X3_3 * X4_5 * X5_2 -
           X1_4 * X2_1 * X3_5 * X4_2 * X5_3 + X1_4 * X2_1 * X3_5 * X4_3 * X5_2 +
           X1_4 * X2_2 * X3_1 * X4_3 * X5_5 - X1_4 * X2_2 * X3_1 * X4_5 * X5_3 -
           X1_4 * X2_2 * X3_3 * X4_1 * X5_5 + X1_4 * X2_2 * X3_3 * X4_5 * X5_1 +
           X1_4 * X2_2 * X3_5 * X4_1 * X5_3 - X1_4 * X2_2 * X3_5 * X4_3 * X5_1 -
           X1_4 * X2_3 * X3_1 * X4_2 * X5_5 + X1_4 * X2_3 * X3_1 * X4_5 * X5_2 +
           X1_4 * X2_3 * X3_2 * X4_1 * X5_5 - X1_4 * X2_3 * X3_2 * X4_5 * X5_1 -
           X1_4 * X2_3 * X3_5 * X4_1 * X5_2 + X1_4 * X2_3 * X3_5 * X4_2 * X5_1 +
           X1_4 * X2_5 * X3_1 * X4_2 * X5_3 - X1_4 * X2_5 * X3_1 * X4_3 * X5_2 -
           X1_4 * X2_5 * X3_2 * X4_1 * X5_3 + X1_4 * X2_5 * X3_2 * X4_3 * X5_1 +
           X1_4 * X2_5 * X3_3 * X4_1 * X5_2 - X1_4 * X2_5 * X3_3 * X4_2 * X5_1 +
           X1_5 * X2_1 * X3_2 * X4_3 * X5_4 - X1_5 * X2_1 * X3_2 * X4_4 * X5_3 -
           X1_5 * X2_1 * X3_3 * X4_2 * X5_4 + X1_5 * X2_1 * X3_3 * X4_4 * X5_2 +
           X1_5 * X2_1 * X3_4 * X4_2 * X5_3 - X1_5 * X2_1 * X3_4 * X4_3 * X5_2 -
           X1_5 * X2_2 * X3_1 * X4_3 * X5_4 + X1_5 * X2_2 * X3_1 * X4_4 * X5_3 +
           X1_5 * X2_2 * X3_3 * X4_1 * X5_4 - X1_5 * X2_2 * X3_3 * X4_4 * X5_1 -
           X1_5 * X2_2 * X3_4 * X4_1 * X5_3 + X1_5 * X2_2 * X3_4 * X4_3 * X5_1 +
           X1_5 * X2_3 * X3_1 * X4_2 * X5_4 - X1_5 * X2_3 * X3_1 * X4_4 * X5_2 -
           X1_5 * X2_3 * X3_2 * X4_1 * X5_4 + X1_5 * X2_3 * X3_2 * X4_4 * X5_1 +
           X1_5 * X2_3 * X3_4 * X4_1 * X5_2 - X1_5 * X2_3 * X3_4 * X4_2 * X5_1 -
           X1_5 * X2_4 * X3_1 * X4_2 * X5_3 + X1_5 * X2_4 * X3_1 * X4_3 * X5_2 +
           X1_5 * X2_4 * X3_2 * X4_1 * X5_3 - X1_5 * X2_4 * X3_2 * X4_3 * X5_1 -
           X1_5 * X2_4 * X3_3 * X4_1 * X5_2 + X1_5 * X2_4 * X3_3 * X4_2 * X5_1;
  }
  const type X1_6 = X(0, 5), X2_6 = X(1, 5), X3_6 = X(2, 5), X4_6 = X(3, 5),
             X5_6 = X(4, 5);
  const type X6_1 = X(5, 0), X6_2 = X(5, 1), X6_3 = X(5, 2), X6_4 = X(5, 3),
             X6_5 = X(5, 4), X6_6 = X(5, 5);
  if (X.n() == 6) {
    return X1_1 * X2_2 * X3_3 * X4_4 * X5_5 * X6_6 -
           X1_1 * X2_2 * X3_3 * X4_4 * X5_6 * X6_5 -
           X1_1 * X2_2 * X3_3 * X4_5 * X5_4 * X6_6 +
           X1_1 * X2_2 * X3_3 * X4_5 * X5_6 * X6_4 +
           X1_1 * X2_2 * X3_3 * X4_6 * X5_4 * X6_5 -
           X1_1 * X2_2 * X3_3 * X4_6 * X5_5 * X6_4 -
           X1_1 * X2_2 * X3_4 * X4_3 * X5_5 * X6_6 +
           X1_1 * X2_2 * X3_4 * X4_3 * X5_6 * X6_5 +
           X1_1 * X2_2 * X3_4 * X4_5 * X5_3 * X6_6 -
           X1_1 * X2_2 * X3_4 * X4_5 * X5_6 * X6_3 -
           X1_1 * X2_2 * X3_4 * X4_6 * X5_3 * X6_5 +
           X1_1 * X2_2 * X3_4 * X4_6 * X5_5 * X6_3 +
           X1_1 * X2_2 * X3_5 * X4_3 * X5_4 * X6_6 -
           X1_1 * X2_2 * X3_5 * X4_3 * X5_6 * X6_4 -
           X1_1 * X2_2 * X3_5 * X4_4 * X5_3 * X6_6 +
           X1_1 * X2_2 * X3_5 * X4_4 * X5_6 * X6_3 +
           X1_1 * X2_2 * X3_5 * X4_6 * X5_3 * X6_4 -
           X1_1 * X2_2 * X3_5 * X4_6 * X5_4 * X6_3 -
           X1_1 * X2_2 * X3_6 * X4_3 * X5_4 * X6_5 +
           X1_1 * X2_2 * X3_6 * X4_3 * X5_5 * X6_4 +
           X1_1 * X2_2 * X3_6 * X4_4 * X5_3 * X6_5 -
           X1_1 * X2_2 * X3_6 * X4_4 * X5_5 * X6_3 -
           X1_1 * X2_2 * X3_6 * X4_5 * X5_3 * X6_4 +
           X1_1 * X2_2 * X3_6 * X4_5 * X5_4 * X6_3 -
           X1_1 * X2_3 * X3_2 * X4_4 * X5_5 * X6_6 +
           X1_1 * X2_3 * X3_2 * X4_4 * X5_6 * X6_5 +
           X1_1 * X2_3 * X3_2 * X4_5 * X5_4 * X6_6 -
           X1_1 * X2_3 * X3_2 * X4_5 * X5_6 * X6_4 -
           X1_1 * X2_3 * X3_2 * X4_6 * X5_4 * X6_5 +
           X1_1 * X2_3 * X3_2 * X4_6 * X5_5 * X6_4 +
           X1_1 * X2_3 * X3_4 * X4_2 * X5_5 * X6_6 -
           X1_1 * X2_3 * X3_4 * X4_2 * X5_6 * X6_5 -
           X1_1 * X2_3 * X3_4 * X4_5 * X5_2 * X6_6 +
           X1_1 * X2_3 * X3_4 * X4_5 * X5_6 * X6_2 +
           X1_1 * X2_3 * X3_4 * X4_6 * X5_2 * X6_5 -
           X1_1 * X2_3 * X3_4 * X4_6 * X5_5 * X6_2 -
           X1_1 * X2_3 * X3_5 * X4_2 * X5_4 * X6_6 +
           X1_1 * X2_3 * X3_5 * X4_2 * X5_6 * X6_4 +
           X1_1 * X2_3 * X3_5 * X4_4 * X5_2 * X6_6 -
           X1_1 * X2_3 * X3_5 * X4_4 * X5_6 * X6_2 -
           X1_1 * X2_3 * X3_5 * X4_6 * X5_2 * X6_4 +
           X1_1 * X2_3 * X3_5 * X4_6 * X5_4 * X6_2 +
           X1_1 * X2_3 * X3_6 * X4_2 * X5_4 * X6_5 -
           X1_1 * X2_3 * X3_6 * X4_2 * X5_5 * X6_4 -
           X1_1 * X2_3 * X3_6 * X4_4 * X5_2 * X6_5 +
           X1_1 * X2_3 * X3_6 * X4_4 * X5_5 * X6_2 +
           X1_1 * X2_3 * X3_6 * X4_5 * X5_2 * X6_4 -
           X1_1 * X2_3 * X3_6 * X4_5 * X5_4 * X6_2 +
           X1_1 * X2_4 * X3_2 * X4_3 * X5_5 * X6_6 -
           X1_1 * X2_4 * X3_2 * X4_3 * X5_6 * X6_5 -
           X1_1 * X2_4 * X3_2 * X4_5 * X5_3 * X6_6 +
           X1_1 * X2_4 * X3_2 * X4_5 * X5_6 * X6_3 +
           X1_1 * X2_4 * X3_2 * X4_6 * X5_3 * X6_5 -
           X1_1 * X2_4 * X3_2 * X4_6 * X5_5 * X6_3 -
           X1_1 * X2_4 * X3_3 * X4_2 * X5_5 * X6_6 +
           X1_1 * X2_4 * X3_3 * X4_2 * X5_6 * X6_5 +
           X1_1 * X2_4 * X3_3 * X4_5 * X5_2 * X6_6 -
           X1_1 * X2_4 * X3_3 * X4_5 * X5_6 * X6_2 -
           X1_1 * X2_4 * X3_3 * X4_6 * X5_2 * X6_5 +
           X1_1 * X2_4 * X3_3 * X4_6 * X5_5 * X6_2 +
           X1_1 * X2_4 * X3_5 * X4_2 * X5_3 * X6_6 -
           X1_1 * X2_4 * X3_5 * X4_2 * X5_6 * X6_3 -
           X1_1 * X2_4 * X3_5 * X4_3 * X5_2 * X6_6 +
           X1_1 * X2_4 * X3_5 * X4_3 * X5_6 * X6_2 +
           X1_1 * X2_4 * X3_5 * X4_6 * X5_2 * X6_3 -
           X1_1 * X2_4 * X3_5 * X4_6 * X5_3 * X6_2 -
           X1_1 * X2_4 * X3_6 * X4_2 * X5_3 * X6_5 +
           X1_1 * X2_4 * X3_6 * X4_2 * X5_5 * X6_3 +
           X1_1 * X2_4 * X3_6 * X4_3 * X5_2 * X6_5 -
           X1_1 * X2_4 * X3_6 * X4_3 * X5_5 * X6_2 -
           X1_1 * X2_4 * X3_6 * X4_5 * X5_2 * X6_3 +
           X1_1 * X2_4 * X3_6 * X4_5 * X5_3 * X6_2 -
           X1_1 * X2_5 * X3_2 * X4_3 * X5_4 * X6_6 +
           X1_1 * X2_5 * X3_2 * X4_3 * X5_6 * X6_4 +
           X1_1 * X2_5 * X3_2 * X4_4 * X5_3 * X6_6 -
           X1_1 * X2_5 * X3_2 * X4_4 * X5_6 * X6_3 -
           X1_1 * X2_5 * X3_2 * X4_6 * X5_3 * X6_4 +
           X1_1 * X2_5 * X3_2 * X4_6 * X5_4 * X6_3 +
           X1_1 * X2_5 * X3_3 * X4_2 * X5_4 * X6_6 -
           X1_1 * X2_5 * X3_3 * X4_2 * X5_6 * X6_4 -
           X1_1 * X2_5 * X3_3 * X4_4 * X5_2 * X6_6 +
           X1_1 * X2_5 * X3_3 * X4_4 * X5_6 * X6_2 +
           X1_1 * X2_5 * X3_3 * X4_6 * X5_2 * X6_4 -
           X1_1 * X2_5 * X3_3 * X4_6 * X5_4 * X6_2 -
           X1_1 * X2_5 * X3_4 * X4_2 * X5_3 * X6_6 +
           X1_1 * X2_5 * X3_4 * X4_2 * X5_6 * X6_3 +
           X1_1 * X2_5 * X3_4 * X4_3 * X5_2 * X6_6 -
           X1_1 * X2_5 * X3_4 * X4_3 * X5_6 * X6_2 -
           X1_1 * X2_5 * X3_4 * X4_6 * X5_2 * X6_3 +
           X1_1 * X2_5 * X3_4 * X4_6 * X5_3 * X6_2 +
           X1_1 * X2_5 * X3_6 * X4_2 * X5_3 * X6_4 -
           X1_1 * X2_5 * X3_6 * X4_2 * X5_4 * X6_3 -
           X1_1 * X2_5 * X3_6 * X4_3 * X5_2 * X6_4 +
           X1_1 * X2_5 * X3_6 * X4_3 * X5_4 * X6_2 +
           X1_1 * X2_5 * X3_6 * X4_4 * X5_2 * X6_3 -
           X1_1 * X2_5 * X3_6 * X4_4 * X5_3 * X6_2 +
           X1_1 * X2_6 * X3_2 * X4_3 * X5_4 * X6_5 -
           X1_1 * X2_6 * X3_2 * X4_3 * X5_5 * X6_4 -
           X1_1 * X2_6 * X3_2 * X4_4 * X5_3 * X6_5 +
           X1_1 * X2_6 * X3_2 * X4_4 * X5_5 * X6_3 +
           X1_1 * X2_6 * X3_2 * X4_5 * X5_3 * X6_4 -
           X1_1 * X2_6 * X3_2 * X4_5 * X5_4 * X6_3 -
           X1_1 * X2_6 * X3_3 * X4_2 * X5_4 * X6_5 +
           X1_1 * X2_6 * X3_3 * X4_2 * X5_5 * X6_4 +
           X1_1 * X2_6 * X3_3 * X4_4 * X5_2 * X6_5 -
           X1_1 * X2_6 * X3_3 * X4_4 * X5_5 * X6_2 -
           X1_1 * X2_6 * X3_3 * X4_5 * X5_2 * X6_4 +
           X1_1 * X2_6 * X3_3 * X4_5 * X5_4 * X6_2 +
           X1_1 * X2_6 * X3_4 * X4_2 * X5_3 * X6_5 -
           X1_1 * X2_6 * X3_4 * X4_2 * X5_5 * X6_3 -
           X1_1 * X2_6 * X3_4 * X4_3 * X5_2 * X6_5 +
           X1_1 * X2_6 * X3_4 * X4_3 * X5_5 * X6_2 +
           X1_1 * X2_6 * X3_4 * X4_5 * X5_2 * X6_3 -
           X1_1 * X2_6 * X3_4 * X4_5 * X5_3 * X6_2 -
           X1_1 * X2_6 * X3_5 * X4_2 * X5_3 * X6_4 +
           X1_1 * X2_6 * X3_5 * X4_2 * X5_4 * X6_3 +
           X1_1 * X2_6 * X3_5 * X4_3 * X5_2 * X6_4 -
           X1_1 * X2_6 * X3_5 * X4_3 * X5_4 * X6_2 -
           X1_1 * X2_6 * X3_5 * X4_4 * X5_2 * X6_3 +
           X1_1 * X2_6 * X3_5 * X4_4 * X5_3 * X6_2 -
           X1_2 * X2_1 * X3_3 * X4_4 * X5_5 * X6_6 +
           X1_2 * X2_1 * X3_3 * X4_4 * X5_6 * X6_5 +
           X1_2 * X2_1 * X3_3 * X4_5 * X5_4 * X6_6 -
           X1_2 * X2_1 * X3_3 * X4_5 * X5_6 * X6_4 -
           X1_2 * X2_1 * X3_3 * X4_6 * X5_4 * X6_5 +
           X1_2 * X2_1 * X3_3 * X4_6 * X5_5 * X6_4 +
           X1_2 * X2_1 * X3_4 * X4_3 * X5_5 * X6_6 -
           X1_2 * X2_1 * X3_4 * X4_3 * X5_6 * X6_5 -
           X1_2 * X2_1 * X3_4 * X4_5 * X5_3 * X6_6 +
           X1_2 * X2_1 * X3_4 * X4_5 * X5_6 * X6_3 +
           X1_2 * X2_1 * X3_4 * X4_6 * X5_3 * X6_5 -
           X1_2 * X2_1 * X3_4 * X4_6 * X5_5 * X6_3 -
           X1_2 * X2_1 * X3_5 * X4_3 * X5_4 * X6_6 +
           X1_2 * X2_1 * X3_5 * X4_3 * X5_6 * X6_4 +
           X1_2 * X2_1 * X3_5 * X4_4 * X5_3 * X6_6 -
           X1_2 * X2_1 * X3_5 * X4_4 * X5_6 * X6_3 -
           X1_2 * X2_1 * X3_5 * X4_6 * X5_3 * X6_4 +
           X1_2 * X2_1 * X3_5 * X4_6 * X5_4 * X6_3 +
           X1_2 * X2_1 * X3_6 * X4_3 * X5_4 * X6_5 -
           X1_2 * X2_1 * X3_6 * X4_3 * X5_5 * X6_4 -
           X1_2 * X2_1 * X3_6 * X4_4 * X5_3 * X6_5 +
           X1_2 * X2_1 * X3_6 * X4_4 * X5_5 * X6_3 +
           X1_2 * X2_1 * X3_6 * X4_5 * X5_3 * X6_4 -
           X1_2 * X2_1 * X3_6 * X4_5 * X5_4 * X6_3 +
           X1_2 * X2_3 * X3_1 * X4_4 * X5_5 * X6_6 -
           X1_2 * X2_3 * X3_1 * X4_4 * X5_6 * X6_5 -
           X1_2 * X2_3 * X3_1 * X4_5 * X5_4 * X6_6 +
           X1_2 * X2_3 * X3_1 * X4_5 * X5_6 * X6_4 +
           X1_2 * X2_3 * X3_1 * X4_6 * X5_4 * X6_5 -
           X1_2 * X2_3 * X3_1 * X4_6 * X5_5 * X6_4 -
           X1_2 * X2_3 * X3_4 * X4_1 * X5_5 * X6_6 +
           X1_2 * X2_3 * X3_4 * X4_1 * X5_6 * X6_5 +
           X1_2 * X2_3 * X3_4 * X4_5 * X5_1 * X6_6 -
           X1_2 * X2_3 * X3_4 * X4_5 * X5_6 * X6_1 -
           X1_2 * X2_3 * X3_4 * X4_6 * X5_1 * X6_5 +
           X1_2 * X2_3 * X3_4 * X4_6 * X5_5 * X6_1 +
           X1_2 * X2_3 * X3_5 * X4_1 * X5_4 * X6_6 -
           X1_2 * X2_3 * X3_5 * X4_1 * X5_6 * X6_4 -
           X1_2 * X2_3 * X3_5 * X4_4 * X5_1 * X6_6 +
           X1_2 * X2_3 * X3_5 * X4_4 * X5_6 * X6_1 +
           X1_2 * X2_3 * X3_5 * X4_6 * X5_1 * X6_4 -
           X1_2 * X2_3 * X3_5 * X4_6 * X5_4 * X6_1 -
           X1_2 * X2_3 * X3_6 * X4_1 * X5_4 * X6_5 +
           X1_2 * X2_3 * X3_6 * X4_1 * X5_5 * X6_4 +
           X1_2 * X2_3 * X3_6 * X4_4 * X5_1 * X6_5 -
           X1_2 * X2_3 * X3_6 * X4_4 * X5_5 * X6_1 -
           X1_2 * X2_3 * X3_6 * X4_5 * X5_1 * X6_4 +
           X1_2 * X2_3 * X3_6 * X4_5 * X5_4 * X6_1 -
           X1_2 * X2_4 * X3_1 * X4_3 * X5_5 * X6_6 +
           X1_2 * X2_4 * X3_1 * X4_3 * X5_6 * X6_5 +
           X1_2 * X2_4 * X3_1 * X4_5 * X5_3 * X6_6 -
           X1_2 * X2_4 * X3_1 * X4_5 * X5_6 * X6_3 -
           X1_2 * X2_4 * X3_1 * X4_6 * X5_3 * X6_5 +
           X1_2 * X2_4 * X3_1 * X4_6 * X5_5 * X6_3 +
           X1_2 * X2_4 * X3_3 * X4_1 * X5_5 * X6_6 -
           X1_2 * X2_4 * X3_3 * X4_1 * X5_6 * X6_5 -
           X1_2 * X2_4 * X3_3 * X4_5 * X5_1 * X6_6 +
           X1_2 * X2_4 * X3_3 * X4_5 * X5_6 * X6_1 +
           X1_2 * X2_4 * X3_3 * X4_6 * X5_1 * X6_5 -
           X1_2 * X2_4 * X3_3 * X4_6 * X5_5 * X6_1 -
           X1_2 * X2_4 * X3_5 * X4_1 * X5_3 * X6_6 +
           X1_2 * X2_4 * X3_5 * X4_1 * X5_6 * X6_3 +
           X1_2 * X2_4 * X3_5 * X4_3 * X5_1 * X6_6 -
           X1_2 * X2_4 * X3_5 * X4_3 * X5_6 * X6_1 -
           X1_2 * X2_4 * X3_5 * X4_6 * X5_1 * X6_3 +
           X1_2 * X2_4 * X3_5 * X4_6 * X5_3 * X6_1 +
           X1_2 * X2_4 * X3_6 * X4_1 * X5_3 * X6_5 -
           X1_2 * X2_4 * X3_6 * X4_1 * X5_5 * X6_3 -
           X1_2 * X2_4 * X3_6 * X4_3 * X5_1 * X6_5 +
           X1_2 * X2_4 * X3_6 * X4_3 * X5_5 * X6_1 +
           X1_2 * X2_4 * X3_6 * X4_5 * X5_1 * X6_3 -
           X1_2 * X2_4 * X3_6 * X4_5 * X5_3 * X6_1 +
           X1_2 * X2_5 * X3_1 * X4_3 * X5_4 * X6_6 -
           X1_2 * X2_5 * X3_1 * X4_3 * X5_6 * X6_4 -
           X1_2 * X2_5 * X3_1 * X4_4 * X5_3 * X6_6 +
           X1_2 * X2_5 * X3_1 * X4_4 * X5_6 * X6_3 +
           X1_2 * X2_5 * X3_1 * X4_6 * X5_3 * X6_4 -
           X1_2 * X2_5 * X3_1 * X4_6 * X5_4 * X6_3 -
           X1_2 * X2_5 * X3_3 * X4_1 * X5_4 * X6_6 +
           X1_2 * X2_5 * X3_3 * X4_1 * X5_6 * X6_4 +
           X1_2 * X2_5 * X3_3 * X4_4 * X5_1 * X6_6 -
           X1_2 * X2_5 * X3_3 * X4_4 * X5_6 * X6_1 -
           X1_2 * X2_5 * X3_3 * X4_6 * X5_1 * X6_4 +
           X1_2 * X2_5 * X3_3 * X4_6 * X5_4 * X6_1 +
           X1_2 * X2_5 * X3_4 * X4_1 * X5_3 * X6_6 -
           X1_2 * X2_5 * X3_4 * X4_1 * X5_6 * X6_3 -
           X1_2 * X2_5 * X3_4 * X4_3 * X5_1 * X6_6 +
           X1_2 * X2_5 * X3_4 * X4_3 * X5_6 * X6_1 +
           X1_2 * X2_5 * X3_4 * X4_6 * X5_1 * X6_3 -
           X1_2 * X2_5 * X3_4 * X4_6 * X5_3 * X6_1 -
           X1_2 * X2_5 * X3_6 * X4_1 * X5_3 * X6_4 +
           X1_2 * X2_5 * X3_6 * X4_1 * X5_4 * X6_3 +
           X1_2 * X2_5 * X3_6 * X4_3 * X5_1 * X6_4 -
           X1_2 * X2_5 * X3_6 * X4_3 * X5_4 * X6_1 -
           X1_2 * X2_5 * X3_6 * X4_4 * X5_1 * X6_3 +
           X1_2 * X2_5 * X3_6 * X4_4 * X5_3 * X6_1 -
           X1_2 * X2_6 * X3_1 * X4_3 * X5_4 * X6_5 +
           X1_2 * X2_6 * X3_1 * X4_3 * X5_5 * X6_4 +
           X1_2 * X2_6 * X3_1 * X4_4 * X5_3 * X6_5 -
           X1_2 * X2_6 * X3_1 * X4_4 * X5_5 * X6_3 -
           X1_2 * X2_6 * X3_1 * X4_5 * X5_3 * X6_4 +
           X1_2 * X2_6 * X3_1 * X4_5 * X5_4 * X6_3 +
           X1_2 * X2_6 * X3_3 * X4_1 * X5_4 * X6_5 -
           X1_2 * X2_6 * X3_3 * X4_1 * X5_5 * X6_4 -
           X1_2 * X2_6 * X3_3 * X4_4 * X5_1 * X6_5 +
           X1_2 * X2_6 * X3_3 * X4_4 * X5_5 * X6_1 +
           X1_2 * X2_6 * X3_3 * X4_5 * X5_1 * X6_4 -
           X1_2 * X2_6 * X3_3 * X4_5 * X5_4 * X6_1 -
           X1_2 * X2_6 * X3_4 * X4_1 * X5_3 * X6_5 +
           X1_2 * X2_6 * X3_4 * X4_1 * X5_5 * X6_3 +
           X1_2 * X2_6 * X3_4 * X4_3 * X5_1 * X6_5 -
           X1_2 * X2_6 * X3_4 * X4_3 * X5_5 * X6_1 -
           X1_2 * X2_6 * X3_4 * X4_5 * X5_1 * X6_3 +
           X1_2 * X2_6 * X3_4 * X4_5 * X5_3 * X6_1 +
           X1_2 * X2_6 * X3_5 * X4_1 * X5_3 * X6_4 -
           X1_2 * X2_6 * X3_5 * X4_1 * X5_4 * X6_3 -
           X1_2 * X2_6 * X3_5 * X4_3 * X5_1 * X6_4 +
           X1_2 * X2_6 * X3_5 * X4_3 * X5_4 * X6_1 +
           X1_2 * X2_6 * X3_5 * X4_4 * X5_1 * X6_3 -
           X1_2 * X2_6 * X3_5 * X4_4 * X5_3 * X6_1 +
           X1_3 * X2_1 * X3_2 * X4_4 * X5_5 * X6_6 -
           X1_3 * X2_1 * X3_2 * X4_4 * X5_6 * X6_5 -
           X1_3 * X2_1 * X3_2 * X4_5 * X5_4 * X6_6 +
           X1_3 * X2_1 * X3_2 * X4_5 * X5_6 * X6_4 +
           X1_3 * X2_1 * X3_2 * X4_6 * X5_4 * X6_5 -
           X1_3 * X2_1 * X3_2 * X4_6 * X5_5 * X6_4 -
           X1_3 * X2_1 * X3_4 * X4_2 * X5_5 * X6_6 +
           X1_3 * X2_1 * X3_4 * X4_2 * X5_6 * X6_5 +
           X1_3 * X2_1 * X3_4 * X4_5 * X5_2 * X6_6 -
           X1_3 * X2_1 * X3_4 * X4_5 * X5_6 * X6_2 -
           X1_3 * X2_1 * X3_4 * X4_6 * X5_2 * X6_5 +
           X1_3 * X2_1 * X3_4 * X4_6 * X5_5 * X6_2 +
           X1_3 * X2_1 * X3_5 * X4_2 * X5_4 * X6_6 -
           X1_3 * X2_1 * X3_5 * X4_2 * X5_6 * X6_4 -
           X1_3 * X2_1 * X3_5 * X4_4 * X5_2 * X6_6 +
           X1_3 * X2_1 * X3_5 * X4_4 * X5_6 * X6_2 +
           X1_3 * X2_1 * X3_5 * X4_6 * X5_2 * X6_4 -
           X1_3 * X2_1 * X3_5 * X4_6 * X5_4 * X6_2 -
           X1_3 * X2_1 * X3_6 * X4_2 * X5_4 * X6_5 +
           X1_3 * X2_1 * X3_6 * X4_2 * X5_5 * X6_4 +
           X1_3 * X2_1 * X3_6 * X4_4 * X5_2 * X6_5 -
           X1_3 * X2_1 * X3_6 * X4_4 * X5_5 * X6_2 -
           X1_3 * X2_1 * X3_6 * X4_5 * X5_2 * X6_4 +
           X1_3 * X2_1 * X3_6 * X4_5 * X5_4 * X6_2 -
           X1_3 * X2_2 * X3_1 * X4_4 * X5_5 * X6_6 +
           X1_3 * X2_2 * X3_1 * X4_4 * X5_6 * X6_5 +
           X1_3 * X2_2 * X3_1 * X4_5 * X5_4 * X6_6 -
           X1_3 * X2_2 * X3_1 * X4_5 * X5_6 * X6_4 -
           X1_3 * X2_2 * X3_1 * X4_6 * X5_4 * X6_5 +
           X1_3 * X2_2 * X3_1 * X4_6 * X5_5 * X6_4 +
           X1_3 * X2_2 * X3_4 * X4_1 * X5_5 * X6_6 -
           X1_3 * X2_2 * X3_4 * X4_1 * X5_6 * X6_5 -
           X1_3 * X2_2 * X3_4 * X4_5 * X5_1 * X6_6 +
           X1_3 * X2_2 * X3_4 * X4_5 * X5_6 * X6_1 +
           X1_3 * X2_2 * X3_4 * X4_6 * X5_1 * X6_5 -
           X1_3 * X2_2 * X3_4 * X4_6 * X5_5 * X6_1 -
           X1_3 * X2_2 * X3_5 * X4_1 * X5_4 * X6_6 +
           X1_3 * X2_2 * X3_5 * X4_1 * X5_6 * X6_4 +
           X1_3 * X2_2 * X3_5 * X4_4 * X5_1 * X6_6 -
           X1_3 * X2_2 * X3_5 * X4_4 * X5_6 * X6_1 -
           X1_3 * X2_2 * X3_5 * X4_6 * X5_1 * X6_4 +
           X1_3 * X2_2 * X3_5 * X4_6 * X5_4 * X6_1 +
           X1_3 * X2_2 * X3_6 * X4_1 * X5_4 * X6_5 -
           X1_3 * X2_2 * X3_6 * X4_1 * X5_5 * X6_4 -
           X1_3 * X2_2 * X3_6 * X4_4 * X5_1 * X6_5 +
           X1_3 * X2_2 * X3_6 * X4_4 * X5_5 * X6_1 +
           X1_3 * X2_2 * X3_6 * X4_5 * X5_1 * X6_4 -
           X1_3 * X2_2 * X3_6 * X4_5 * X5_4 * X6_1 +
           X1_3 * X2_4 * X3_1 * X4_2 * X5_5 * X6_6 -
           X1_3 * X2_4 * X3_1 * X4_2 * X5_6 * X6_5 -
           X1_3 * X2_4 * X3_1 * X4_5 * X5_2 * X6_6 +
           X1_3 * X2_4 * X3_1 * X4_5 * X5_6 * X6_2 +
           X1_3 * X2_4 * X3_1 * X4_6 * X5_2 * X6_5 -
           X1_3 * X2_4 * X3_1 * X4_6 * X5_5 * X6_2 -
           X1_3 * X2_4 * X3_2 * X4_1 * X5_5 * X6_6 +
           X1_3 * X2_4 * X3_2 * X4_1 * X5_6 * X6_5 +
           X1_3 * X2_4 * X3_2 * X4_5 * X5_1 * X6_6 -
           X1_3 * X2_4 * X3_2 * X4_5 * X5_6 * X6_1 -
           X1_3 * X2_4 * X3_2 * X4_6 * X5_1 * X6_5 +
           X1_3 * X2_4 * X3_2 * X4_6 * X5_5 * X6_1 +
           X1_3 * X2_4 * X3_5 * X4_1 * X5_2 * X6_6 -
           X1_3 * X2_4 * X3_5 * X4_1 * X5_6 * X6_2 -
           X1_3 * X2_4 * X3_5 * X4_2 * X5_1 * X6_6 +
           X1_3 * X2_4 * X3_5 * X4_2 * X5_6 * X6_1 +
           X1_3 * X2_4 * X3_5 * X4_6 * X5_1 * X6_2 -
           X1_3 * X2_4 * X3_5 * X4_6 * X5_2 * X6_1 -
           X1_3 * X2_4 * X3_6 * X4_1 * X5_2 * X6_5 +
           X1_3 * X2_4 * X3_6 * X4_1 * X5_5 * X6_2 +
           X1_3 * X2_4 * X3_6 * X4_2 * X5_1 * X6_5 -
           X1_3 * X2_4 * X3_6 * X4_2 * X5_5 * X6_1 -
           X1_3 * X2_4 * X3_6 * X4_5 * X5_1 * X6_2 +
           X1_3 * X2_4 * X3_6 * X4_5 * X5_2 * X6_1 -
           X1_3 * X2_5 * X3_1 * X4_2 * X5_4 * X6_6 +
           X1_3 * X2_5 * X3_1 * X4_2 * X5_6 * X6_4 +
           X1_3 * X2_5 * X3_1 * X4_4 * X5_2 * X6_6 -
           X1_3 * X2_5 * X3_1 * X4_4 * X5_6 * X6_2 -
           X1_3 * X2_5 * X3_1 * X4_6 * X5_2 * X6_4 +
           X1_3 * X2_5 * X3_1 * X4_6 * X5_4 * X6_2 +
           X1_3 * X2_5 * X3_2 * X4_1 * X5_4 * X6_6 -
           X1_3 * X2_5 * X3_2 * X4_1 * X5_6 * X6_4 -
           X1_3 * X2_5 * X3_2 * X4_4 * X5_1 * X6_6 +
           X1_3 * X2_5 * X3_2 * X4_4 * X5_6 * X6_1 +
           X1_3 * X2_5 * X3_2 * X4_6 * X5_1 * X6_4 -
           X1_3 * X2_5 * X3_2 * X4_6 * X5_4 * X6_1 -
           X1_3 * X2_5 * X3_4 * X4_1 * X5_2 * X6_6 +
           X1_3 * X2_5 * X3_4 * X4_1 * X5_6 * X6_2 +
           X1_3 * X2_5 * X3_4 * X4_2 * X5_1 * X6_6 -
           X1_3 * X2_5 * X3_4 * X4_2 * X5_6 * X6_1 -
           X1_3 * X2_5 * X3_4 * X4_6 * X5_1 * X6_2 +
           X1_3 * X2_5 * X3_4 * X4_6 * X5_2 * X6_1 +
           X1_3 * X2_5 * X3_6 * X4_1 * X5_2 * X6_4 -
           X1_3 * X2_5 * X3_6 * X4_1 * X5_4 * X6_2 -
           X1_3 * X2_5 * X3_6 * X4_2 * X5_1 * X6_4 +
           X1_3 * X2_5 * X3_6 * X4_2 * X5_4 * X6_1 +
           X1_3 * X2_5 * X3_6 * X4_4 * X5_1 * X6_2 -
           X1_3 * X2_5 * X3_6 * X4_4 * X5_2 * X6_1 +
           X1_3 * X2_6 * X3_1 * X4_2 * X5_4 * X6_5 -
           X1_3 * X2_6 * X3_1 * X4_2 * X5_5 * X6_4 -
           X1_3 * X2_6 * X3_1 * X4_4 * X5_2 * X6_5 +
           X1_3 * X2_6 * X3_1 * X4_4 * X5_5 * X6_2 +
           X1_3 * X2_6 * X3_1 * X4_5 * X5_2 * X6_4 -
           X1_3 * X2_6 * X3_1 * X4_5 * X5_4 * X6_2 -
           X1_3 * X2_6 * X3_2 * X4_1 * X5_4 * X6_5 +
           X1_3 * X2_6 * X3_2 * X4_1 * X5_5 * X6_4 +
           X1_3 * X2_6 * X3_2 * X4_4 * X5_1 * X6_5 -
           X1_3 * X2_6 * X3_2 * X4_4 * X5_5 * X6_1 -
           X1_3 * X2_6 * X3_2 * X4_5 * X5_1 * X6_4 +
           X1_3 * X2_6 * X3_2 * X4_5 * X5_4 * X6_1 +
           X1_3 * X2_6 * X3_4 * X4_1 * X5_2 * X6_5 -
           X1_3 * X2_6 * X3_4 * X4_1 * X5_5 * X6_2 -
           X1_3 * X2_6 * X3_4 * X4_2 * X5_1 * X6_5 +
           X1_3 * X2_6 * X3_4 * X4_2 * X5_5 * X6_1 +
           X1_3 * X2_6 * X3_4 * X4_5 * X5_1 * X6_2 -
           X1_3 * X2_6 * X3_4 * X4_5 * X5_2 * X6_1 -
           X1_3 * X2_6 * X3_5 * X4_1 * X5_2 * X6_4 +
           X1_3 * X2_6 * X3_5 * X4_1 * X5_4 * X6_2 +
           X1_3 * X2_6 * X3_5 * X4_2 * X5_1 * X6_4 -
           X1_3 * X2_6 * X3_5 * X4_2 * X5_4 * X6_1 -
           X1_3 * X2_6 * X3_5 * X4_4 * X5_1 * X6_2 +
           X1_3 * X2_6 * X3_5 * X4_4 * X5_2 * X6_1 -
           X1_4 * X2_1 * X3_2 * X4_3 * X5_5 * X6_6 +
           X1_4 * X2_1 * X3_2 * X4_3 * X5_6 * X6_5 +
           X1_4 * X2_1 * X3_2 * X4_5 * X5_3 * X6_6 -
           X1_4 * X2_1 * X3_2 * X4_5 * X5_6 * X6_3 -
           X1_4 * X2_1 * X3_2 * X4_6 * X5_3 * X6_5 +
           X1_4 * X2_1 * X3_2 * X4_6 * X5_5 * X6_3 +
           X1_4 * X2_1 * X3_3 * X4_2 * X5_5 * X6_6 -
           X1_4 * X2_1 * X3_3 * X4_2 * X5_6 * X6_5 -
           X1_4 * X2_1 * X3_3 * X4_5 * X5_2 * X6_6 +
           X1_4 * X2_1 * X3_3 * X4_5 * X5_6 * X6_2 +
           X1_4 * X2_1 * X3_3 * X4_6 * X5_2 * X6_5 -
           X1_4 * X2_1 * X3_3 * X4_6 * X5_5 * X6_2 -
           X1_4 * X2_1 * X3_5 * X4_2 * X5_3 * X6_6 +
           X1_4 * X2_1 * X3_5 * X4_2 * X5_6 * X6_3 +
           X1_4 * X2_1 * X3_5 * X4_3 * X5_2 * X6_6 -
           X1_4 * X2_1 * X3_5 * X4_3 * X5_6 * X6_2 -
           X1_4 * X2_1 * X3_5 * X4_6 * X5_2 * X6_3 +
           X1_4 * X2_1 * X3_5 * X4_6 * X5_3 * X6_2 +
           X1_4 * X2_1 * X3_6 * X4_2 * X5_3 * X6_5 -
           X1_4 * X2_1 * X3_6 * X4_2 * X5_5 * X6_3 -
           X1_4 * X2_1 * X3_6 * X4_3 * X5_2 * X6_5 +
           X1_4 * X2_1 * X3_6 * X4_3 * X5_5 * X6_2 +
           X1_4 * X2_1 * X3_6 * X4_5 * X5_2 * X6_3 -
           X1_4 * X2_1 * X3_6 * X4_5 * X5_3 * X6_2 +
           X1_4 * X2_2 * X3_1 * X4_3 * X5_5 * X6_6 -
           X1_4 * X2_2 * X3_1 * X4_3 * X5_6 * X6_5 -
           X1_4 * X2_2 * X3_1 * X4_5 * X5_3 * X6_6 +
           X1_4 * X2_2 * X3_1 * X4_5 * X5_6 * X6_3 +
           X1_4 * X2_2 * X3_1 * X4_6 * X5_3 * X6_5 -
           X1_4 * X2_2 * X3_1 * X4_6 * X5_5 * X6_3 -
           X1_4 * X2_2 * X3_3 * X4_1 * X5_5 * X6_6 +
           X1_4 * X2_2 * X3_3 * X4_1 * X5_6 * X6_5 +
           X1_4 * X2_2 * X3_3 * X4_5 * X5_1 * X6_6 -
           X1_4 * X2_2 * X3_3 * X4_5 * X5_6 * X6_1 -
           X1_4 * X2_2 * X3_3 * X4_6 * X5_1 * X6_5 +
           X1_4 * X2_2 * X3_3 * X4_6 * X5_5 * X6_1 +
           X1_4 * X2_2 * X3_5 * X4_1 * X5_3 * X6_6 -
           X1_4 * X2_2 * X3_5 * X4_1 * X5_6 * X6_3 -
           X1_4 * X2_2 * X3_5 * X4_3 * X5_1 * X6_6 +
           X1_4 * X2_2 * X3_5 * X4_3 * X5_6 * X6_1 +
           X1_4 * X2_2 * X3_5 * X4_6 * X5_1 * X6_3 -
           X1_4 * X2_2 * X3_5 * X4_6 * X5_3 * X6_1 -
           X1_4 * X2_2 * X3_6 * X4_1 * X5_3 * X6_5 +
           X1_4 * X2_2 * X3_6 * X4_1 * X5_5 * X6_3 +
           X1_4 * X2_2 * X3_6 * X4_3 * X5_1 * X6_5 -
           X1_4 * X2_2 * X3_6 * X4_3 * X5_5 * X6_1 -
           X1_4 * X2_2 * X3_6 * X4_5 * X5_1 * X6_3 +
           X1_4 * X2_2 * X3_6 * X4_5 * X5_3 * X6_1 -
           X1_4 * X2_3 * X3_1 * X4_2 * X5_5 * X6_6 +
           X1_4 * X2_3 * X3_1 * X4_2 * X5_6 * X6_5 +
           X1_4 * X2_3 * X3_1 * X4_5 * X5_2 * X6_6 -
           X1_4 * X2_3 * X3_1 * X4_5 * X5_6 * X6_2 -
           X1_4 * X2_3 * X3_1 * X4_6 * X5_2 * X6_5 +
           X1_4 * X2_3 * X3_1 * X4_6 * X5_5 * X6_2 +
           X1_4 * X2_3 * X3_2 * X4_1 * X5_5 * X6_6 -
           X1_4 * X2_3 * X3_2 * X4_1 * X5_6 * X6_5 -
           X1_4 * X2_3 * X3_2 * X4_5 * X5_1 * X6_6 +
           X1_4 * X2_3 * X3_2 * X4_5 * X5_6 * X6_1 +
           X1_4 * X2_3 * X3_2 * X4_6 * X5_1 * X6_5 -
           X1_4 * X2_3 * X3_2 * X4_6 * X5_5 * X6_1 -
           X1_4 * X2_3 * X3_5 * X4_1 * X5_2 * X6_6 +
           X1_4 * X2_3 * X3_5 * X4_1 * X5_6 * X6_2 +
           X1_4 * X2_3 * X3_5 * X4_2 * X5_1 * X6_6 -
           X1_4 * X2_3 * X3_5 * X4_2 * X5_6 * X6_1 -
           X1_4 * X2_3 * X3_5 * X4_6 * X5_1 * X6_2 +
           X1_4 * X2_3 * X3_5 * X4_6 * X5_2 * X6_1 +
           X1_4 * X2_3 * X3_6 * X4_1 * X5_2 * X6_5 -
           X1_4 * X2_3 * X3_6 * X4_1 * X5_5 * X6_2 -
           X1_4 * X2_3 * X3_6 * X4_2 * X5_1 * X6_5 +
           X1_4 * X2_3 * X3_6 * X4_2 * X5_5 * X6_1 +
           X1_4 * X2_3 * X3_6 * X4_5 * X5_1 * X6_2 -
           X1_4 * X2_3 * X3_6 * X4_5 * X5_2 * X6_1 +
           X1_4 * X2_5 * X3_1 * X4_2 * X5_3 * X6_6 -
           X1_4 * X2_5 * X3_1 * X4_2 * X5_6 * X6_3 -
           X1_4 * X2_5 * X3_1 * X4_3 * X5_2 * X6_6 +
           X1_4 * X2_5 * X3_1 * X4_3 * X5_6 * X6_2 +
           X1_4 * X2_5 * X3_1 * X4_6 * X5_2 * X6_3 -
           X1_4 * X2_5 * X3_1 * X4_6 * X5_3 * X6_2 -
           X1_4 * X2_5 * X3_2 * X4_1 * X5_3 * X6_6 +
           X1_4 * X2_5 * X3_2 * X4_1 * X5_6 * X6_3 +
           X1_4 * X2_5 * X3_2 * X4_3 * X5_1 * X6_6 -
           X1_4 * X2_5 * X3_2 * X4_3 * X5_6 * X6_1 -
           X1_4 * X2_5 * X3_2 * X4_6 * X5_1 * X6_3 +
           X1_4 * X2_5 * X3_2 * X4_6 * X5_3 * X6_1 +
           X1_4 * X2_5 * X3_3 * X4_1 * X5_2 * X6_6 -
           X1_4 * X2_5 * X3_3 * X4_1 * X5_6 * X6_2 -
           X1_4 * X2_5 * X3_3 * X4_2 * X5_1 * X6_6 +
           X1_4 * X2_5 * X3_3 * X4_2 * X5_6 * X6_1 +
           X1_4 * X2_5 * X3_3 * X4_6 * X5_1 * X6_2 -
           X1_4 * X2_5 * X3_3 * X4_6 * X5_2 * X6_1 -
           X1_4 * X2_5 * X3_6 * X4_1 * X5_2 * X6_3 +
           X1_4 * X2_5 * X3_6 * X4_1 * X5_3 * X6_2 +
           X1_4 * X2_5 * X3_6 * X4_2 * X5_1 * X6_3 -
           X1_4 * X2_5 * X3_6 * X4_2 * X5_3 * X6_1 -
           X1_4 * X2_5 * X3_6 * X4_3 * X5_1 * X6_2 +
           X1_4 * X2_5 * X3_6 * X4_3 * X5_2 * X6_1 -
           X1_4 * X2_6 * X3_1 * X4_2 * X5_3 * X6_5 +
           X1_4 * X2_6 * X3_1 * X4_2 * X5_5 * X6_3 +
           X1_4 * X2_6 * X3_1 * X4_3 * X5_2 * X6_5 -
           X1_4 * X2_6 * X3_1 * X4_3 * X5_5 * X6_2 -
           X1_4 * X2_6 * X3_1 * X4_5 * X5_2 * X6_3 +
           X1_4 * X2_6 * X3_1 * X4_5 * X5_3 * X6_2 +
           X1_4 * X2_6 * X3_2 * X4_1 * X5_3 * X6_5 -
           X1_4 * X2_6 * X3_2 * X4_1 * X5_5 * X6_3 -
           X1_4 * X2_6 * X3_2 * X4_3 * X5_1 * X6_5 +
           X1_4 * X2_6 * X3_2 * X4_3 * X5_5 * X6_1 +
           X1_4 * X2_6 * X3_2 * X4_5 * X5_1 * X6_3 -
           X1_4 * X2_6 * X3_2 * X4_5 * X5_3 * X6_1 -
           X1_4 * X2_6 * X3_3 * X4_1 * X5_2 * X6_5 +
           X1_4 * X2_6 * X3_3 * X4_1 * X5_5 * X6_2 +
           X1_4 * X2_6 * X3_3 * X4_2 * X5_1 * X6_5 -
           X1_4 * X2_6 * X3_3 * X4_2 * X5_5 * X6_1 -
           X1_4 * X2_6 * X3_3 * X4_5 * X5_1 * X6_2 +
           X1_4 * X2_6 * X3_3 * X4_5 * X5_2 * X6_1 +
           X1_4 * X2_6 * X3_5 * X4_1 * X5_2 * X6_3 -
           X1_4 * X2_6 * X3_5 * X4_1 * X5_3 * X6_2 -
           X1_4 * X2_6 * X3_5 * X4_2 * X5_1 * X6_3 +
           X1_4 * X2_6 * X3_5 * X4_2 * X5_3 * X6_1 +
           X1_4 * X2_6 * X3_5 * X4_3 * X5_1 * X6_2 -
           X1_4 * X2_6 * X3_5 * X4_3 * X5_2 * X6_1 +
           X1_5 * X2_1 * X3_2 * X4_3 * X5_4 * X6_6 -
           X1_5 * X2_1 * X3_2 * X4_3 * X5_6 * X6_4 -
           X1_5 * X2_1 * X3_2 * X4_4 * X5_3 * X6_6 +
           X1_5 * X2_1 * X3_2 * X4_4 * X5_6 * X6_3 +
           X1_5 * X2_1 * X3_2 * X4_6 * X5_3 * X6_4 -
           X1_5 * X2_1 * X3_2 * X4_6 * X5_4 * X6_3 -
           X1_5 * X2_1 * X3_3 * X4_2 * X5_4 * X6_6 +
           X1_5 * X2_1 * X3_3 * X4_2 * X5_6 * X6_4 +
           X1_5 * X2_1 * X3_3 * X4_4 * X5_2 * X6_6 -
           X1_5 * X2_1 * X3_3 * X4_4 * X5_6 * X6_2 -
           X1_5 * X2_1 * X3_3 * X4_6 * X5_2 * X6_4 +
           X1_5 * X2_1 * X3_3 * X4_6 * X5_4 * X6_2 +
           X1_5 * X2_1 * X3_4 * X4_2 * X5_3 * X6_6 -
           X1_5 * X2_1 * X3_4 * X4_2 * X5_6 * X6_3 -
           X1_5 * X2_1 * X3_4 * X4_3 * X5_2 * X6_6 +
           X1_5 * X2_1 * X3_4 * X4_3 * X5_6 * X6_2 +
           X1_5 * X2_1 * X3_4 * X4_6 * X5_2 * X6_3 -
           X1_5 * X2_1 * X3_4 * X4_6 * X5_3 * X6_2 -
           X1_5 * X2_1 * X3_6 * X4_2 * X5_3 * X6_4 +
           X1_5 * X2_1 * X3_6 * X4_2 * X5_4 * X6_3 +
           X1_5 * X2_1 * X3_6 * X4_3 * X5_2 * X6_4 -
           X1_5 * X2_1 * X3_6 * X4_3 * X5_4 * X6_2 -
           X1_5 * X2_1 * X3_6 * X4_4 * X5_2 * X6_3 +
           X1_5 * X2_1 * X3_6 * X4_4 * X5_3 * X6_2 -
           X1_5 * X2_2 * X3_1 * X4_3 * X5_4 * X6_6 +
           X1_5 * X2_2 * X3_1 * X4_3 * X5_6 * X6_4 +
           X1_5 * X2_2 * X3_1 * X4_4 * X5_3 * X6_6 -
           X1_5 * X2_2 * X3_1 * X4_4 * X5_6 * X6_3 -
           X1_5 * X2_2 * X3_1 * X4_6 * X5_3 * X6_4 +
           X1_5 * X2_2 * X3_1 * X4_6 * X5_4 * X6_3 +
           X1_5 * X2_2 * X3_3 * X4_1 * X5_4 * X6_6 -
           X1_5 * X2_2 * X3_3 * X4_1 * X5_6 * X6_4 -
           X1_5 * X2_2 * X3_3 * X4_4 * X5_1 * X6_6 +
           X1_5 * X2_2 * X3_3 * X4_4 * X5_6 * X6_1 +
           X1_5 * X2_2 * X3_3 * X4_6 * X5_1 * X6_4 -
           X1_5 * X2_2 * X3_3 * X4_6 * X5_4 * X6_1 -
           X1_5 * X2_2 * X3_4 * X4_1 * X5_3 * X6_6 +
           X1_5 * X2_2 * X3_4 * X4_1 * X5_6 * X6_3 +
           X1_5 * X2_2 * X3_4 * X4_3 * X5_1 * X6_6 -
           X1_5 * X2_2 * X3_4 * X4_3 * X5_6 * X6_1 -
           X1_5 * X2_2 * X3_4 * X4_6 * X5_1 * X6_3 +
           X1_5 * X2_2 * X3_4 * X4_6 * X5_3 * X6_1 +
           X1_5 * X2_2 * X3_6 * X4_1 * X5_3 * X6_4 -
           X1_5 * X2_2 * X3_6 * X4_1 * X5_4 * X6_3 -
           X1_5 * X2_2 * X3_6 * X4_3 * X5_1 * X6_4 +
           X1_5 * X2_2 * X3_6 * X4_3 * X5_4 * X6_1 +
           X1_5 * X2_2 * X3_6 * X4_4 * X5_1 * X6_3 -
           X1_5 * X2_2 * X3_6 * X4_4 * X5_3 * X6_1 +
           X1_5 * X2_3 * X3_1 * X4_2 * X5_4 * X6_6 -
           X1_5 * X2_3 * X3_1 * X4_2 * X5_6 * X6_4 -
           X1_5 * X2_3 * X3_1 * X4_4 * X5_2 * X6_6 +
           X1_5 * X2_3 * X3_1 * X4_4 * X5_6 * X6_2 +
           X1_5 * X2_3 * X3_1 * X4_6 * X5_2 * X6_4 -
           X1_5 * X2_3 * X3_1 * X4_6 * X5_4 * X6_2 -
           X1_5 * X2_3 * X3_2 * X4_1 * X5_4 * X6_6 +
           X1_5 * X2_3 * X3_2 * X4_1 * X5_6 * X6_4 +
           X1_5 * X2_3 * X3_2 * X4_4 * X5_1 * X6_6 -
           X1_5 * X2_3 * X3_2 * X4_4 * X5_6 * X6_1 -
           X1_5 * X2_3 * X3_2 * X4_6 * X5_1 * X6_4 +
           X1_5 * X2_3 * X3_2 * X4_6 * X5_4 * X6_1 +
           X1_5 * X2_3 * X3_4 * X4_1 * X5_2 * X6_6 -
           X1_5 * X2_3 * X3_4 * X4_1 * X5_6 * X6_2 -
           X1_5 * X2_3 * X3_4 * X4_2 * X5_1 * X6_6 +
           X1_5 * X2_3 * X3_4 * X4_2 * X5_6 * X6_1 +
           X1_5 * X2_3 * X3_4 * X4_6 * X5_1 * X6_2 -
           X1_5 * X2_3 * X3_4 * X4_6 * X5_2 * X6_1 -
           X1_5 * X2_3 * X3_6 * X4_1 * X5_2 * X6_4 +
           X1_5 * X2_3 * X3_6 * X4_1 * X5_4 * X6_2 +
           X1_5 * X2_3 * X3_6 * X4_2 * X5_1 * X6_4 -
           X1_5 * X2_3 * X3_6 * X4_2 * X5_4 * X6_1 -
           X1_5 * X2_3 * X3_6 * X4_4 * X5_1 * X6_2 +
           X1_5 * X2_3 * X3_6 * X4_4 * X5_2 * X6_1 -
           X1_5 * X2_4 * X3_1 * X4_2 * X5_3 * X6_6 +
           X1_5 * X2_4 * X3_1 * X4_2 * X5_6 * X6_3 +
           X1_5 * X2_4 * X3_1 * X4_3 * X5_2 * X6_6 -
           X1_5 * X2_4 * X3_1 * X4_3 * X5_6 * X6_2 -
           X1_5 * X2_4 * X3_1 * X4_6 * X5_2 * X6_3 +
           X1_5 * X2_4 * X3_1 * X4_6 * X5_3 * X6_2 +
           X1_5 * X2_4 * X3_2 * X4_1 * X5_3 * X6_6 -
           X1_5 * X2_4 * X3_2 * X4_1 * X5_6 * X6_3 -
           X1_5 * X2_4 * X3_2 * X4_3 * X5_1 * X6_6 +
           X1_5 * X2_4 * X3_2 * X4_3 * X5_6 * X6_1 +
           X1_5 * X2_4 * X3_2 * X4_6 * X5_1 * X6_3 -
           X1_5 * X2_4 * X3_2 * X4_6 * X5_3 * X6_1 -
           X1_5 * X2_4 * X3_3 * X4_1 * X5_2 * X6_6 +
           X1_5 * X2_4 * X3_3 * X4_1 * X5_6 * X6_2 +
           X1_5 * X2_4 * X3_3 * X4_2 * X5_1 * X6_6 -
           X1_5 * X2_4 * X3_3 * X4_2 * X5_6 * X6_1 -
           X1_5 * X2_4 * X3_3 * X4_6 * X5_1 * X6_2 +
           X1_5 * X2_4 * X3_3 * X4_6 * X5_2 * X6_1 +
           X1_5 * X2_4 * X3_6 * X4_1 * X5_2 * X6_3 -
           X1_5 * X2_4 * X3_6 * X4_1 * X5_3 * X6_2 -
           X1_5 * X2_4 * X3_6 * X4_2 * X5_1 * X6_3 +
           X1_5 * X2_4 * X3_6 * X4_2 * X5_3 * X6_1 +
           X1_5 * X2_4 * X3_6 * X4_3 * X5_1 * X6_2 -
           X1_5 * X2_4 * X3_6 * X4_3 * X5_2 * X6_1 +
           X1_5 * X2_6 * X3_1 * X4_2 * X5_3 * X6_4 -
           X1_5 * X2_6 * X3_1 * X4_2 * X5_4 * X6_3 -
           X1_5 * X2_6 * X3_1 * X4_3 * X5_2 * X6_4 +
           X1_5 * X2_6 * X3_1 * X4_3 * X5_4 * X6_2 +
           X1_5 * X2_6 * X3_1 * X4_4 * X5_2 * X6_3 -
           X1_5 * X2_6 * X3_1 * X4_4 * X5_3 * X6_2 -
           X1_5 * X2_6 * X3_2 * X4_1 * X5_3 * X6_4 +
           X1_5 * X2_6 * X3_2 * X4_1 * X5_4 * X6_3 +
           X1_5 * X2_6 * X3_2 * X4_3 * X5_1 * X6_4 -
           X1_5 * X2_6 * X3_2 * X4_3 * X5_4 * X6_1 -
           X1_5 * X2_6 * X3_2 * X4_4 * X5_1 * X6_3 +
           X1_5 * X2_6 * X3_2 * X4_4 * X5_3 * X6_1 +
           X1_5 * X2_6 * X3_3 * X4_1 * X5_2 * X6_4 -
           X1_5 * X2_6 * X3_3 * X4_1 * X5_4 * X6_2 -
           X1_5 * X2_6 * X3_3 * X4_2 * X5_1 * X6_4 +
           X1_5 * X2_6 * X3_3 * X4_2 * X5_4 * X6_1 +
           X1_5 * X2_6 * X3_3 * X4_4 * X5_1 * X6_2 -
           X1_5 * X2_6 * X3_3 * X4_4 * X5_2 * X6_1 -
           X1_5 * X2_6 * X3_4 * X4_1 * X5_2 * X6_3 +
           X1_5 * X2_6 * X3_4 * X4_1 * X5_3 * X6_2 +
           X1_5 * X2_6 * X3_4 * X4_2 * X5_1 * X6_3 -
           X1_5 * X2_6 * X3_4 * X4_2 * X5_3 * X6_1 -
           X1_5 * X2_6 * X3_4 * X4_3 * X5_1 * X6_2 +
           X1_5 * X2_6 * X3_4 * X4_3 * X5_2 * X6_1 -
           X1_6 * X2_1 * X3_2 * X4_3 * X5_4 * X6_5 +
           X1_6 * X2_1 * X3_2 * X4_3 * X5_5 * X6_4 +
           X1_6 * X2_1 * X3_2 * X4_4 * X5_3 * X6_5 -
           X1_6 * X2_1 * X3_2 * X4_4 * X5_5 * X6_3 -
           X1_6 * X2_1 * X3_2 * X4_5 * X5_3 * X6_4 +
           X1_6 * X2_1 * X3_2 * X4_5 * X5_4 * X6_3 +
           X1_6 * X2_1 * X3_3 * X4_2 * X5_4 * X6_5 -
           X1_6 * X2_1 * X3_3 * X4_2 * X5_5 * X6_4 -
           X1_6 * X2_1 * X3_3 * X4_4 * X5_2 * X6_5 +
           X1_6 * X2_1 * X3_3 * X4_4 * X5_5 * X6_2 +
           X1_6 * X2_1 * X3_3 * X4_5 * X5_2 * X6_4 -
           X1_6 * X2_1 * X3_3 * X4_5 * X5_4 * X6_2 -
           X1_6 * X2_1 * X3_4 * X4_2 * X5_3 * X6_5 +
           X1_6 * X2_1 * X3_4 * X4_2 * X5_5 * X6_3 +
           X1_6 * X2_1 * X3_4 * X4_3 * X5_2 * X6_5 -
           X1_6 * X2_1 * X3_4 * X4_3 * X5_5 * X6_2 -
           X1_6 * X2_1 * X3_4 * X4_5 * X5_2 * X6_3 +
           X1_6 * X2_1 * X3_4 * X4_5 * X5_3 * X6_2 +
           X1_6 * X2_1 * X3_5 * X4_2 * X5_3 * X6_4 -
           X1_6 * X2_1 * X3_5 * X4_2 * X5_4 * X6_3 -
           X1_6 * X2_1 * X3_5 * X4_3 * X5_2 * X6_4 +
           X1_6 * X2_1 * X3_5 * X4_3 * X5_4 * X6_2 +
           X1_6 * X2_1 * X3_5 * X4_4 * X5_2 * X6_3 -
           X1_6 * X2_1 * X3_5 * X4_4 * X5_3 * X6_2 +
           X1_6 * X2_2 * X3_1 * X4_3 * X5_4 * X6_5 -
           X1_6 * X2_2 * X3_1 * X4_3 * X5_5 * X6_4 -
           X1_6 * X2_2 * X3_1 * X4_4 * X5_3 * X6_5 +
           X1_6 * X2_2 * X3_1 * X4_4 * X5_5 * X6_3 +
           X1_6 * X2_2 * X3_1 * X4_5 * X5_3 * X6_4 -
           X1_6 * X2_2 * X3_1 * X4_5 * X5_4 * X6_3 -
           X1_6 * X2_2 * X3_3 * X4_1 * X5_4 * X6_5 +
           X1_6 * X2_2 * X3_3 * X4_1 * X5_5 * X6_4 +
           X1_6 * X2_2 * X3_3 * X4_4 * X5_1 * X6_5 -
           X1_6 * X2_2 * X3_3 * X4_4 * X5_5 * X6_1 -
           X1_6 * X2_2 * X3_3 * X4_5 * X5_1 * X6_4 +
           X1_6 * X2_2 * X3_3 * X4_5 * X5_4 * X6_1 +
           X1_6 * X2_2 * X3_4 * X4_1 * X5_3 * X6_5 -
           X1_6 * X2_2 * X3_4 * X4_1 * X5_5 * X6_3 -
           X1_6 * X2_2 * X3_4 * X4_3 * X5_1 * X6_5 +
           X1_6 * X2_2 * X3_4 * X4_3 * X5_5 * X6_1 +
           X1_6 * X2_2 * X3_4 * X4_5 * X5_1 * X6_3 -
           X1_6 * X2_2 * X3_4 * X4_5 * X5_3 * X6_1 -
           X1_6 * X2_2 * X3_5 * X4_1 * X5_3 * X6_4 +
           X1_6 * X2_2 * X3_5 * X4_1 * X5_4 * X6_3 +
           X1_6 * X2_2 * X3_5 * X4_3 * X5_1 * X6_4 -
           X1_6 * X2_2 * X3_5 * X4_3 * X5_4 * X6_1 -
           X1_6 * X2_2 * X3_5 * X4_4 * X5_1 * X6_3 +
           X1_6 * X2_2 * X3_5 * X4_4 * X5_3 * X6_1 -
           X1_6 * X2_3 * X3_1 * X4_2 * X5_4 * X6_5 +
           X1_6 * X2_3 * X3_1 * X4_2 * X5_5 * X6_4 +
           X1_6 * X2_3 * X3_1 * X4_4 * X5_2 * X6_5 -
           X1_6 * X2_3 * X3_1 * X4_4 * X5_5 * X6_2 -
           X1_6 * X2_3 * X3_1 * X4_5 * X5_2 * X6_4 +
           X1_6 * X2_3 * X3_1 * X4_5 * X5_4 * X6_2 +
           X1_6 * X2_3 * X3_2 * X4_1 * X5_4 * X6_5 -
           X1_6 * X2_3 * X3_2 * X4_1 * X5_5 * X6_4 -
           X1_6 * X2_3 * X3_2 * X4_4 * X5_1 * X6_5 +
           X1_6 * X2_3 * X3_2 * X4_4 * X5_5 * X6_1 +
           X1_6 * X2_3 * X3_2 * X4_5 * X5_1 * X6_4 -
           X1_6 * X2_3 * X3_2 * X4_5 * X5_4 * X6_1 -
           X1_6 * X2_3 * X3_4 * X4_1 * X5_2 * X6_5 +
           X1_6 * X2_3 * X3_4 * X4_1 * X5_5 * X6_2 +
           X1_6 * X2_3 * X3_4 * X4_2 * X5_1 * X6_5 -
           X1_6 * X2_3 * X3_4 * X4_2 * X5_5 * X6_1 -
           X1_6 * X2_3 * X3_4 * X4_5 * X5_1 * X6_2 +
           X1_6 * X2_3 * X3_4 * X4_5 * X5_2 * X6_1 +
           X1_6 * X2_3 * X3_5 * X4_1 * X5_2 * X6_4 -
           X1_6 * X2_3 * X3_5 * X4_1 * X5_4 * X6_2 -
           X1_6 * X2_3 * X3_5 * X4_2 * X5_1 * X6_4 +
           X1_6 * X2_3 * X3_5 * X4_2 * X5_4 * X6_1 +
           X1_6 * X2_3 * X3_5 * X4_4 * X5_1 * X6_2 -
           X1_6 * X2_3 * X3_5 * X4_4 * X5_2 * X6_1 +
           X1_6 * X2_4 * X3_1 * X4_2 * X5_3 * X6_5 -
           X1_6 * X2_4 * X3_1 * X4_2 * X5_5 * X6_3 -
           X1_6 * X2_4 * X3_1 * X4_3 * X5_2 * X6_5 +
           X1_6 * X2_4 * X3_1 * X4_3 * X5_5 * X6_2 +
           X1_6 * X2_4 * X3_1 * X4_5 * X5_2 * X6_3 -
           X1_6 * X2_4 * X3_1 * X4_5 * X5_3 * X6_2 -
           X1_6 * X2_4 * X3_2 * X4_1 * X5_3 * X6_5 +
           X1_6 * X2_4 * X3_2 * X4_1 * X5_5 * X6_3 +
           X1_6 * X2_4 * X3_2 * X4_3 * X5_1 * X6_5 -
           X1_6 * X2_4 * X3_2 * X4_3 * X5_5 * X6_1 -
           X1_6 * X2_4 * X3_2 * X4_5 * X5_1 * X6_3 +
           X1_6 * X2_4 * X3_2 * X4_5 * X5_3 * X6_1 +
           X1_6 * X2_4 * X3_3 * X4_1 * X5_2 * X6_5 -
           X1_6 * X2_4 * X3_3 * X4_1 * X5_5 * X6_2 -
           X1_6 * X2_4 * X3_3 * X4_2 * X5_1 * X6_5 +
           X1_6 * X2_4 * X3_3 * X4_2 * X5_5 * X6_1 +
           X1_6 * X2_4 * X3_3 * X4_5 * X5_1 * X6_2 -
           X1_6 * X2_4 * X3_3 * X4_5 * X5_2 * X6_1 -
           X1_6 * X2_4 * X3_5 * X4_1 * X5_2 * X6_3 +
           X1_6 * X2_4 * X3_5 * X4_1 * X5_3 * X6_2 +
           X1_6 * X2_4 * X3_5 * X4_2 * X5_1 * X6_3 -
           X1_6 * X2_4 * X3_5 * X4_2 * X5_3 * X6_1 -
           X1_6 * X2_4 * X3_5 * X4_3 * X5_1 * X6_2 +
           X1_6 * X2_4 * X3_5 * X4_3 * X5_2 * X6_1 -
           X1_6 * X2_5 * X3_1 * X4_2 * X5_3 * X6_4 +
           X1_6 * X2_5 * X3_1 * X4_2 * X5_4 * X6_3 +
           X1_6 * X2_5 * X3_1 * X4_3 * X5_2 * X6_4 -
           X1_6 * X2_5 * X3_1 * X4_3 * X5_4 * X6_2 -
           X1_6 * X2_5 * X3_1 * X4_4 * X5_2 * X6_3 +
           X1_6 * X2_5 * X3_1 * X4_4 * X5_3 * X6_2 +
           X1_6 * X2_5 * X3_2 * X4_1 * X5_3 * X6_4 -
           X1_6 * X2_5 * X3_2 * X4_1 * X5_4 * X6_3 -
           X1_6 * X2_5 * X3_2 * X4_3 * X5_1 * X6_4 +
           X1_6 * X2_5 * X3_2 * X4_3 * X5_4 * X6_1 +
           X1_6 * X2_5 * X3_2 * X4_4 * X5_1 * X6_3 -
           X1_6 * X2_5 * X3_2 * X4_4 * X5_3 * X6_1 -
           X1_6 * X2_5 * X3_3 * X4_1 * X5_2 * X6_4 +
           X1_6 * X2_5 * X3_3 * X4_1 * X5_4 * X6_2 +
           X1_6 * X2_5 * X3_3 * X4_2 * X5_1 * X6_4 -
           X1_6 * X2_5 * X3_3 * X4_2 * X5_4 * X6_1 -
           X1_6 * X2_5 * X3_3 * X4_4 * X5_1 * X6_2 +
           X1_6 * X2_5 * X3_3 * X4_4 * X5_2 * X6_1 +
           X1_6 * X2_5 * X3_4 * X4_1 * X5_2 * X6_3 -
           X1_6 * X2_5 * X3_4 * X4_1 * X5_3 * X6_2 -
           X1_6 * X2_5 * X3_4 * X4_2 * X5_1 * X6_3 +
           X1_6 * X2_5 * X3_4 * X4_2 * X5_3 * X6_1 +
           X1_6 * X2_5 * X3_4 * X4_3 * X5_1 * X6_2 -
           X1_6 * X2_5 * X3_4 * X4_3 * X5_2 * X6_1;
  }
  std::cout << "unsupported size " << X.m() << std::endl;
  NOT_POSSIBLE;
  return 0;
}

template <typename type>
type det(const symd<type>& X) {
  ASSERT(X.m() == X.n());
  if (X.m() == 1) return X(0, 0);
  if (X.m() == 2) return X(1, 1) * X(0, 0) - X(0, 1) * X(1, 0);
  if (X.m() == 3)
    return X(0, 0) * (X(2, 2) * X(1, 1) - X(2, 1) * X(1, 2)) -
           X(1, 0) * (X(2, 2) * X(0, 1) - X(2, 1) * X(0, 2)) +
           X(2, 0) * (X(1, 2) * X(0, 1) - X(1, 1) * X(0, 2));
  const type X1_1 = X(0, 0), X1_2 = X(0, 1), X1_3 = X(0, 2), X1_4 = X(0, 3);
  const type X2_1 = X(1, 0), X2_2 = X(1, 1), X2_3 = X(1, 2), X2_4 = X(1, 3);
  const type X3_1 = X(2, 0), X3_2 = X(2, 1), X3_3 = X(2, 2), X3_4 = X(2, 3);
  const type X4_1 = X(3, 0), X4_2 = X(3, 1), X4_3 = X(3, 2), X4_4 = X(3, 3);
  if (X.m() == 4) {
    return X1_1 * X2_2 * X3_3 * X4_4 - X1_1 * X2_2 * X3_4 * X4_3 -
           X1_1 * X2_3 * X3_2 * X4_4 + X1_1 * X2_3 * X3_4 * X4_2 +
           X1_1 * X2_4 * X3_2 * X4_3 - X1_1 * X2_4 * X3_3 * X4_2 -
           X1_2 * X2_1 * X3_3 * X4_4 + X1_2 * X2_1 * X3_4 * X4_3 +
           X1_2 * X2_3 * X3_1 * X4_4 - X1_2 * X2_3 * X3_4 * X4_1 -
           X1_2 * X2_4 * X3_1 * X4_3 + X1_2 * X2_4 * X3_3 * X4_1 +
           X1_3 * X2_1 * X3_2 * X4_4 - X1_3 * X2_1 * X3_4 * X4_2 -
           X1_3 * X2_2 * X3_1 * X4_4 + X1_3 * X2_2 * X3_4 * X4_1 +
           X1_3 * X2_4 * X3_1 * X4_2 - X1_3 * X2_4 * X3_2 * X4_1 -
           X1_4 * X2_1 * X3_2 * X4_3 + X1_4 * X2_1 * X3_3 * X4_2 +
           X1_4 * X2_2 * X3_1 * X4_3 - X1_4 * X2_2 * X3_3 * X4_1 -
           X1_4 * X2_3 * X3_1 * X4_2 + X1_4 * X2_3 * X3_2 * X4_1;
  }
  std::cout << "unsupported size " << X.m() << std::endl;
  NOT_POSSIBLE;
  return 0;
}

template <int N, typename T>
T det(const mats<N, N, T>& X) {
  if (N == 1) return X(0, 0);
  if (N == 2) return X(1, 1) * X(0, 0) - X(0, 1) * X(1, 0);
  if (N == 3)
    return X(0, 0) * (X(2, 2) * X(1, 1) - X(2, 1) * X(1, 2)) -
           X(1, 0) * (X(2, 2) * X(0, 1) - X(2, 1) * X(0, 2)) +
           X(2, 0) * (X(1, 2) * X(0, 1) - X(1, 1) * X(0, 2));
  const T X1_1 = X(0, 0), X1_2 = X(0, 1), X1_3 = X(0, 2), X1_4 = X(0, 3);
  const T X2_1 = X(1, 0), X2_2 = X(1, 1), X2_3 = X(1, 2), X2_4 = X(1, 3);
  const T X3_1 = X(2, 0), X3_2 = X(2, 1), X3_3 = X(2, 2), X3_4 = X(2, 3);
  const T X4_1 = X(3, 0), X4_2 = X(3, 1), X4_3 = X(3, 2), X4_4 = X(3, 3);
  if (N == 4) {
    return X1_1 * X2_2 * X3_3 * X4_4 - X1_1 * X2_2 * X3_4 * X4_3 -
           X1_1 * X2_3 * X3_2 * X4_4 + X1_1 * X2_3 * X3_4 * X4_2 +
           X1_1 * X2_4 * X3_2 * X4_3 - X1_1 * X2_4 * X3_3 * X4_2 -
           X1_2 * X2_1 * X3_3 * X4_4 + X1_2 * X2_1 * X3_4 * X4_3 +
           X1_2 * X2_3 * X3_1 * X4_4 - X1_2 * X2_3 * X3_4 * X4_1 -
           X1_2 * X2_4 * X3_1 * X4_3 + X1_2 * X2_4 * X3_3 * X4_1 +
           X1_3 * X2_1 * X3_2 * X4_4 - X1_3 * X2_1 * X3_4 * X4_2 -
           X1_3 * X2_2 * X3_1 * X4_4 + X1_3 * X2_2 * X3_4 * X4_1 +
           X1_3 * X2_4 * X3_1 * X4_2 - X1_3 * X2_4 * X3_2 * X4_1 -
           X1_4 * X2_1 * X3_2 * X4_3 + X1_4 * X2_1 * X3_3 * X4_2 +
           X1_4 * X2_2 * X3_1 * X4_3 - X1_4 * X2_2 * X3_3 * X4_1 -
           X1_4 * X2_3 * X3_1 * X4_2 + X1_4 * X2_3 * X3_2 * X4_1;
  }
  std::cout << "unsupported size " << N << std::endl;
  NOT_POSSIBLE;
  return 0;
}

template <typename type>
symd<type> expm(const symd<type>& m) {
  std::pair<vecd<type>, matd<type>> decomp = m.eig();
  for (int k = 0; k < m.n(); k++) decomp.first(k) = ::exp(decomp.first(k));
  return symd<type>(decomp);
}

template <typename type>
symd<type> logm(const symd<type>& m) {
  std::pair<vecd<type>, matd<type>> decomp = m.eig();
  for (int k = 0; k < m.n(); k++) decomp.first(k) = ::log(decomp.first(k));
  return symd<type>(decomp);
}

template <typename type>
symd<type> powm(const symd<type>& m, double p) {
  std::pair<vecd<type>, matd<type>> decomp = m.eig();
  for (int k = 0; k < m.n(); k++) decomp.first(k) = ::pow(decomp.first(k), p);
  return symd<type>(decomp);
}

template <typename type>
symd<type> sqrtm(const symd<type>& m) {
  std::pair<vecd<type>, matd<type>> decomp = m.eig();
  for (int k = 0; k < m.n(); k++) decomp.first(k) = ::sqrt(decomp.first(k));
  return symd<type>(decomp);
}

template <int N, typename T>
syms<N, T> expm(const syms<N, T>& m) {
  std::pair<vecs<N, T>, mats<N, N, T>> decomp = m.eig();
  for (int k = 0; k < N; k++) decomp.first(k) = ::exp(decomp.first(k));
  return syms<N, T>(decomp);
}

template <int N, typename T>
syms<N, T> logm(const syms<N, T>& m) {
  std::pair<vecs<N, T>, mats<N, N, T>> decomp = m.eig();
  for (int k = 0; k < N; k++) decomp.first(k) = ::log(decomp.first(k));
  return syms<N, T>(decomp);
}

template <int N, typename T>
syms<N, T> powm(const syms<N, T>& m, double p) {
  std::pair<vecs<N, T>, mats<N, N, T>> decomp = m.eig();
  for (int k = 0; k < N; k++) decomp.first(k) = ::pow(decomp.first(k), p);
  return syms<N, T>(decomp);
}

template <int N, typename T>
syms<N, T> sqrtm(const syms<N, T>& m) {
  std::pair<vecs<N, T>, mats<N, N, T>> decomp = m.eig();
  for (int k = 0; k < N; k++) decomp.first(k) = ::sqrt(decomp.first(k));
  return syms<N, T>(decomp);
}

template void solveLUP(const matd<double>&, const vecd<double>&, vecd<double>&);
template void inverseLUP(const matd<double>&, matd<double>&);

#define INSTANTIATE_TRANSPOSE(T) template matd<T> transpose(const matd<T>&);
INSTANTIATE_TRANSPOSE(double)
INSTANTIATE_TRANSPOSE(float)
#undef INSTANTIATE_TRANSPOSE

#define INSTANTIATE_DIAG(T) template matd<T> diag(const vecd<T>&);
INSTANTIATE_DIAG(double)
INSTANTIATE_DIAG(float)
#undef INSTANTIATE_DIAG

#define INSTANTIATE_DIAGS(N, T) template mats<N, N, T> diag(const vecs<N, T>&);
INSTANTIATE_DIAGS(2, double)
INSTANTIATE_DIAGS(3, double)
#undef INSTANTIATE_DIAG

#define INSTANTIATE_INV(T)                  \
  template matd<T> inverse(const matd<T>&); \
  template symd<T> inverse(const symd<T>&);
INSTANTIATE_INV(double)
INSTANTIATE_INV(float)
#undef INSTANTIATE_INV

#define INSTANTIATE_EIG(T)                                     \
  template void eig(const symd<T>& m, vecd<T>& L, matd<T>& Q); \
  template std::pair<vecd<T>, matd<T>> eig(const symd<T>& m);
INSTANTIATE_EIG(double)
#undef INSTANTIATE_EIG

#define INSTANTIATE_EIGS(N, T)                                             \
  template void eig(const syms<N, T>& m, vecs<N, T>& L, mats<N, N, T>& Q); \
  template std::pair<vecs<N, T>, mats<N, N, T>> eig(const syms<N, T>& m);
INSTANTIATE_EIGS(2, double)
INSTANTIATE_EIGS(3, double)
#undef INSTANTIATE_EIGS

template int det(const matd<int>& X);
template double det(const matd<double>& X);
template double det(const symd<double>& X);

#define INSTANTIATE_SYMD_FUNC(T)                      \
  template symd<T> logm(const symd<T>&);              \
  template symd<T> expm(const symd<T>&);              \
  template symd<T> sqrtm(const symd<T>&);             \
  template symd<T> powm(const symd<T>&, double);      \
  template symd<T> interp(const std::vector<double>&, \
                          const std::vector<symd<T>>&);
INSTANTIATE_SYMD_FUNC(double)
INSTANTIATE_SYMD_FUNC(float)
#undef INSTANTIATE_SYMD_FUNC

#define INSTANTIATE_SYMS_FUNC(N, T)                      \
  template syms<N, T> logm(const syms<N, T>&);           \
  template syms<N, T> expm(const syms<N, T>&);           \
  template syms<N, T> sqrtm(const syms<N, T>&);          \
  template syms<N, T> powm(const syms<N, T>&, double);   \
  template syms<N, T> interp(const std::vector<double>&, \
                             const std::vector<syms<N, T>>&);
INSTANTIATE_SYMS_FUNC(2, double)
INSTANTIATE_SYMS_FUNC(3, double)
#undef INSTANTIATE_SYMS_FUNC

template double det(const mats<3, 3, double>&);
template double det(const mats<4, 4, double>&);
template float det(const mats<4, 4, float>&);

template mats<2, 2, double> inverse(const mats<2, 2, double>&);
template mats<3, 3, double> inverse(const mats<3, 3, double>&);
template mats<4, 4, double> inverse(const mats<4, 4, double>&);

}  // namespace vortex
