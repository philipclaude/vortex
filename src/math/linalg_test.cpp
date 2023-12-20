#include "linalg.h"

#include "mat.h"
#include "mat.hpp"
#include "sym.h"
#include "sym.hpp"
#include "tester.h"

using namespace vortex;

matd<double> random_matrix(int n) {
  matd<double> X(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) X(i, j) = double(rand()) / double(RAND_MAX);
  return X;
}

symd<double> random_tensor(int n) {
  matd<double> X(n, n);
  matd<double> Xt(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      X(i, j) = double(rand()) / double(RAND_MAX);
      Xt(j, i) = X(i, j);
    }

  symd<double> N = Xt * X;

  // add n to diagonal for spd-ness
  for (int i = 0; i < n; i++) N(i, i) += n;
  return N;
}

template <int N> syms<N, double> random_tensor() {
  matd<double> X(N, N);
  matd<double> Xt(N, N);
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) {
      X(i, j) = double(rand()) / double(RAND_MAX);
      Xt(j, i) = X(i, j);
    }

  symd<double> Y = Xt * X;

  // add n to diagonal for spd-ness
  for (int i = 0; i < N; i++) Y(i, i) += N;

  syms<N, double> Z;
  for (int i = 0; i < Y.nb(); i++) Z.data(i) = Y.data(i);
  return Z;
}

UT_TEST_SUITE(linear_algebra_tests)

UT_TEST_CASE(determinant_tests) {
  // hopefully the determinants are not zero when checking the inverse
  double tol = 1e-12;

  mats<2, 2, double> A2;
  A2(0, 0) = 1;
  A2(0, 1) = 2;
  A2(1, 0) = 3;
  A2(1, 1) = 4;
  UT_ASSERT_NEAR(det(A2), -2.0, tol);

  mats<2, 2, double> A2inv = inverse(A2);
  mats<2, 2, double> I2 = A2inv * A2;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      if (i == j)
        UT_ASSERT_NEAR(I2(i, j), 1.0, tol);
      else
        UT_ASSERT_NEAR(I2(i, j), 0.0, tol);
    }
  }

  mats<3, 3, double> A3;
  A3(0, 0) = 1;
  A3(0, 1) = 2;
  A3(0, 2) = 3;
  A3(1, 0) = 4;
  A3(1, 1) = 5;
  A3(1, 2) = 6;
  A3(2, 0) = 7;
  A3(2, 1) = 8;
  A3(2, 2) = 10;
  UT_ASSERT_NEAR(det(A3), -3., tol);

  mats<3, 3, double> A3inv = inverse(A3);
  mats<3, 3, double> I3 = A3inv * A3;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j)
        UT_ASSERT_NEAR(I3(i, j), 1.0, tol);
      else
        UT_ASSERT_NEAR(I3(i, j), 0.0, tol);
    }
  }

  mats<4, 4, double> A4;
  A4(0, 0) = 1;
  A4(0, 1) = 2;
  A4(0, 2) = 6;
  A4(0, 3) = 8;
  A4(1, 0) = 9;
  A4(1, 1) = 2;
  A4(1, 2) = 4;
  A4(1, 3) = 1;
  A4(2, 0) = 7;
  A4(2, 1) = 4;
  A4(2, 2) = 3;
  A4(2, 3) = 2;
  A4(3, 0) = 1;
  A4(3, 1) = 1;
  A4(3, 2) = 1;
  A4(3, 3) = 1;
  UT_ASSERT_NEAR(det(A4), -29., tol);

  mats<4, 4, double> A4inv = inverse(A4);
  mats<4, 4, double> I4 = A4inv * A4;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      if (i == j)
        UT_ASSERT_NEAR(I4(i, j), 1.0, tol);
      else
        UT_ASSERT_NEAR(I4(i, j), 0.0, tol);
    }
  }

  mats<5, 5, double> A5, A5inv;
  UT_CATCH_EXCEPTION(A5inv = inverse(A5));
}
UT_TEST_CASE_END(determinant_tests)

UT_TEST_CASE(inverse_tests) {
  double tol = 1e-10;
  int ntests = 10;

  for (int n = 1; n <= 4; n++) {
    matd<double> I(n, n);
    I.eye();

    for (int k = 0; k < ntests; k++) {
      matd<double> A(n, n);
      A = random_matrix(n);

      if (std::fabs(det(A) < 1e-12)) continue;

      matd<double> Ainv(n, n);
      Ainv = inverse(A);

      matd<double> B(n, n);
      B = A * Ainv;

      matd<double> zero = (A * Ainv - I);

      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) UT_ASSERT_NEAR(zero(i, j), 0.0, tol);

      symd<double> T = random_tensor(n);
      symd<double> Tinv = inverse(T);

      symd<double> Is(n);
      Is.eye();
      symd<double> zeros = T * Tinv - Is;

      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) UT_ASSERT_NEAR(zeros(i, j), 0.0, tol);

      zeros = Is - T * Tinv;
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) UT_ASSERT_NEAR(zeros(i, j), 0.0, tol);
    }

    symd<double> S5(5);
    UT_CATCH_EXCEPTION(inverse(S5));
    UT_CATCH_EXCEPTION(det(S5));
  }

  matd<double> A(5, 5), Ainv(5, 5), B(5, 4);
  UT_CATCH_EXCEPTION(Ainv = inverse(A));
  UT_CATCH_EXCEPTION(Ainv = inverse(B));

  // inverse LUP
  for (int n = 5; n <= 10; n++) {
    matd<double> I(n, n);
    I.eye();

    for (int k = 0; k < ntests; k++) {
      matd<double> A(n, n);
      A = random_matrix(n);

      matd<double> Ainv(n, n);
      inverseLUP(A, Ainv);

      matd<double> B(n, n);
      B = A * Ainv;

      matd<double> zero = (A * Ainv - I);

      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) UT_ASSERT_NEAR(zero(i, j), 0.0, tol);
    }
  }
}
UT_TEST_CASE_END(inverse_tests)

UT_TEST_CASE(solveLUP_tests) {
  double tol = 1e-12;
  int ntests = 10;

  // inverse LUP
  for (int n = 5; n <= 10; n++) {
    for (int k = 0; k < ntests; k++) {
      matd<double> A(n, n);
      A = random_matrix(n);

      vecd<double> b(n);
      for (int i = 0; i < n; i++) b[i] = double(rand()) / double(RAND_MAX);

      vecd<double> x(n);
      solveLUP(A, b, x);

      vecd<double> y = A * x;
      for (int i = 0; i < n; i++) UT_ASSERT_NEAR(y[i], b[i], tol);
    }
  }
}
UT_TEST_CASE_END(solveLUP_tests)

UT_TEST_CASE(symd_eign_test) {
  double tol = 1e-10;

  for (int n = 2; n < 6; n++) {
    symd<double> A = random_tensor(n);

    std::pair<vecd<double>, matd<double> > e = eig(A);
    const matd<double>& Q = e.second;
    const vecd<double>& L = e.first;

    symd<double> A1(n);
    A1.from_eig(L, Q);

    // A1.print();
    // A.print();

    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) UT_ASSERT_NEAR(A(i, j), A1(i, j), tol);

    matd<double> Q2(n, n);
    vecd<double> L2(n);
    eig(A, L2, Q2);

    UT_ASSERT_EQUALS(L2.m(), n);
    UT_ASSERT_EQUALS(Q2.m(), n);
    UT_ASSERT_EQUALS(Q2.n(), n);
    for (int i = 0; i < n; i++) {
      UT_ASSERT_NEAR(L2(i), L(i), tol);
      for (int j = 0; j < n; j++) {
        UT_ASSERT_NEAR(Q2(i, j), Q(i, j), tol);
      }
    }
  }
}
UT_TEST_CASE_END(symd_eign_test)

UT_TEST_CASE(symd_functions_test) {
  double tol = 1e-10;

  for (int n = 2; n <= 6; n++) {
    symd<double> A = random_tensor(n);

    std::pair<vecd<double>, matd<double> > e = eig(A);
    const matd<double>& Q = e.second;
    const vecd<double>& L0 = e.first;

    vecd<double> L(n);
    symd<double> X(n);

    // exponent
    for (int i = 0; i < n; i++) L(i) = std::exp(L0(i));
    X.from_eig(L, Q);

    symd<double> expA = expm(A);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) UT_ASSERT_NEAR(expA(i, j), X(i, j), tol);

    // logarithm
    symd<double> logexpA = logm(expA);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) UT_ASSERT_NEAR(logexpA(i, j), A(i, j), tol);

    // square root
    for (int i = 0; i < n; i++) L(i) = std::sqrt(L0(i));
    X.from_eig(L, Q);

    symd<double> sqrtA = sqrtm(A);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) UT_ASSERT_NEAR(sqrtA(i, j), X(i, j), tol);

    // power
    symd<double> squaresqrtA = powm(sqrtA, 2.0);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        UT_ASSERT_NEAR(squaresqrtA(i, j), A(i, j), tol);
  }
}
UT_TEST_CASE_END(symd_functions_test)

UT_TEST_CASE(syms_functions_test) {
  double tol = 1e-10;
  typedef double T;

  {
    static const int n = 2;

    syms<n, double> A = random_tensor<n>();

    std::pair<vecs<n, T>, mats<n, n, T> > e = eig(A);
    const mats<n, n, double>& Q = e.second;
    const vecs<n, double>& L0 = e.first;

    vecs<n, double> L;
    syms<n, double> X;

    // exponent
    for (int i = 0; i < n; i++) L(i) = std::exp(L0(i));
    X.from_eig(L, Q);

    syms<n, double> expA = expm(A);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) UT_ASSERT_NEAR(expA(i, j), X(i, j), tol);

    // logarithm
    syms<n, double> logexpA = logm(expA);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) UT_ASSERT_NEAR(logexpA(i, j), A(i, j), tol);

    // square root
    for (int i = 0; i < n; i++) L(i) = std::sqrt(L0(i));
    X.from_eig(L, Q);

    syms<n, double> sqrtA = sqrtm(A);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) UT_ASSERT_NEAR(sqrtA(i, j), X(i, j), tol);

    // power
    syms<n, double> squaresqrtA = powm(sqrtA, 2.0);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        UT_ASSERT_NEAR(squaresqrtA(i, j), A(i, j), tol);
  }

  {
    static const int n = 3;

    syms<n, double> A = random_tensor<n>();

    std::pair<vecs<n, T>, mats<n, n, T> > e = eig(A);
    const mats<n, n, double>& Q = e.second;
    const vecs<n, double>& L0 = e.first;

    vecs<n, double> L;
    syms<n, double> X;

    // exponent
    for (int i = 0; i < n; i++) L(i) = std::exp(L0(i));
    X.from_eig(L, Q);

    syms<n, double> expA = expm(A);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) UT_ASSERT_NEAR(expA(i, j), X(i, j), tol);

    // logarithm
    syms<n, double> logexpA = logm(expA);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) UT_ASSERT_NEAR(logexpA(i, j), A(i, j), tol);

    // square root
    for (int i = 0; i < n; i++) L(i) = std::sqrt(L0(i));
    X.from_eig(L, Q);

    syms<n, double> sqrtA = sqrtm(A);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) UT_ASSERT_NEAR(sqrtA(i, j), X(i, j), tol);

    // power
    syms<n, double> squaresqrtA = powm(sqrtA, 2.0);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        UT_ASSERT_NEAR(squaresqrtA(i, j), A(i, j), tol);
  }
}
UT_TEST_CASE_END(syms_functions_test)

UT_TEST_SUITE_END(linear_algebra_tests)
