#include "sym.h"

#include "linalg.h"
#include "mat.hpp"
#include "sym.hpp"
#include "tester.h"
#include "vec.h"

using namespace vortex;

UT_TEST_SUITE(symmetric_matrix_test_suite)

template <typename T>
class randsymd : public symd<T> {
 public:
  randsymd(int n) : symd<T>(n) {
    matd<double> X(n, n);
    matd<double> Xt(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
        X(i, j) = double(rand()) / double(RAND_MAX);
        Xt(j, i) = X(i, j);
      }

    symd<T> S = Xt * X;
    this->set(S);

    // add n to diagonal for spd-ness
    for (int i = 0; i < n; i++) (*this)(i, i) = S(i, i) + n;
  }
};

template <int m, typename T>
class randvecs : public vecs<m, T> {
 public:
  randvecs() : vecs<m, T>() {
    for (int i = 0; i < m; i++) (*this)(i) = T(rand()) / T(RAND_MAX);
  }
};

template <int n, typename T>
class randsyms : public syms<n, T> {
 public:
  randsyms() : syms<n, T>() {
    mats<n, n, T> X;
    mats<n, n, T> Xt;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
        X(i, j) = double(rand()) / double(RAND_MAX);
        Xt(j, i) = X(i, j);
      }

    syms<n, T> S = Xt * X;
    this->set(S);

    // add n to diagonal for spd-ness
    for (int i = 0; i < n; i++) (*this)(i, i) = S(i, i) + n;
  }
};

template <typename T>
class randvecd : public vecd<T> {
 public:
  randvecd(int m) : vecd<T>(m) {
    for (int i = 0; i < m; i++) (*this)(i) = T(rand()) / T(RAND_MAX);
  }
};

UT_TEST_CASE(symd_test) {
  double dtol = 1e-12;
  double ftol = 1e-5;

  symd<double> X;
  UT_ASSERT_EQUALS(X.m(), 0);
  UT_ASSERT_EQUALS(X.n(), 0);

  X.resize(3);
  UT_ASSERT_EQUALS(X.m(), 3);
  UT_ASSERT_EQUALS(X.n(), 3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) UT_ASSERT_EQUALS(X(i, j), 0);

  // constructor
  symd<double> m1(3);
  UT_ASSERT_EQUALS(m1.m(), 3);
  UT_ASSERT_EQUALS(m1.n(), 3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) UT_ASSERT_EQUALS(m1(i, j), 0);

  // constructor from another matrix
  m1(1, 2) = 1.2;
  symd<double> m2(m1);
  UT_ASSERT_EQUALS(m2.m(), 3);
  UT_ASSERT_EQUALS(m2.n(), 3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) UT_ASSERT_NEAR(m2(i, j), m1(i, j), dtol);

  symd<double> m3 = m1;
  UT_ASSERT_EQUALS(m3.m(), 3);
  UT_ASSERT_EQUALS(m3.n(), 3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) UT_ASSERT_NEAR(m3(i, j), m1(i, j), dtol);

  symd<float> m4 = m2;
  UT_ASSERT_EQUALS(m4.m(), 3);
  UT_ASSERT_EQUALS(m4.n(), 3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      UT_ASSERT_NEAR(m4(i, j), m2(i, j), ftol);  // m5 is a matrix of floats

  // operators
  randsymd<double> A(3), B(3);

  matd<double> C = A * B;
  UT_ASSERT_EQUALS(C.m(), 3);
  UT_ASSERT_EQUALS(C.n(), 3);

  (A * B).print("AB");
  C.print("C");

  UT_ASSERT_NEAR(
      C(0, 0), A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0), dtol);
  UT_ASSERT_NEAR(
      C(0, 1), A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1), dtol);
  UT_ASSERT_NEAR(
      C(0, 2), A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2), dtol);

  UT_ASSERT_NEAR(
      C(1, 0), A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0), dtol);
  UT_ASSERT_NEAR(
      C(1, 1), A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1), dtol);
  UT_ASSERT_NEAR(
      C(1, 2), A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2), dtol);

  UT_ASSERT_NEAR(
      C(2, 0), A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0), dtol);
  UT_ASSERT_NEAR(
      C(2, 1), A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1), dtol);
  UT_ASSERT_NEAR(
      C(2, 2), A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2), dtol);

  C.print("C");

  symd<double> C0 = C;
  UT_ASSERT_EQUALS(C0.m(), 3);
  UT_ASSERT_EQUALS(C0.n(), 3);
  C0.print();

  randsymd<double> D(3);
  C = C0 + D;
  C.print("C");
  UT_ASSERT_EQUALS(C.m(), 3);
  UT_ASSERT_EQUALS(C.n(), 3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      UT_ASSERT_NEAR(C(i, j), C0(i, j) + D(i, j), dtol);
    }

  C = +C0 + D - D;
  UT_ASSERT_EQUALS(C.m(), 3);
  UT_ASSERT_EQUALS(C.n(), 3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) UT_ASSERT_NEAR(C(i, j), C0(i, j), dtol);

  C = -C;
  UT_ASSERT_EQUALS(C.m(), 3);
  UT_ASSERT_EQUALS(C.n(), 3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) UT_ASSERT_NEAR(C(i, j), -C0(i, j), dtol);

  C.zero();
  UT_ASSERT_EQUALS(C.m(), 3);
  UT_ASSERT_EQUALS(C.n(), 3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) UT_ASSERT_EQUALS(C(i, j), 0.0);

  C.print();

  m1 = 2.0;
  m1.print();
  for (int i = 0; i < m1.nb(); i++) UT_ASSERT_NEAR(m1.data(i), 2.0, dtol);

  symd<double> Z(3);
  Z.eye();
  symd<double> Y = interp<double>({0.5, 0.5}, {Z, Z});
  for (int i = 0; i < 3; i++) {
    UT_ASSERT_NEAR(Y(i, i), 1.0, dtol);
    for (int j = i + 1; j < 3; j++) UT_ASSERT_NEAR(Y(i, j), 0.0, dtol);
  }
}
UT_TEST_CASE_END(symd_test)

UT_TEST_CASE(syms_test) {
  double dtol = 1e-12;

  {
    static const int N = 3;

    // constructor
    syms<N, double> m1;
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) UT_ASSERT_EQUALS(m1(i, j), 0);

    // constructor from another matrix
    m1(1, 2) = 1.2;
    syms<N, double> m2(m1);
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) UT_ASSERT_NEAR(m2(i, j), m1(i, j), dtol);

    syms<N, double> m3 = m1;
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) UT_ASSERT_NEAR(m3(i, j), m1(i, j), dtol);

    // operators
    randsyms<N, double> A, B;

    mats<N, N, double> C = A * B;

    (A * B).print();
    C.print();

    UT_ASSERT_NEAR(C(0, 0),
                   A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0),
                   dtol);
    UT_ASSERT_NEAR(C(0, 1),
                   A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1),
                   dtol);
    UT_ASSERT_NEAR(C(0, 2),
                   A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2),
                   dtol);

    UT_ASSERT_NEAR(C(1, 0),
                   A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0),
                   dtol);
    UT_ASSERT_NEAR(C(1, 1),
                   A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1),
                   dtol);
    UT_ASSERT_NEAR(C(1, 2),
                   A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2),
                   dtol);

    UT_ASSERT_NEAR(C(2, 0),
                   A(2, 0) * B(0, 0) + A(2, 1) * B(1, 0) + A(2, 2) * B(2, 0),
                   dtol);
    UT_ASSERT_NEAR(C(2, 1),
                   A(2, 0) * B(0, 1) + A(2, 1) * B(1, 1) + A(2, 2) * B(2, 1),
                   dtol);
    UT_ASSERT_NEAR(C(2, 2),
                   A(2, 0) * B(0, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(2, 2),
                   dtol);

    C.print();

    syms<N, double> C0 = C;
    C0.print();

    randsyms<N, double> D;
    C = C0 + D;
    C.print();
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) {
        UT_ASSERT_NEAR(C(i, j), C0(i, j) + D(i, j), dtol);
      }

    C = +C0 + D - D;
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) UT_ASSERT_NEAR(C(i, j), C0(i, j), dtol);

    C = -C;
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) UT_ASSERT_NEAR(C(i, j), -C0(i, j), dtol);

    C.zero();
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) UT_ASSERT_EQUALS(C(i, j), 0.0);

    C.print();

    C0 = 2.0;
    for (int i = 0; i < C0.nb(); i++) UT_ASSERT_NEAR(C0.data(i), 2.0, dtol);

    syms<N, double> Y = C;
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) UT_ASSERT_NEAR(Y(i, j), C(i, j), dtol);

    Y.eye();
    for (int i = 0; i < N; i++) {
      UT_ASSERT_NEAR(Y(i, i), 1.0, dtol);
      for (int j = i + 1; j < N; j++) UT_ASSERT_EQUALS(Y(i, j), 0.0);
    }

    syms<N, double> Z;
    Z.eye();
    Y = interp<N, double>({0.5, 0.5}, {Z, Z});
    for (int i = 0; i < 3; i++) {
      UT_ASSERT_NEAR(Y(i, i), 1.0, dtol);
      for (int j = i + 1; j < 3; j++) UT_ASSERT_NEAR(Y(i, j), 0.0, dtol);
    }
  }
}
UT_TEST_CASE_END(syms_test)

UT_TEST_SUITE_END(symmetric_matrix_test_suite)
