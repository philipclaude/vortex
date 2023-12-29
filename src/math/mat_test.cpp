#include "mat.hpp"

#include "linalg.h"
#include "tester.h"
#include "vec.h"
#include "vec.hpp"

using namespace vortex;

UT_TEST_SUITE(matrix_test_suite)

template <typename T>
class randmatd : public matd<T> {
 public:
  randmatd(int m, int n) : matd<T>(m, n) {
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++) (*this)(i, j) = T(rand()) / T(RAND_MAX);
  }
};

template <typename T>
class randvecd : public vecd<T> {
 public:
  randvecd(int m) : vecd<T>(m) {
    for (int i = 0; i < m; i++) (*this)(i) = T(rand()) / T(RAND_MAX);
  }
};

template <int M, int N, typename T>
class randmats : public mats<M, N, T> {
 public:
  randmats() : mats<M, N, T>() {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++) (*this)(i, j) = T(rand()) / T(RAND_MAX);
  }
};

UT_TEST_CASE(matd_test) {
  double dtol = 1e-12;
  double ftol = 1e-7;

  matd<double> X;
  UT_ASSERT_EQUALS(X.m(), 0);
  UT_ASSERT_EQUALS(X.n(), 0);

  X.resize(2, 3);
  UT_ASSERT_EQUALS(X.m(), 2);
  UT_ASSERT_EQUALS(X.n(), 3);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 3; j++) UT_ASSERT_EQUALS(X(i, j), 0);

  // constructor for a square matrix
  matd<double> m1(3);
  UT_ASSERT_EQUALS(m1.m(), 3);
  UT_ASSERT_EQUALS(m1.n(), 3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) UT_ASSERT_EQUALS(m1(i, j), 0);

  // constructor for a rectangular matrix
  matd<double> m2(3, 5);
  UT_ASSERT_EQUALS(m2.m(), 3);
  UT_ASSERT_EQUALS(m2.n(), 5);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 5; j++) UT_ASSERT_EQUALS(m2(i, j), 0);

  // constructor from another matrix
  m2(1, 2) = 1.2;
  matd<double> m3(m2);
  UT_ASSERT_EQUALS(m3.m(), 3);
  UT_ASSERT_EQUALS(m3.n(), 5);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 5; j++) UT_ASSERT_NEAR(m3(i, j), m2(i, j), dtol);

  matd<double> m4 = m2;
  UT_ASSERT_EQUALS(m4.m(), 3);
  UT_ASSERT_EQUALS(m4.n(), 5);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 5; j++) UT_ASSERT_NEAR(m4(i, j), m2(i, j), dtol);

  matd<float> m5 = m2;
  UT_ASSERT_EQUALS(m5.m(), 3);
  UT_ASSERT_EQUALS(m5.n(), 5);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 5; j++)
      UT_ASSERT_NEAR(m5(i, j), m2(i, j), ftol);  // m5 is a matrix of floats

  // operators
  randmatd<double> A(2, 3), B(3, 4);

  matd<double> C = A * B;
  UT_ASSERT_EQUALS(C.m(), 2);
  UT_ASSERT_EQUALS(C.n(), 4);

  (A * B).print("AB");
  C.print("C");

  UT_ASSERT_NEAR(
      C(0, 0), A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0), dtol);
  UT_ASSERT_NEAR(
      C(0, 1), A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1), dtol);
  UT_ASSERT_NEAR(
      C(0, 2), A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2), dtol);
  UT_ASSERT_NEAR(
      C(0, 3), A(0, 0) * B(0, 3) + A(0, 1) * B(1, 3) + A(0, 2) * B(2, 3), dtol);

  UT_ASSERT_NEAR(
      C(1, 0), A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0), dtol);
  UT_ASSERT_NEAR(
      C(1, 1), A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1), dtol);
  UT_ASSERT_NEAR(
      C(1, 2), A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2), dtol);
  UT_ASSERT_NEAR(
      C(1, 3), A(1, 0) * B(0, 3) + A(1, 1) * B(1, 3) + A(1, 2) * B(2, 3), dtol);

  C.print("C");

  matd<double> C0 = C;
  UT_ASSERT_EQUALS(C0.m(), 2);
  UT_ASSERT_EQUALS(C0.n(), 4);

  randmatd<double> D(2, 4);
  D.print("D");
  (C0 + D).print("C0+D");
  C = C0 + D;
  C.print("C");
  UT_ASSERT_EQUALS(C.m(), 2);
  UT_ASSERT_EQUALS(C.n(), 4);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 4; j++) {
      UT_ASSERT_NEAR(C(i, j), C0(i, j) + D(i, j), dtol);
    }

  C = +C0 + D - D;
  UT_ASSERT_EQUALS(C.m(), 2);
  UT_ASSERT_EQUALS(C.n(), 4);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 4; j++) UT_ASSERT_NEAR(C(i, j), C0(i, j), dtol);

  C = -C;
  UT_ASSERT_EQUALS(C.m(), 2);
  UT_ASSERT_EQUALS(C.n(), 4);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 4; j++) UT_ASSERT_NEAR(C(i, j), -C0(i, j), dtol);

  C.zero();
  UT_ASSERT_EQUALS(C.m(), 2);
  UT_ASSERT_EQUALS(C.n(), 4);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 4; j++) UT_ASSERT_EQUALS(C(i, j), 0.0);

  C.print();
}
UT_TEST_CASE_END(matd_test)

UT_TEST_CASE(mats_test) {
  double dtol = 1e-15;
  double ftol = 1e-7;

  // constructor for a square matrix
  mats<3, 3, double> m1;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) UT_ASSERT_EQUALS(m1(i, j), 0);

  // constructor for a rectangular matrix
  mats<3, 5, double> m2;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 5; j++) UT_ASSERT_EQUALS(m2(i, j), 0);

  // constructor from another matrix
  m2(1, 2) = 1.2;
  mats<3, 5, double> m3(m2);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 5; j++) UT_ASSERT_NEAR(m3(i, j), m2(i, j), dtol);

  mats<3, 5, double> m4 = m2;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 5; j++) UT_ASSERT_NEAR(m4(i, j), m2(i, j), dtol);

  mats<3, 5, float> m5 = m2;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 5; j++)
      UT_ASSERT_NEAR(m5(i, j), m2(i, j), ftol);  // m5 is a matrix of floats

  // operators
  randmats<2, 3, double> A;
  randmats<3, 4, double> B;

  mats<2, 4, double> C = A * B;

  UT_ASSERT_NEAR(
      C(0, 0), A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0) + A(0, 2) * B(2, 0), dtol);
  UT_ASSERT_NEAR(
      C(0, 1), A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1) + A(0, 2) * B(2, 1), dtol);
  UT_ASSERT_NEAR(
      C(0, 2), A(0, 0) * B(0, 2) + A(0, 1) * B(1, 2) + A(0, 2) * B(2, 2), dtol);
  UT_ASSERT_NEAR(
      C(0, 3), A(0, 0) * B(0, 3) + A(0, 1) * B(1, 3) + A(0, 2) * B(2, 3), dtol);

  UT_ASSERT_NEAR(
      C(1, 0), A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0) + A(1, 2) * B(2, 0), dtol);
  UT_ASSERT_NEAR(
      C(1, 1), A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1) + A(1, 2) * B(2, 1), dtol);
  UT_ASSERT_NEAR(
      C(1, 2), A(1, 0) * B(0, 2) + A(1, 1) * B(1, 2) + A(1, 2) * B(2, 2), dtol);
  UT_ASSERT_NEAR(
      C(1, 3), A(1, 0) * B(0, 3) + A(1, 1) * B(1, 3) + A(1, 2) * B(2, 3), dtol);

  mats<2, 4, double> C0 = C;

  randmats<2, 4, double> D;
  C = C + D;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 4; j++)
      UT_ASSERT_NEAR(C(i, j), C0(i, j) + D(i, j), dtol);

  C = +C0 + D - D;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 4; j++) UT_ASSERT_NEAR(C(i, j), C0(i, j), dtol);

  C = -C;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 4; j++) UT_ASSERT_NEAR(C(i, j), -C0(i, j), dtol);

  C.zero();
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 4; j++) UT_ASSERT_EQUALS(C(i, j), 0.0);

  C.print();
}
UT_TEST_CASE_END(mats_test)

UT_TEST_CASE(mvm_test) {
  // matrix-vector-multiplication tests
  double dtol = 1e-15;

  randmatd<double> A(3, 3);
  randvecd<double> b(3);

  b.print("b");

  vecd<double> c = A * b;

  UT_ASSERT_NEAR(c(0), A(0, 0) * b(0) + A(0, 1) * b(1) + A(0, 2) * b(2), dtol);
  UT_ASSERT_NEAR(c(1), A(1, 0) * b(0) + A(1, 1) * b(1) + A(1, 2) * b(2), dtol);
  UT_ASSERT_NEAR(c(2), A(2, 0) * b(0) + A(2, 1) * b(1) + A(2, 2) * b(2), dtol);

  randvecd<double> b2(5);
  UT_CATCH_EXCEPTION(A * b2);

  randmatd<double> A2(3, 5);
  vecd<double> c2 = A2 * b2;

  UT_ASSERT_EQUALS(c2.m(), 3);

  UT_ASSERT_NEAR(c2(0),
                 A2(0, 0) * b2(0) + A2(0, 1) * b2(1) + A2(0, 2) * b2(2) +
                     A2(0, 3) * b2(3) + A2(0, 4) * b2(4),
                 dtol);
  UT_ASSERT_NEAR(c2(1),
                 A2(1, 0) * b2(0) + A2(1, 1) * b2(1) + A2(1, 2) * b2(2) +
                     A2(1, 3) * b2(3) + A2(1, 4) * b2(4),
                 dtol);
  UT_ASSERT_NEAR(c2(2),
                 A2(2, 0) * b2(0) + A2(2, 1) * b2(1) + A2(2, 2) * b2(2) +
                     A2(2, 3) * b2(3) + A2(2, 4) * b2(4),
                 dtol);

  transpose(b).print();
  double x = transpose(b) * A * b;
  std::cout << "x = " << x << std::endl;
  UT_ASSERT_NEAR(x, b(0) * c(0) + b(1) * c(1) + b(2) * c(2), dtol);
}
UT_TEST_CASE_END(mvm_test)

UT_TEST_CASE(quadratic_form_test) {
  vec4d x = {1, 2, 3, 4};

  mat4d m;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) m(i, j) = 1.0;

  // m * x is {10,10,10,10}, so this should be 10 + 20 + 30 + 40
  double e = transpose(x) * m * x;

  UT_ASSERT_NEAR(e, 100.0, 1e-12);
}
UT_TEST_CASE_END(quadratic_form_test)

UT_TEST_SUITE_END(matrix_test_suite)
