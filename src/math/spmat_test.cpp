#include "spmat.h"

#include "tester.h"
#include "vec.h"

using namespace vortex;

UT_TEST_SUITE(spmat_test_suite)

UT_TEST_CASE(test1) {
  spmat<double> A(2, 2);
  A(0, 0) = 1.0;
  A(0, 1) = 2.0;

  A(1, 0) = 3.0;
  A(1, 1) = 4.0;

  vecd<double> b(2);
  b(0) = 5.0;
  b(1) = 6.0;

  UT_ASSERT_EQUALS(A.nb_rows(), 2);
  UT_ASSERT_EQUALS(A.nb_cols(), 2);

  UT_ASSERT_EQUALS(A.nb_nnz(), 4);

  vecd<double> x(2);
  A.solve_nl(b, x);

  x.print();

  double tol = 1e-9;

  double X = x(0);
  double Y = x(1);

  UT_ASSERT_NEAR(X, -4.0, tol);
  UT_ASSERT_NEAR(Y, 4.5, tol);
}
UT_TEST_CASE_END(test1)

#define RHS_SINE 0  // we cannot represent the sin( pi * x) solution exactly

void get_poisson1d_system(int n, spmat<double>& A, vecd<double>& f) {
  ASSERT(f.m() == n + 1);
  double h = 1.0 / n;

  A(0, 0) = 1;
  f(0) = 0;
  for (int i = 1; i < n; i++) {
    A(i, i - 1) = -1 / (h * h);
    A(i, i) = 2 / (h * h);
    A(i, i + 1) = -1 / (h * h);

#if RHS_SINE
    double x = double(i) / n;
    f(i) = M_PI * M_PI * sin(M_PI * x);
#else
    f(i) = 1;
#endif
  }
  A(n, n) = 1.0;
  f(n) = 0.0;
}

UT_TEST_CASE(jacobi_test) {
  int n = 1e1;
  spmat<double> A(n + 1, n + 1);
  vecd<double> f(n + 1);

  get_poisson1d_system(n, A, f);
  A.print_full();

  double tol = 1e-10;

  vecd<double> x(n + 1);
  x.zero();

  A.solve_jacobi(f, x, tol, 1e4, true);

  double error = 0.0;
  for (int i = 0; i < n + 1; i++) {
    double xi = double(i) / n;
#if RHS_SINE
    double xa = sin(M_PI * xi);
#else
    double xa = 0.5 * xi * (1.0 - xi);
#endif
    error += (xa - x(i)) * (xa - x(i));
  }
  error = std::sqrt(error);
  LOG << "error = " << error;
}
UT_TEST_CASE_END(jacobi_test)

UT_TEST_CASE(gauss_siedel_test) {
  int n = 1e1;
  spmat<double> A(n + 1, n + 1);
  vecd<double> f(n + 1);

  get_poisson1d_system(n, A, f);
  A.print_full();

  double tol = 1e-10;

  vecd<double> x(n + 1);
  x.zero();

  int iter = 0;
  double e = norm(A * x - f);
  while (e > tol && iter++ < 1e4) {
    // perform the gauss-siedel step for each row
    for (int i = 1; i < n; i++) {
      double sigma = A(i, i - 1) * x(i - 1) + A(i, i + 1) * x(i + 1);
      x(i) = (f(i) - sigma) / A(i, i);
    }

    // recompute and log the error
    e = norm(A * x - f);
    LOG << fmt::format("iter {}, e = {}", iter, e);
  }
  LOG << fmt::format("converged to {} in {} iterations\n", e, iter);

  double error = 0.0;
  for (int i = 0; i < n + 1; i++) {
    double xi = double(i) / n;
#if RHS_SINE
    double xa = sin(M_PI * xi);
#else
    double xa = 0.5 * xi * (1.0 - xi);
#endif
    error += (xa - x(i)) * (xa - x(i));
  }
  error = std::sqrt(error);
  std::cout << "error = " << error << std::endl;
}
UT_TEST_CASE_END(gauss_siedel_test)

UT_TEST_CASE(cg_test) {
  int n = 1e2;
  spmat<double> A(n + 1, n + 1);
  vecd<double> b(n + 1);

  // get the 1d poisson system
  get_poisson1d_system(n, A, b);

  // solve the system
  vecd<double> x(n + 1);
  x.zero();
  A.solve_nl(b, x, false);
  x.print();

  // check the error
  double error = 0.0;
  for (int i = 0; i < n + 1; i++) {
    double xi = double(i) / n;
#if RHS_SINE
    double xa = sin(M_PI * double(i) / n);
#else
    double xa = 0.5 * xi * (1.0 - xi);
#endif
    error += (xa - x(i)) * (xa - x(i));
  }
  error = std::sqrt(error);
  std::cout << "error = " << error << std::endl;
  UT_ASSERT(error < 1e-10);
}
UT_TEST_CASE_END(cg_test)

UT_TEST_SUITE_END(spmat_test_suite)
