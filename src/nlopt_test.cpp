#include <fmt/format.h>

#if VORTEX_WITH_NLOPT
#include <nlopt.hpp>
#endif

#include "log.h"
#include "tester.h"

UT_TEST_SUITE(nlopt_test_suite)

#if VORTEX_WITH_NLOPT
struct nlopt_data {
  double k;
  double p;
  int iter{0};
};

// for Rosenbrock function (2D)
struct rosen_data {
  double a, b;
  int iter{0};
};

// for Rosenbrock function (3D)
struct rosen3d_data {
  double a, b;
  int iter{0};
};

// for booth function
struct booth_data {
  int iter{0};
};

// for circle constraint
struct circ_data {
  double h, k, r;
};

// for sphere constraint
struct sphere_data {
  double a, b, c, r;
};

// Energy function: F(x, y) = k(x^p + y^p)
double energy(unsigned n, const double *w, double *de_dw, void *data0) {
  // convert to the data type we need
  nlopt_data &data = *static_cast<nlopt_data *>(data0);

  LOG << fmt::format("iter. {}, guessed x = {}, y = {}", data.iter, w[0], w[1]);
  data.iter++;

  // if the algorithm requires the gradient of the function (here, energy)
  if (de_dw) {
    de_dw[0] = data.k * data.p * std::pow(w[0], data.p - 1);
    de_dw[1] = data.k * data.p * std::pow(w[1], data.p - 1);
  }

  // return the value of the function
  return data.k * (std::pow(w[0], data.p) + std::pow(w[1], data.p));
}

// Rosenbrock function: 2D --> F(x, y) = (a - x)^2 + b(y - x^2)^2
double rosen(unsigned n, const double *X, double *grad, void *rd) {
  // data type conversion
  rosen_data &data_r = *static_cast<rosen_data *>(rd);

  // record current guess
  LOG << fmt::format("iter. {}, guessed ({}, {})", data_r.iter, X[0], X[1]);
  data_r.iter++;

  // if function gradient is needed
  if (grad) {
    grad[0] =
        -2 * (data_r.a - X[0]) - (4 * X[0] * data_r.b * (X[1] - pow(X[0], 2)));
    grad[1] = 2 * data_r.b * (X[1] - pow(X[0], 2));
  }

  // return value of function
  return pow((data_r.a - X[0]), 2) + (data_r.b * pow((X[1] - pow(X[0], 2)), 2));
}

// Rosenbrock function: 3D --> F(X) = Σ((a - x_i)^2 + b(x_i+1 - x_i^2)^2)
double rosen3d(unsigned n, const double *X, double *grad, void *r3d) {
  // data type conversion
  rosen3d_data &data_r3d = *static_cast<rosen3d_data *>(r3d);

  // record current guess
  LOG << fmt::format("iter. {}, guessed ({}, {}, {})", data_r3d.iter, X[0],
                     X[1], X[2]);
  data_r3d.iter++;

  // if function gradient is needed
  if (grad) {
    grad[0] = (-2 * (data_r3d.a - X[0])) -
              (4 * data_r3d.b * X[0] * (X[1] - pow(X[0], 2)));
    grad[1] = (2 * data_r3d.b * (X[1] - pow(X[0], 2))) -
              (2 * (data_r3d.a - X[1])) -
              (4 * data_r3d.b * X[1] * (X[2] - pow(X[1], 2)));
    grad[2] = (2 * data_r3d.b * (X[2] - pow(X[1], 2)));
  }

  return pow(data_r3d.a - X[0], 2) +
         (data_r3d.b * pow(X[1] - pow(X[0], 2), 2)) +
         pow(data_r3d.a - X[1], 2) + (data_r3d.b * pow(X[2] - pow(X[1], 2), 2));
}

// return val of 3D Rosenbrock function at specified coords (check that
// constants are consistent with rosen3d) -- used for discretized min check
double rosen3d_val(double x, double y, double z) {
  return pow((1 - x), 2) + 100 * pow(y - pow(x, 2), 2) + pow(1 - y, 2) +
         100 * (pow(z - pow(y, 2), 2));
}

// Booth function: F(x, y) = (x + 2y - 7)^2 + (2x + y - 5)^2
double booth(unsigned n, const double *X, double *grad, void *bd) {
  // data type conversion
  booth_data &data_b = *static_cast<booth_data *>(bd);

  // record current guess
  LOG << fmt::format("iter. {}, guessed ({}, {})", data_b.iter, X[0], X[1]);
  data_b.iter++;

  // if function gradient is needed
  if (grad) {
    grad[0] = (2 * (X[0] + (2 * X[1]) - 7)) + (4 * ((2 * X[0]) + X[1] - 5));
    grad[1] = (4 * (X[0] + (2 * X[1]) - 7)) + (2 * ((2 * X[0]) + X[1] - 5));
  }

  // return value of function
  return pow(X[0] + (2 * X[1]) - 7, 2) + pow((2 * X[0]) + X[1] - 5, 2);
}

// Circle constraint (used for Rosenbrock and Booth function optimizations)
double circ_constr(unsigned n, const double *X, double *grad, void *cd) {
  // data type conversion
  circ_data &data_c = *static_cast<circ_data *>(cd);

  // if constraint gradient is needed
  if (grad) {
    grad[0] = 2 * (X[0] - data_c.h);
    grad[1] = 2 * (X[1] - data_c.k);
  }

  // return value of constraint
  return pow((X[0] - data_c.h), 2) + pow((X[1] - data_c.k), 2) -
         pow(data_c.r, 2);
}

// Sphere constraint
double sphere_constr(unsigned n, const double *X, double *grad, void *sd) {
  // data type conversion
  sphere_data &data_s = *static_cast<sphere_data *>(sd);

  // if constraint gradient is needed
  if (grad) {
    grad[0] = 2 * (X[0] - data_s.a);
    grad[1] = 2 * (X[1] - data_s.b);
    grad[2] = 2 * (X[2] - data_s.c);
  }

  // return vale of constraint
  return pow(X[0] - data_s.a, 2) + pow(X[1] - data_s.b, 2) +
         pow(X[2] - data_s.c, 2) - pow(data_s.r, 2);
}

// Optimizing Energy function (no constraint)
UT_TEST_CASE(test1) {
  int n = 2;  // two variables to optimize (x, y)
  nlopt::opt opt(nlopt::LD_LBFGS, n);

  // initialize the data used in the function and set the function to minimize
  nlopt_data data = {10, 2, 0};
  opt.set_min_objective(&energy, static_cast<void *>(&data));

  // set some optimization parameters
  opt.set_xtol_rel(1e-12);
  opt.set_ftol_rel(1e-12);
  opt.set_maxeval(100);

  // set the lower and upper bounds on the weights
  std::vector<double> lower_bound(n, -10.0);
  std::vector<double> upper_bound(n, 10.0);
  opt.set_lower_bounds(lower_bound);
  opt.set_upper_bounds(upper_bound);

  double f_opt;
  std::vector<double> x = {1, 1};  // initial guess
  nlopt::result result = opt.optimize(x, f_opt);
  LOG << fmt::format(
      "Optimized min (LD_LBFGS): ({}, {}), Energy Func. val. = {}", x[0], x[1],
      f_opt);
  LOG << fmt::format("Iterations: {}", data.iter);

  UT_ASSERT_EQUALS(result, nlopt::SUCCESS);
  UT_ASSERT_NEAR(x[0], 0, 1e-10);
  UT_ASSERT_NEAR(x[1], 0, 1e-10);
}
UT_TEST_CASE_END(test1)

// Optimizing Rosenbrock function (w/ a circle equality constraint)
UT_TEST_CASE(test2) {
  int n = 2;  // two variables to optimize (x, y)
  nlopt::opt opt_rosen(nlopt::LD_SLSQP, n);

  // initialize data used in the function, set opt_rosen to minimize
  rosen_data data_r = {1, 100, 0};
  circ_data data_c = {3, 4, 2};
  opt_rosen.set_min_objective(&rosen, static_cast<void *>(&data_r));

  // add the circle constraint
  opt_rosen.add_equality_constraint(&circ_constr, &data_c, 1e-8);

  // set the lower and upper bounds
  std::vector<double> lower_bound(n, -10.0);
  std::vector<double> upper_bound(n, 10.0);
  opt_rosen.set_lower_bounds(lower_bound);
  opt_rosen.set_upper_bounds(upper_bound);

  // optimization parameters
  opt_rosen.set_xtol_rel(1e-12);
  opt_rosen.set_ftol_rel(1e-12);
  opt_rosen.set_maxeval(100);

  double rosen_opt;
  std::vector<double> X = {0, 0};  // initial guess
  nlopt::result result = opt_rosen.optimize(X, rosen_opt);

  // discretized min check
  int N = 1000000;
  double min = std::numeric_limits<double>::max();
  double delta_theta = (2 * M_PI) / N;
  double x, y;
  double xm, ym;
  for (int i = 0; i < N; i++) {
    x = data_c.h + (data_c.r * cos(i * delta_theta));
    y = data_c.k + (data_c.r * sin(i * delta_theta));
    double curr = pow((data_r.a - x), 2) + (data_r.b * pow((y - pow(x, 2)), 2));
    if (curr < min) {
      xm = x;
      ym = y;
      min = curr;
    }
  }
  LOG << fmt::format("Optimized min (LD_SLSQP): ({}, {}), Rosenbrock val. = {}",
                     X[0], X[1], rosen_opt);
  LOG << fmt::format("Discretized min: ({}, {}), Rosenbrock val. = {}", xm, ym,
                     min);
  LOG << fmt::format("Iterations: {}", data_r.iter);

  // check constraint equals 0 at min
  UT_ASSERT_NEAR(
      0., pow(X[0] - data_c.h, 2) + pow(X[1] - data_c.k, 2) - pow(data_c.r, 2),
      1e-8);

  // check that the minimum val and x1, x2 coords on the constraint is roughly
  // equal to the value returned from optimization
  UT_ASSERT_NEAR(min, rosen_opt, 1e-8);
  UT_ASSERT_NEAR(X[0], xm, 1e-4);
  UT_ASSERT_NEAR(X[1], ym, 1e-4);

  // UT_ASSERT_EQUALS(result, nlopt::SUCCESS);
}
UT_TEST_CASE_END(test2)

// Optimizing a 3D Rosenbrock function, constrained on a sphere
UT_TEST_CASE(test3) {
  int n = 3;  // three variables to optimize (x1, x2, x3)
  nlopt::opt opt_rosen3d(nlopt::LD_SLSQP, n);

  // initialize the data used in the function, set opt_rosen3d to minimize
  rosen3d_data data_r3d = {1, 100, 0};
  sphere_data data_s = {0, 0, 1e-15, 1};  // unit sphere
  opt_rosen3d.set_min_objective(&rosen3d, static_cast<void *>(&data_r3d));

  // add sphere constraint
  opt_rosen3d.add_equality_constraint(&sphere_constr, &data_s, 1e-8);

  // set the lower and upper bounds
  std::vector<double> lower_bound(n, -10.0);
  std::vector<double> upper_bound(n, 10.0);
  opt_rosen3d.set_lower_bounds(lower_bound);
  opt_rosen3d.set_upper_bounds(upper_bound);

  // optimization parameters
  opt_rosen3d.set_xtol_rel(1e-12);
  opt_rosen3d.set_ftol_rel(1e-12);
  opt_rosen3d.set_maxeval(100);

  double rosen3d_opt;
  std::vector<double> X = {0, 0, 0};  // initial guess
  nlopt::result result = opt_rosen3d.optimize(X, rosen3d_opt);

  // discretized min check
  int N = 5000;
  double min = std::numeric_limits<double>::max();
  double delta_theta = (2 * M_PI) / N;
  double delta_phi = M_PI / N;
  double x, y, z;
  double xm, ym, zm;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      x = data_s.a + data_s.r * cos(j * delta_theta) * sin(i * delta_phi);
      y = data_s.b + data_s.r * sin(j * delta_theta) * sin(i * delta_phi);
      z = data_s.c + data_s.r * cos(i * delta_phi);
      double curr = rosen3d_val(x, y, z);
      if (curr < min) {
        min = curr;
        xm = x;
        ym = y;
        zm = z;
      }
    }
  }

  LOG << fmt::format(
      "Optimized min (LD_SLSQP): ({}, {}, {}), Rosenbrock (3D) val. = {}", X[0],
      X[1], X[2], rosen3d_opt);
  LOG << fmt::format("Discretized min: ({}, {}, {}), Rosenbrock (3D) val. = {}",
                     xm, ym, zm, min);
  LOG << fmt::format("Iterations: {}", data_r3d.iter);

  // check constraint equals 0 at min
  UT_ASSERT_NEAR(0.,
                 pow(X[0] - data_s.a, 2) + pow(X[1] - data_s.b, 2) +
                     pow(X[2] - data_s.c, 2) - pow(data_s.r, 2),
                 1e-8);

  // check that the minimum val and x1, x2, x3 coords on the constraint is
  // roughly equal to the value returned from optimization
  UT_ASSERT_NEAR(min, rosen3d_opt, 1e-4);
  UT_ASSERT_NEAR(X[0], xm, 1e-3);
  UT_ASSERT_NEAR(X[1], ym, 1e-3);
  UT_ASSERT_NEAR(X[2], zm, 1e-3);

  // UT_ASSERT_EQUALS(result, nlopt::SUCCESS);
}
UT_TEST_CASE_END(test3)

// Optimizing the Booth function (w/ or w/out a circle constraint)
UT_TEST_CASE(test4) {
  int n = 2;                                 // two variables to optimize (x, y)
  nlopt::opt opt_booth(nlopt::LD_LBFGS, n);  // defaults to unconstrained

  // initialize the data used in the function
  booth_data data_b = {0};

  // set opt_booth to minimize
  opt_booth.set_min_objective(&booth, static_cast<void *>(&data_b));

  // set the lower and upper bounds
  std::vector<double> lower_bound(n, -10.0);
  std::vector<double> upper_bound(n, 10.0);
  opt_booth.set_lower_bounds(lower_bound);
  opt_booth.set_upper_bounds(upper_bound);

  // optimization parameters
  opt_booth.set_xtol_rel(1e-12);
  opt_booth.set_ftol_rel(1e-12);
  opt_booth.set_maxeval(100);

  double booth_opt;
  std::vector<double> X = {0, 0};  // initial guess
  nlopt::result result = opt_booth.optimize(X, booth_opt);

  LOG << fmt::format("Optimized min (LD_LBFGS): ({}, {}), Booth val. = {}",
                     X[0], X[1], booth_opt);
  UT_ASSERT_NEAR(X[0], 1, 1e-8);
  UT_ASSERT_NEAR(X[1], 3, 1e-8);

  LOG << fmt::format("Iterations: {}", data_b.iter);

  UT_ASSERT_EQUALS(result, nlopt::SUCCESS);
}
UT_TEST_CASE_END(test4)
#else
UT_TEST_CASE(test1) {}
UT_TEST_CASE_END(test1)
#endif
UT_TEST_SUITE_END(nlopt_test_suite)