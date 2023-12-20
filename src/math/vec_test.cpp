#include "vec.h"

#include "linear_algebra.h"
#include "tester.h"
#include "vec.hpp"

using namespace vortex;

UT_TEST_SUITE(vector_suite)

UT_TEST_CASE(vecd_test) {
  double tol = 1e-12;

  vecd<double> x({1, 2, 3});
  vecd<double> y({4, 5, 6});

  UT_ASSERT_NEAR(dot(x, x), 14., tol);
  UT_ASSERT_NEAR(length(x), std::sqrt(14), tol);

  x = x + y;

  UT_ASSERT_EQUALS(x[0], 5);
  UT_ASSERT_EQUALS(x[1], 7);
  UT_ASSERT_EQUALS(x[2], 9);

  UT_ASSERT_NEAR(dot(x, x), 155., tol);
  UT_ASSERT_NEAR(length(x), std::sqrt(155), tol);

  vecd<double> z = 3.0 * y;
  UT_ASSERT_EQUALS(z[0], 12);
  UT_ASSERT_EQUALS(z[1], 15);
  UT_ASSERT_EQUALS(z[2], 18);

  z = y * 3.0;
  UT_ASSERT_EQUALS(z[0], 12);
  UT_ASSERT_EQUALS(z[1], 15);
  UT_ASSERT_EQUALS(z[2], 18);

  x = x - y;
  UT_ASSERT_EQUALS(x[0], 1);
  UT_ASSERT_EQUALS(x[1], 2);
  UT_ASSERT_EQUALS(x[2], 3);

  x.set(y);
  UT_ASSERT_EQUALS(x[0], y[0]);
  UT_ASSERT_EQUALS(x[1], y[1]);
  UT_ASSERT_EQUALS(x[2], y[2]);
}
UT_TEST_CASE_END(vecd_test)

UT_TEST_CASE(vecs_test) {
  double tol = 1e-12;

  // testing double assignment (note the . so that we explicitly use double)
  vec3d x = {1., 2., 3.};
  UT_ASSERT_NEAR(dot(x, x), 14., tol);

  // testing assignment of an integer vector to a double vector
  vecs<3, int> yi = {1, 2, 3};
  vec3d yd = yi;
  vec3d yd2 = {1, 2, 3};
  for (int i = 0; i < 3; i++) {
    UT_ASSERT_EQUALS(yd[i], i + 1);
    UT_ASSERT_EQUALS(yd2[i], i + 1);
  }

  vec3d z = x / 2.0;
  UT_ASSERT_NEAR(z[0], 0.5, tol);
  UT_ASSERT_NEAR(z[1], 1.0, tol);
  UT_ASSERT_NEAR(z[2], 1.5, tol);

  const double& z0 = z[0];
  UT_ASSERT_NEAR(z0, 0.5, tol);
}
UT_TEST_CASE_END(vecs_test)

UT_TEST_SUITE_END(vector_suite)
