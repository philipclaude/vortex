#include "vec.h"

#include "linalg.h"
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

  // testing constructors
  // default constructor
  vecs<3, double> v1;
  // constructor from an array of values
  double a[3];
  a[0] = .05;
  a[1] = 99;
  a[2] = -3.78;
  vecs<3, double> v2(a);
  // constructor setting all entries to an int. value
  vecs<3, int> v3(3);
  // constructor copying components of one vector to another (same type)
  vecs<3, double> v4(v2);
  // constructor copying components of one vector to another (different type)
  vecs<3, double> v5(v3);
  // constructor copying components of a M x 1 matrix
  mats<3, 1, double> m;
  vecs<3, double> v6(m);
  // using the default constructor and "=" operator to create a vector by
  // assigning equality to another vector
  v1 = v2;
  // constructors using initializer lists
  vecs<3, double> v7 = {1.3, 2.6645, -31.903};
  vecs<3, double> v8 = {5, 9, 22};
  std::initializer_list<double> l = {3.5, 6.2, 9.};
  vecs<3, double> v9(l);
  // constructor setting all entries of a vector equal to a scalar value
  vecs<3, int> v10 = 13;

  for (int i = 0; i < 3; i++) {
    UT_ASSERT_EQUALS(a[i], v2[i]);
    UT_ASSERT_EQUALS(3, v3[i]);
    UT_ASSERT_EQUALS(v2[i], v4[i]);
    UT_ASSERT_EQUALS(v3[i], v5[i]);
    UT_ASSERT_EQUALS(m(i, 0), v6[i]);
    UT_ASSERT_EQUALS(v2[i], v1[i]);
    UT_ASSERT_EQUALS(v10[i], 13);
  }
  // testing for vectors v7 and v8
  UT_ASSERT_EQUALS(v7[0], 1.3);
  UT_ASSERT_EQUALS(v7[1], 2.6645);
  UT_ASSERT_EQUALS(v7[2], -31.903);
  UT_ASSERT_EQUALS(v8[0], 5);
  UT_ASSERT_EQUALS(v8[1], 9);
  UT_ASSERT_EQUALS(v8[2], 22);
  UT_ASSERT_EQUALS(v9[0], 3.5);
  UT_ASSERT_EQUALS(v9[1], 6.2);
  UT_ASSERT_EQUALS(v9[2], 9.);

  // testing length
  // 3-dimensional vectors with int. entries (for testing's sake)
  vecs<3, int> vi1 = {99, 34, 76};
  vecs<3, int> vi2 = {3, 4, 5};
  vecs<3, int> vi3 = {17, 2, 22};

  double l1 = length(vi1);
  double l2 = length(vi2);
  double l3 = length(vi3);

  UT_ASSERT_NEAR(l1, floor(sqrt(16733)), tol);
  UT_ASSERT_NEAR(l2, floor(sqrt(50)), tol);
  UT_ASSERT_NEAR(l3, floor(sqrt(777)), tol);

  // 3-dimensional vectors with double entries
  vecs<3, double> vd1 = {2.789, 0.384, 59.2};
  vecs<3, double> vd2 = {4.8, 3.7, 29.00064};
  vecs<3, double> vd3 = {3.55, -26.09, 54.54};
  // 3d vectors with z-value 0
  vecs<3, double> vd4 = {65.33, -12.97, 0.};
  vecs<3, double> vd5 = {34.92, 15.09917, 0.};
  // 4d vectors
  vecs<4, double> vd6 = {-793.672, 29.292929, 0., -0.0023};
  vecs<4, double> vd7 = {-3.922, -27., -89.98, -55.2};
  vecs<4, double> vd8 = {5., 5., 5., 5.};

  double l4 = length(vd1);
  double l5 = length(vd2);
  double l6 = length(vd3);
  double l7 = length(vd4);
  double l8 = length(vd5);
  double l9 = length(vd6);
  double l10 = length(vd7);
  double l11 = length(vd8);

  UT_ASSERT_NEAR(l4, (sqrt(3512565977) / 1000), tol);
  UT_ASSERT_NEAR(l5, (sqrt(34287778141) / 6250), tol);
  UT_ASSERT_NEAR(l6, (sqrt(36679022) / 100), tol);
  UT_ASSERT_NEAR(l7, (sqrt(44362298) / 100), tol);
  UT_ASSERT_NEAR(l8, (sqrt(14473913346889) / 100000), tol);
  UT_ASSERT_NEAR(l9, (sqrt(630773319278689041) / 1000000), tol);
  UT_ASSERT_NEAR(l10, (sqrt(2971955621) / 500), tol);
  UT_ASSERT_NEAR(l11, 10, tol);

  // testing normalization of 3 & 4-dimensional vectors with double entries
  auto ul1 = normalize(vd1);
  auto ul2 = normalize(vd2);
  auto ul3 = normalize(vd3);
  auto ul4 = normalize(vd4);
  auto ul5 = normalize(vd5);
  auto ul6 = normalize(vd6);

  UT_ASSERT_NEAR(1., length(ul1), tol);
  UT_ASSERT_NEAR(1., length(ul2), tol);
  UT_ASSERT_NEAR(1., length(ul3), tol);
  UT_ASSERT_NEAR(1., length(ul4), tol);
  UT_ASSERT_NEAR(1., length(ul5), tol);
  UT_ASSERT_NEAR(1., length(ul6), tol);

  // testing cross product for 3-dimensional vectors with double entries
  vecs<3, double> cp1 = cross(vd1, vd2);
  vecs<3, double> cp2 = cross(vd2, vd1);
  vecs<3, double> cp3 = cross(vd2, vd3);
  vecs<3, double> cp4 = cross(vd3, vd2);
  vecs<3, double> cp5 = cross(vd1, vd3);
  vecs<3, double> cp6 = cross(vd3, vd1);
  vecs<3, double> cp7 = cross(vd4, vd5);
  vecs<3, double> cp8 = cross(vd5, vd4);
  // corresponding cross-product solutions to check
  vecs<3, double> a1 = {-207.90375424, 203.27721504, 8.4761};
  vecs<3, double> a2 = {207.90375424, -203.27721504, -8.4761};
  vecs<3, double> a3 = {958.4246976, -158.839728, -138.367};
  vecs<3, double> a4 = {-958.4246976, 158.839728, 138.367};
  vecs<3, double> a5 = {1565.47136, 58.04794, -74.12821};
  vecs<3, double> a6 = {-1565.47136, -58.04794, 74.12821};
  vecs<3, double> a7 = {0., 0., 1439.3411761};
  vecs<3, double> a8 = {0., 0., -1439.3411761};

  for (int i = 0; i < 3; i++) {
    UT_ASSERT_NEAR(cp1[i], a1[i], tol);
    UT_ASSERT_NEAR(cp2[i], a2[i], tol);
    UT_ASSERT_NEAR(cp3[i], a3[i], tol);
    UT_ASSERT_NEAR(cp4[i], a4[i], tol);
    UT_ASSERT_NEAR(cp5[i], a5[i], tol);
    UT_ASSERT_NEAR(cp6[i], a6[i], tol);
    UT_ASSERT_NEAR(cp7[i], a7[i], tol);
    UT_ASSERT_NEAR(cp8[i], a8[i], tol);
  }

  // testing vec-vec (+, -, (-), *) and vec-sca (*, /) operators
  // vec-vec operations
  vecs<3, double> s1 = vd1 + vd2;
  vecs<3, double> s2 = vd2 + vd3;
  vecs<3, double> s3 = vd3 + vd1;
  vecs<3, double> d1 = vd1 - vd2;
  vecs<3, double> d2 = vd2 - vd3;
  vecs<3, double> d3 = vd3 - vd1;
  vecs<3, double> n1 = -vd1;
  vecs<3, double> n2 = -vd2;
  vecs<3, double> n3 = -vd3;
  vecs<3, double> m1 = vd1 * vd2;
  vecs<3, double> m2 = vd2 * vd3;
  vecs<3, double> m3 = vd3 * vd1;
  vecs<3, double> m4 = vd5 * vd4;

  // vec-sca operations
  vecs<3, double> sc1 = 3.5 * vd1;
  vecs<3, double> sc2 = vd2 * 57.0158;
  vecs<3, double> sc3 = -19.1 * vd3;
  vecs<3, double> div1 = vd1 / 8.;
  vecs<3, double> div2 = vd2 / -6.23;
  vecs<3, double> div3 = vd3 / 0.5;

  for (int i = 0; i < 3; i++) {
    UT_ASSERT_NEAR(s1[i], vd1[i] + vd2[i], tol);
    UT_ASSERT_NEAR(s2[i], vd2[i] + vd3[i], tol);
    UT_ASSERT_NEAR(s3[i], vd3[i] + vd1[i], tol);
    UT_ASSERT_NEAR(d1[i], vd1[i] - vd2[i], tol);
    UT_ASSERT_NEAR(d2[i], vd2[i] - vd3[i], tol);
    UT_ASSERT_NEAR(d3[i], vd3[i] - vd1[i], tol);
    UT_ASSERT_NEAR(n1[i], -vd1[i], tol);
    UT_ASSERT_NEAR(n2[i], -vd2[i], tol);
    UT_ASSERT_NEAR(n3[i], -vd3[i], tol);
    UT_ASSERT_NEAR(m1[i], vd1[i] * vd2[i], tol);
    UT_ASSERT_NEAR(m2[i], vd2[i] * vd3[i], tol);
    UT_ASSERT_NEAR(m3[i], vd3[i] * vd1[i], tol);
    UT_ASSERT_NEAR(m4[i], vd5[i] * vd4[i], tol);
    UT_ASSERT_NEAR(sc1[i], 3.5 * vd1[i], tol);
    UT_ASSERT_NEAR(sc2[i], 57.0158 * vd2[i], tol);
    UT_ASSERT_NEAR(sc3[i], -19.1 * vd3[i], tol);
    UT_ASSERT_NEAR(div1[i], vd1[i] / 8., tol);
    UT_ASSERT_NEAR(div2[i], vd2[i] / -6.23, tol);
    UT_ASSERT_NEAR(div3[i], vd3[i] / 0.5, tol);
  }

  // Testing outter product
  vecs<3, double> uo1 = {1., 2., 3.};
  vecs<3, double> vo1 = {4., 5., 6.};
  mats<3, 3, double> o1 = outer(uo1, vo1);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      UT_ASSERT_EQUALS(o1(i, j), uo1(i) * vo1(j));
    }
  }

  vecs<3, double> uo2 = {23.556, 2.79, 0.9332};
  vecs<6, double> vo2 = {-4.6, -5.442, 6.804, 13.21, 16.781, 100.100};
  mats<3, 6, double> o2 = outer(uo2, vo2);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 6; j++) {
      UT_ASSERT_EQUALS(o2(i, j), uo2(i) * vo2(j));
    }
  }

  vecs<4, double> uo3 = {1.21, 4.54, 7.87, 89.89};
  vecs<2, double> vo3 = {-10000.6, -312.1688883};
  mats<4, 2, double> o3 = outer(uo3, vo3);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 2; j++) {
      UT_ASSERT_EQUALS(o3(i, j), uo3(i) * vo3(j));
    }
  }
  // Manual tests
  vecs<2, double> uo4 = {2., 2.};
  vecs<2, double> vo4 = {5., 5.};
  mats<2, 2, double> o4 = outer(uo4, vo4);
  UT_ASSERT_EQUALS(o1(0, 0), 4.);
  UT_ASSERT_EQUALS(o1(0, 1), 5.);
  UT_ASSERT_EQUALS(o1(0, 2), 6.);
  UT_ASSERT_EQUALS(o1(1, 0), 8.);
  UT_ASSERT_EQUALS(o1(1, 1), 10.);
  UT_ASSERT_EQUALS(o1(1, 2), 12.);
  UT_ASSERT_EQUALS(o1(2, 0), 12.);
  UT_ASSERT_EQUALS(o1(2, 1), 15.);
  UT_ASSERT_EQUALS(o1(2, 2), 18.);
  UT_ASSERT_EQUALS(o4(0, 0), 10.);
  UT_ASSERT_EQUALS(o4(0, 1), 10.);
  UT_ASSERT_EQUALS(o4(1, 0), 10.);
  UT_ASSERT_EQUALS(o4(1, 1), 10.);
}
UT_TEST_CASE_END(vecs_test)

UT_TEST_SUITE_END(vector_suite)
