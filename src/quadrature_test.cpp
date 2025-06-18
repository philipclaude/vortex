//
//  vortex: Voronoi mesher and fluid simulator for the Earth's oceans and
//  atmosphere.
//
//  Copyright 2023 - 2025 Philip Claude Caplan
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
#include "quadrature.hpp"

#include "library.h"
#include "math/vec.hpp"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(quadrature_test_suite)

UT_TEST_CASE(planar_triangle_tests) {
  TriangleQuadrature<Triangle> quad(5);

  Grid<Triangle> mesh({4, 4});
  LOG << fmt::format("# triangles = {}", mesh.triangles().n());

  auto check_integration = [&quad, &mesh](const auto& fn, double value,
                                          double tol) {
    double integral = 0;
    for (size_t k = 0; k < mesh.triangles().n(); k++) {
      const auto* t = mesh.triangles()[k];
      const auto* pa = mesh.vertices()[t[0]];
      const auto* pb = mesh.vertices()[t[1]];
      const auto* pc = mesh.vertices()[t[2]];
      integral += quad.integrate(fn, pa, pb, pc);
    }
    LOG << fmt::format("integral = {}, analytic = {}, error = {:1.16e}",
                       integral, value, std::fabs(integral - value));
    return std::fabs(integral - value) < tol;
  };

  auto fn0 = [](const vec3d& x) -> double { return 1.0; };
  UT_ASSERT(check_integration(fn0, 1.0, 1e-15));

  vec3d z = {0.5, 0.5, 0.0};
  auto fn1 = [&z](const vec3d& x) -> double {
    return std::pow(length(x - z), 2);
  };
  UT_ASSERT(check_integration(fn1, 1.0 / 6.0, 1e-15));

  auto fn2 = [](const vec3d& x) -> double {
    return std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]);
  };
  UT_ASSERT(check_integration(fn2, 4.0 / (M_PI * M_PI), 1e-6));
}
UT_TEST_CASE_END(planar_triangle_tests)

UT_TEST_CASE(spherical_triangle_tests) {
  TriangleQuadrature<SphericalTriangle> quad(4);

  SubdividedIcosahedron mesh(1);
  LOG << fmt::format("# triangles = {}", mesh.triangles().n());

  auto check_integration = [&quad, &mesh](const auto& fn, double value) {
    double integral = 0;
    for (size_t k = 0; k < mesh.triangles().n(); k++) {
      const auto* t = mesh.triangles()[k];
      const auto* pa = mesh.vertices()[t[0]];
      const auto* pb = mesh.vertices()[t[1]];
      const auto* pc = mesh.vertices()[t[2]];
      integral += quad.integrate(fn, pa, pb, pc);
    }
    LOG << fmt::format("integral = {}, analytic = {}, error = {:1.16e}",
                       integral, value, std::fabs(integral - value));
    return std::fabs(integral - value) < 1e-12;
  };

  auto fn0 = [](const vec3d& x) -> double { return 1.0; };
  UT_ASSERT(check_integration(fn0, 4 * M_PI));

  vec3d z = {1, 0, 0};
  auto fn1 = [&z](const vec3d& x) -> double {
    return std::pow(length(x - z), 2);
  };
  UT_ASSERT(check_integration(fn1, 8 * M_PI));

  auto fn2 = [](const vec3d& x) -> double { return x[0] * x[2] * x[2]; };
  UT_ASSERT(check_integration(fn2, 0.0));

  auto fn3 = [](const vec3d& x) -> double { return x[0] + x[2] * x[2]; };
  UT_ASSERT(check_integration(fn3, 4.0 * M_PI / 3.0));
}
UT_TEST_CASE_END(spherical_triangle_tests)

UT_TEST_SUITE_END(quadrature_test_suite)