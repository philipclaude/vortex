//
//  vortex: Voronoi mesher and fluid simulator for the Earth's oceans and
//  atmosphere.
//
//  Copyright 2023 - 2024 Philip Claude Caplan
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
#include "library.h"
#include "math/vec.hpp"
#include "quadrature.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(elements_test_suite)

UT_TEST_CASE(spherical_triangle_jacobian_tests) {
  SubdividedIcosahedron mesh(0);
  TriangleQuadrature<SphericalTriangle> quad(10);

  double area_approx = 0.0;
  double area_analytic = 0.0;
  for (size_t elem = 0; elem < mesh.triangles().n(); elem++) {
    const auto* t = mesh.triangles()[elem];
    const auto* a = mesh.vertices()[t[0]];
    const auto* b = mesh.vertices()[t[1]];
    const auto* c = mesh.vertices()[t[2]];
    for (size_t k = 0; k < quad.points().size(); k++) {
      auto& q = quad.points()[k];
      double da = SphericalTriangle::jacobian(a, b, c, &q[0]);
      double w = quad.weights()[k];
      area_approx += da * w;
    }
    area_analytic += SphericalTriangle::area(a, b, c);
  }
  double area_exact = 4 * M_PI;
  double error_approx = std::fabs(area_exact - area_approx);
  LOG << fmt::format("area_approx = {}, area_analytic = {}, error = {}",
                     area_approx, area_analytic, error_approx);
  UT_ASSERT(error_approx < 1e-4);
}
UT_TEST_CASE_END(spherical_triangle_jacobian_tests)

UT_TEST_SUITE_END(elements_test_suite)