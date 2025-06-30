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

#include "util.h"

#include <limits>
#include <vector>

#include "io.h"
#include "library.h"
#include "log.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(util_test_suite)

UT_TEST_CASE(test_get_bounding_box) {
  const int dim = 2;
  const int64_t n_points = 7;
  std::vector<double> points = {1.0, 5.0, 3.0,  7.0, -2.0, 4.0, 0.5, -0.5,
                                0.0, 0.0, -0.6, 3.0, 0.0,  0.0

  };

  double xmin[dim], xmax[dim];

  get_bounding_box(points.data(), n_points, dim, xmin, xmax);

  UT_ASSERT_EQUALS(xmin[0], -2.0);
  UT_ASSERT_EQUALS(xmin[1], -0.5);

  UT_ASSERT_EQUALS(xmax[0], 3.0);
  UT_ASSERT_EQUALS(xmax[1], 7.0);
}

UT_TEST_CASE_END(test_get_bounding_box)

UT_TEST_CASE(test_sort_points_on_zcurve_one_z) {
  const int dim = 2;
  const int64_t n_points = 4;
  std::vector<double> points = {

      0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0,

  };

  std::vector<index_t> order(n_points);
  sort_points_on_zcurve(points.data(), n_points, dim, order);

  Mesh mesh(dim);
  coord_t x[dim];

  for (int i = 0; i < n_points; i++) {
    int idx = order[i];
    x[0] = points[dim * idx + 0];  // x
    x[1] = points[dim * idx + 1];  // y
    mesh.vertices().add(x);
  }

  for (int i = 1; i < n_points; i++) {
    index_t line[2] = {static_cast<index_t>(i - 1), static_cast<index_t>(i)};
    mesh.lines().add(line, 2);
  }

  meshb::write(mesh, "zcurve.meshb");

  std::vector<index_t> expected_order = {0, 1, 2, 3};
  UT_ASSERT(order.size() == expected_order.size());
  for (size_t i = 0; i < order.size(); ++i) {
    UT_ASSERT(order[i] == expected_order[i]);
  }

  // Prints out the actual order of the z-curve
  LOG << "Z-curve order: ";
  for (auto idx : order) {
    std::cout << idx << " ";
  }
  std::cout << std::endl;
}

UT_TEST_CASE_END(test_sort_points_on_zcurve_one_z)

UT_TEST_CASE(test_sort_points_on_z_curve_grid) {
  const int dim = 3;
  const int nx = 4;
  const int ny = 4;

  Grid<Quad> grid({nx, ny});

  const int64_t n_points = grid.vertices().n();
  std::vector<double> points(n_points * dim);

  for (int64_t i = 0; i < n_points; ++i) {
    const double* v = grid.vertices()[i];
    points[dim * i + 0] = v[0];
    points[dim * i + 1] = v[1];
  }

  std::vector<index_t> order(n_points);
  sort_points_on_zcurve(points.data(), n_points, dim, order);

  Mesh mesh(dim);
  coord_t x[dim];
  for (int64_t i = 0; i < n_points; ++i) {
    int idx = order[i];
    x[0] = points[dim * idx + 0];  // x
    x[1] = points[dim * idx + 1];  // y
    mesh.vertices().add(x);
  }

  for (int i = 1; i < n_points; ++i) {
    index_t line[2] = {static_cast<index_t>(i - 1), static_cast<index_t>(i)};
    mesh.lines().add(line, 2);
  }

  meshb::write(mesh, "zcurve_grid.meshb");
}
UT_TEST_CASE_END(test_sort_points_on_z_curve_grid)

UT_TEST_SUITE_END(util_test_suite)