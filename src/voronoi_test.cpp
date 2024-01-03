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
#include "voronoi.h"

#include <Predicates_psm.h>

#include <cmath>

#include "graphics.h"
#include "io.h"
#include "library.h"
#include "math/vec.hpp"
#include "tester.h"
#include "util.h"

using namespace vortex;

UT_TEST_SUITE(voronoi_test_suite)

UT_TEST_CASE(test_pool) {
  // basic tests
  {
    std::vector<int> x_pool(300);

    pool<int> x(x_pool.data(), x_pool.capacity());
    UT_ASSERT_EQUALS(x.size(), 0);
    UT_ASSERT(x.empty());
    x.push_back(1);
    UT_ASSERT_EQUALS(x.capacity(), 300);
    UT_ASSERT_EQUALS(x.size(), 1);

    for (int i = 0; i < 200; i++) {
      x.push_back(i + 2);
    }
    UT_ASSERT_EQUALS(x.size(), 201);
    UT_ASSERT_EQUALS(x.capacity(), 300);

    for (size_t i = 0; i < x.size(); i++) {
      UT_ASSERT_EQUALS(size_t(x[i]), i + 1);
    }

    UT_ASSERT_EQUALS(x.back(), 201);

    x.clear();
    UT_ASSERT(x.empty());
    UT_ASSERT_EQUALS(x.capacity(), 300);
  }

  // swap test
  {
    std::vector<int> x_pool(3), y_pool(4);
    pool<int> x(x_pool.data(), 4);
    pool<int> y(y_pool.data(), 4);
    x.set_size(3);
    y.set_size(4);
    x[0] = 1;
    x[1] = 4;
    x[2] = 9;
    y[0] = 2;
    y[1] = 10;
    y[2] = 0;
    y[3] = 5;
    x.swap(y);

    UT_ASSERT_EQUALS(x.size(), 4);
    UT_ASSERT_EQUALS(y.size(), 3);
    UT_ASSERT_EQUALS(x[0], 2);
    UT_ASSERT_EQUALS(x[1], 10);
    UT_ASSERT_EQUALS(x[2], 0);
    UT_ASSERT_EQUALS(x[3], 5);
    UT_ASSERT_EQUALS(y[0], 1);
    UT_ASSERT_EQUALS(y[1], 4);
    UT_ASSERT_EQUALS(y[2], 9);
  }
}
UT_TEST_CASE_END(test_pool)

UT_TEST_CASE(test_square) {
  const double tol = 1e-12;
  static const int dim = 4;
  size_t n_sites = 1e3;
  std::vector<coord_t> sites(n_sites * dim, 0.0);
  for (size_t k = 0; k < n_sites; k++) {
    sites[k * dim + 0] = double(rand()) / double(RAND_MAX);
    sites[k * dim + 1] = double(rand()) / double(RAND_MAX);
    sites[k * dim + 2] = 0.0;
    if (dim > 3) sites[k * dim + 3] = 0.0;
  }

  std::vector<index_t> order(n_sites);
  sort_points_on_zcurve(sites.data(), n_sites, dim, order);

  Vertices vertices(dim);
  vertices.reserve(n_sites);
  coord_t x[dim];
  for (size_t i = 0; i < n_sites; i++) {
    for (int d = 0; d < dim; d++) x[d] = sites[dim * order[i] + d];
    vertices.add(x);
  }

  SquareDomain domain;
  VoronoiDiagram voronoi(dim, vertices[0], n_sites);
  VoronoiDiagramOptions options;
  options.n_neighbors = 75;
  options.allow_reattempt = false;
  options.parallel = true;
  int n_iter = 20;
  auto& weights = voronoi.weights();
  weights.resize(n_sites, 0.0);
  for (int iter = 1; iter <= n_iter; ++iter) {
    options.store_mesh = iter == n_iter;
    options.verbose = (iter == 1 || iter == n_iter - 1);
    voronoi.vertices().clear();
    voronoi.vertices().set_dim(3);
    voronoi.polygons().clear();
    voronoi.triangles().clear();
    voronoi.compute(domain, options);

    // move each site to the centroid of the corresponding cell
    voronoi.smooth(vertices);
    auto props = voronoi.analyze();
    LOG << fmt::format("iter = {}, area = {}", iter, props.area);
  }
  for (size_t k = 0; k < n_sites; k++) {
    // coordinates relative to square center (0.5, 0.5)
    double x = vertices[k][0] - 0.5;
    double y = vertices[k][1] - 0.5;
    double r = std::sqrt(x * x + y * y);
    weights[k] = 2e-1 * std::pow(r - 0.5, 2);
    if (r > 0.5) weights[k] = 0;
    // weights[k] = 2e-1 * std::exp(0.5 - r);
  }
  lift_sites(vertices, voronoi.weights());

  voronoi.vertices().clear();
  voronoi.vertices().set_dim(3);
  voronoi.polygons().clear();
  voronoi.triangles().clear();
  options.store_mesh = true;
  options.verbose = true;
  voronoi.compute(domain, options);
  auto props = voronoi.analyze();
  LOG << fmt::format("power diagram area = {}", props.area);
  UT_ASSERT_NEAR(props.area, 1.0, tol);

  // check the power diagram
  UT_ASSERT_EQUALS(voronoi.vertices().n() - n_sites, voronoi.triangles().n());
  for (size_t k = 0; k < voronoi.triangles().n(); k++) {
    auto* t = voronoi.triangles()[k];

    // only check interior triangles
    if (voronoi.triangles().group(k) < 0) continue;

    // coordinates of voronoi vertex
    auto* v = voronoi.vertices()[k + n_sites];
    vec3d c(v);

    // coordinates of 3 vertices of delaunay triangle
    vec3d zi(vertices[t[0]]);
    vec3d zj(vertices[t[1]]);
    vec3d zk(vertices[t[2]]);
    ASSERT(t[0] != t[1]) << t[0] << ", " << t[1];
    ASSERT(t[1] != t[2]) << t[1] << ", " << t[2];
    ASSERT(t[0] != t[2]) << t[0] << ", " << t[2];

    vec3d n = normalize(cross(zj - zi, zk - zi));
    UT_ASSERT_NEAR(n[2], 1.0, tol);

    // power distance to sites
    double di = std::pow(length(zi - c), 2) - weights[t[0]];
    double dj = std::pow(length(zj - c), 2) - weights[t[1]];
    double dk = std::pow(length(zk - c), 2) - weights[t[2]];
    UT_ASSERT_NEAR(di, dj, tol);
    UT_ASSERT_NEAR(di, dk, tol);
    UT_ASSERT_NEAR(dj, dk, tol);
  }

  voronoi.merge();
  LOG << fmt::format("writing {} polygons", voronoi.polygons().n());
  if (voronoi.polygons().n() > 0) meshb::write(voronoi, "square.meshb");
}
UT_TEST_CASE_END(test_square)

UT_TEST_CASE(test_sphere) {
  auto irand = [](int min, int max) {
    return min + double(rand()) / (double(RAND_MAX) + 1.0) * (max - min);
  };
  const double tol = 1e-12;
  static const int dim = 4;
  size_t n_sites = 1e4;
  std::vector<coord_t> sites(n_sites * dim, 0.0);
  for (size_t k = 0; k < n_sites; k++) {
    coord_t theta = 2.0 * M_PI * irand(0, 1);
    coord_t phi = acos(2.0 * irand(0, 1) - 1.0);
    sites[k * dim + 0] = cos(theta) * sin(phi);
    sites[k * dim + 1] = sin(theta) * sin(phi);
    sites[k * dim + 2] = cos(phi);
    if (dim > 3) sites[k * dim + 3] = 0.0;
  }

  std::vector<index_t> order(n_sites);
  sort_points_on_zcurve(sites.data(), n_sites, dim, order);

  Vertices vertices(dim);
  vertices.reserve(n_sites);
  coord_t x[dim];
  for (size_t i = 0; i < n_sites; i++) {
    for (int d = 0; d < dim; d++) x[d] = sites[dim * order[i] + d];
    vertices.add(x);
  }

  SphereDomain domain(1.0);
  VoronoiDiagram voronoi(dim, vertices[0], n_sites);
  VoronoiDiagramOptions options;
  options.n_neighbors = 75;
  options.allow_reattempt = false;
  options.parallel = true;
  int n_iter = 20;
  auto& weights = voronoi.weights();
  weights.resize(n_sites, 0.0);
  for (int iter = 1; iter <= n_iter; ++iter) {
    options.store_mesh = iter == n_iter;
    options.verbose = (iter == 1 || iter == n_iter - 1);
    voronoi.vertices().clear();
    // voronoi.vertices().set_dim(3);
    voronoi.polygons().clear();
    voronoi.triangles().clear();
    voronoi.compute(domain, options);

    // move each site to the centroid of the corresponding cell
    voronoi.smooth(vertices);
    auto props = voronoi.analyze();
    LOG << fmt::format("iter = {}, area = {}", iter, props.area);
  }

  for (size_t k = 0; k < n_sites; k++) {
    double x = vertices[k][0];
    double y = vertices[k][1];
    double r = std::sqrt(x * x + y * y);
    // weights[k] = 4e-1 * std::pow(r - 0.5, 2);
    //  if (r > 0.5) weights[k] = 0;
    weights[k] = 2e-1 * std::exp(0.5 - r);
  }
  lift_sites(vertices, voronoi.weights());

  voronoi.vertices().clear();
  // voronoi.vertices().set_dim(3);
  voronoi.polygons().clear();
  voronoi.triangles().clear();
  voronoi.compute(domain, options);
  auto props = voronoi.analyze();
  LOG << fmt::format("power diagram area = {}", props.area);
  UT_ASSERT_NEAR(props.area, 4 * M_PI, tol);

  // check the power diagram
  UT_ASSERT_EQUALS(voronoi.vertices().n() - n_sites, voronoi.triangles().n());
  for (size_t k = 0; k < voronoi.triangles().n(); k++) {
    auto* t = voronoi.triangles()[k];

    // only check interior triangles
    if (voronoi.triangles().group(k) < 0) continue;

    // coordinates of voronoi vertex
    auto* v = voronoi.vertices()[k + n_sites];
    vec3d c(v);

    // coordinates of 3 vertices of delaunay triangle
    vec3d zi(vertices[t[0]]);
    vec3d zj(vertices[t[1]]);
    vec3d zk(vertices[t[2]]);
    ASSERT(t[0] != t[1]) << t[0] << ", " << t[1];
    ASSERT(t[1] != t[2]) << t[1] << ", " << t[2];
    ASSERT(t[0] != t[2]) << t[0] << ", " << t[2];

    // power distance to sites
    double di = std::pow(length(zi - c), 2) - weights[t[0]];
    double dj = std::pow(length(zj - c), 2) - weights[t[1]];
    double dk = std::pow(length(zk - c), 2) - weights[t[2]];
    UT_ASSERT_NEAR(di, dj, tol);
    UT_ASSERT_NEAR(di, dk, tol);
    UT_ASSERT_NEAR(dj, dk, tol);
  }

  voronoi.merge();
  LOG << fmt::format("writing {} polygons", voronoi.polygons().n());
  if (voronoi.polygons().n() > 0) meshb::write(voronoi, "sphere.meshb");
}
UT_TEST_CASE_END(test_sphere)

UT_TEST_CASE(test_sphere_triangulation) {
  Sphere sphere(2);
  static const int dim = 3;
  size_t n_sites = 1e4;
  Vertices data(3);
  sample_surface(sphere, data, n_sites);
  auto& sites = data.data();

  std::vector<index_t> order(n_sites);
  sort_points_on_zcurve(sites.data(), n_sites, dim, order);
  Vertices vertices(dim);
  vertices.reserve(n_sites);
  coord_t x[dim];
  for (size_t i = 0; i < n_sites; i++) {
    for (int d = 0; d < dim; d++) x[d] = sites[dim * order[i] + d];
    vertices.add(x);
  }

  TriangulationDomain domain(sphere.vertices()[0], sphere.vertices().n(),
                             sphere.triangles()[0], sphere.triangles().n());
  VoronoiDiagram voronoi(dim, vertices[0], n_sites);
  VoronoiDiagramOptions options;
  options.n_neighbors = 75;
  options.allow_reattempt = false;
  options.parallel = true;
  int n_iter = 5;
  for (int iter = 1; iter <= n_iter; ++iter) {
    options.store_mesh = iter == n_iter;
    options.verbose = (iter == 1 || iter == n_iter - 1);
    voronoi.vertices().clear();
    voronoi.vertices().set_dim(3);
    voronoi.triangles().clear();
    voronoi.polygons().clear();
    voronoi.compute(domain, options);

    // move each site to the centroid of the corresponding cell
    voronoi.smooth(vertices);
    auto props = voronoi.analyze();
    LOG << fmt::format("iter = {}, area = {}", iter, props.area);
  }

  // randomize the colors a bit, otherwise neighboring cells
  // will have similar colors and won't visually stand out
  size_t n_colors = 20;
  std::vector<int> site2color(n_sites);
  for (size_t k = 0; k < n_sites; k++)
    site2color[k] = int(n_colors * double(rand()) / double(RAND_MAX));
  for (size_t k = 0; k < voronoi.polygons().n(); k++) {
    int group = voronoi.polygons().group(k);  // the group is the site
    voronoi.polygons().set_group(k, site2color[group]);
  }

  voronoi.merge();

  // voronoi.triangles().clear();  // not implemented yet
  LOG << fmt::format("writing {} polygons", voronoi.polygons().n());
  if (voronoi.polygons().n() > 0)
    meshb::write(voronoi, "sphere_triangulation.meshb");
}
UT_TEST_CASE_END(test_sphere_triangulation)

UT_TEST_SUITE_END(voronoi_test_suite)