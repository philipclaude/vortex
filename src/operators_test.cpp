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
#include "operators.h"

#include <cmath>

#include "graphics.h"
#include "io.h"
#include "library.h"
#include "math/mat.hpp"
#include "math/vec.hpp"
#include "tester.h"
#include "util.h"
#include "voronoi.h"

using namespace vortex;

UT_TEST_SUITE(voronoi_operators_test_suite)

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

  typedef SquareDomain Domain_t;
  Domain_t domain;
  VoronoiDiagram voronoi(dim, vertices[0], n_sites);
  VoronoiDiagramOptions options;
  options.n_neighbors = 75;
  options.parallel = true;
  // options.neighbor_algorithm = NearestNeighborAlgorithm::kKdtree;
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
    voronoi.smooth(vertices, false);
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
  UT_ASSERT_NEAR(props.area, domain.area(), tol);

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

  std::vector<double> f(n_sites, 0.0), grad_f(n_sites * 3, 0.0);

  for (int i = 0; i < n_sites; i++) {
    vec3d p(vertices[i]);
    f[i] = 0.5 * p[0] + 1.25 * p[1];
  }

  VoronoiOperators<Domain_t> ops(voronoi);
  ops.calculate_gradient(f.data(), grad_f.data());

  for (int i = 0; i < n_sites; i++) {
    // skip boundary polygons since the gradient is not accurate for these
    // (missing boundary term)
    vec3d gi(&grad_f[3 * i]);
    if (gi[0] >= 1e20) continue;
    UT_ASSERT_NEAR(gi[0], 0.5, 1e-12);
    UT_ASSERT_NEAR(gi[1], 1.25, 1e-12);
  }

  std::vector<double> u(n_sites * 3, 0.0), div_u(n_sites, 0.0);
  for (int i = 0; i < n_sites; i++) {
    vec3d p(vertices[i]);
    u[3 * i] = p[0] * 0.5;
    u[3 * i + 1] = p[1] * 1.25;
    u[3 * i + 2] = 0.0;
  }
  ops.calculate_divergence(u.data(), div_u.data());
  for (int i = 0; i < n_sites; i++) {
    if (div_u[i] >= 1e20) continue;  // skip boundary polygons
    UT_ASSERT_NEAR(div_u[i], 1.75, 1e-12);
  }

  LOG << fmt::format("writing {} polygons", voronoi.polygons().n());
  if (voronoi.polygons().n() > 0) meshb::write(voronoi, "square.meshb");
}
UT_TEST_CASE_END(test_square)

UT_TEST_CASE(test_plane) {
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

  vec3 p_u = {1, 0, 0};
  vec3 p_v = {0, 1, 0.5};
  vec3 normal = unit_vector(cross(p_u, p_v));
  LOG << fmt::format("normal = {}, {}, {}", normal[0], normal[1], normal[2]);
  typedef SquareDomain Domain_t;
  Domain_t domain({0, 0, 0}, p_u, p_v);
  VoronoiDiagram voronoi(dim, vertices[0], n_sites);
  VoronoiDiagramOptions options;
  options.n_neighbors = 75;
  options.parallel = true;
  // options.neighbor_algorithm = NearestNeighborAlgorithm::kKdtree;
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
    voronoi.smooth(vertices, false);
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
  UT_ASSERT_NEAR(props.area, domain.area(), tol);

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
    UT_ASSERT_NEAR(n[2], normal[2], tol);

    // power distance to sites
    double di = std::pow(length(zi - c), 2) - weights[t[0]];
    double dj = std::pow(length(zj - c), 2) - weights[t[1]];
    double dk = std::pow(length(zk - c), 2) - weights[t[2]];
    UT_ASSERT_NEAR(di, dj, tol);
    UT_ASSERT_NEAR(di, dk, tol);
    UT_ASSERT_NEAR(dj, dk, tol);
  }

  std::vector<double> f(n_sites, 0.0), grad_f(n_sites * 3, 0.0);

  for (int i = 0; i < n_sites; i++) {
    vec3d p(vertices[i]);
    f[i] = 1.5 * p[0] + 1.25 * p[1] + 2.75 * p[2];
  }

  VoronoiOperators<Domain_t> ops(voronoi);
  ops.calculate_gradient(f.data(), grad_f.data());

  mats<3, 3, double> projection;
  projection.eye();
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      projection(i, j) -= normal[i] * normal[j];
    }
  }
  vec3d ga = {1.5, 1.25, 2.75};
  ga = projection * ga;
  LOG << fmt::format("ga = {}, {}, {}", ga[0], ga[1], ga[2]);

  for (int i = 0; i < n_sites; i++) {
    // skip boundary polygons since the gradient is not accurate for these
    // (missing boundary term)
    vec3d gi(&grad_f[3 * i]);
    if (gi[0] >= 1e20) continue;

    UT_ASSERT_NEAR(gi[0], ga[0], 1e-10);
    UT_ASSERT_NEAR(gi[1], ga[1], 1e-10);
  }

  std::vector<double> u(n_sites * 3, 0.0), div_u(n_sites, 0.0);
  for (int i = 0; i < n_sites; i++) {
    vec3d p(vertices[i]);
    u[3 * i] = p[0] * ga[0];
    u[3 * i + 1] = p[1] * ga[1];
    u[3 * i + 2] = p[2] * ga[2];
  }
  projection.print();
  ops.calculate_divergence(u.data(), div_u.data());
  double div_a = ga[0] * (1 - normal[0] * normal[0]) +
                 ga[1] * (1 - normal[1] * normal[1]) +
                 ga[2] * (1 - normal[2] * normal[2]);
  LOG << fmt::format("div_a = {}", div_a);
  for (int i = 0; i < 10; i++) {
    if (div_u[i] >= 1e20) continue;  // skip boundary polygons
    UT_ASSERT_NEAR(div_u[i], div_a, 1e-12);
  }

  LOG << fmt::format("writing {} polygons", voronoi.polygons().n());
  if (voronoi.polygons().n() > 0) meshb::write(voronoi, "plane.meshb");
}
UT_TEST_CASE_END(test_plane)

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

  typedef SphereDomain Domain_t;
  Domain_t domain;
  domain.set_initialization_fraction(0.7);
  VoronoiDiagram voronoi(dim, vertices[0], n_sites);
  VoronoiDiagramOptions options;
  // options.neighbor_algorithm = NearestNeighborAlgorithm::kKdtree;
  options.n_neighbors = 75;
  options.parallel = true;
  options.store_facet_data = true;
  int n_iter = 2;
  auto& weights = voronoi.weights();
  weights.resize(n_sites, 0);
  for (int iter = 1; iter <= n_iter; ++iter) {
    options.store_mesh = iter == n_iter;
    options.verbose = (iter == 1 || iter == n_iter - 1);
    voronoi.vertices().clear();
    voronoi.polygons().clear();
    voronoi.triangles().clear();
    voronoi.compute(domain, options);

    // move each site to the centroid of the corresponding cell
    voronoi.smooth(vertices, true);
    auto props = voronoi.analyze();
    LOG << fmt::format("iter = {}, area = {}", iter, props.area);
  }

  for (size_t k = 0; k < weights.size(); k++) {
    // double x = vertices[k][0];
    // double y = vertices[k][1];
    // double r = std::sqrt(x * x + y * y);
    // weights[k] = 4e-1 * std::pow(r - 0.5, 2);
    //  if (r > 0.5) weights[k] = 0;
    weights[k] = 0.0;  // 2e-1 * std::exp(0.5 - r);
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

  std::vector<double> f(n_sites, 0.0), grad_f(n_sites * 3, 0.0);

  for (int i = 0; i < n_sites; i++) {
    vec3d p(vertices[i]);
    f[i] = 0.5 * p[0] + 1.25 * p[1] + 2.75 * p[2];

    // double phi = atan2(p[1], p[0]) + M_PI;
    // double theta = acos(p[2]);
    //  f[i] = 1.25 * phi;  // 12 * sin(3 * phi) * pow(sin(theta), 3.0);
  }

  VoronoiOperators<Domain_t> ops(voronoi);
  ops.calculate_gradient(f.data(), grad_f.data());

  double emin = 1e20;
  double emax = 0.0;
  for (int i = 0; i < n_sites; i++) {
    vec3d gi(&grad_f[3 * i]);
    UT_ASSERT(gi[0] < 1e20);  // should be closed

    vec3d p(vertices[i]);
    double phi = atan2(p[1], p[0]) + M_PI;
    double theta = acos(p[2]);
    vec3d gradient = {0.5, 1.25, 2.75};
    // vec3d gradient = {0.0, 36 * sin(3 * phi) * pow(sin(theta), 2),
    //                   36 * cos(3 * phi) * pow(sin(theta), 2)};
    // vec3d gradient = {0.0, 0.0, 1.25 / sin(theta)};
    mats<3, 3, double> m;
    m(0, 0) = sin(theta) * cos(phi);
    m(0, 1) = sin(theta) * sin(phi);
    m(0, 2) = cos(theta);
    m(1, 0) = cos(theta) * cos(phi);
    m(1, 1) = cos(theta) * sin(phi);
    m(1, 2) = -sin(theta);
    m(2, 0) = -sin(phi);
    m(2, 1) = cos(phi);
    m(2, 2) = 0.0;
    vec3d ga = gradient;
    Domain_t::project(vertices[i], &gradient[0], &ga[0]);
    double ei = length(ga - gi);
    if (ei < emin) emin = ei;
    if (ei > emax) emax = ei;
  }
  LOG << fmt::format("emin = {}, emax = {}", emin, emax);

  double a = 1;
  double b = 1;
  double c = 1;
  std::vector<double> u(n_sites * 3, 0.0), div_u(n_sites, 0.0);
  emin = 1e20;
  emax = 0.0;
  for (int i = 0; i < n_sites; i++) {
    vec3d p(vertices[i]);
    u[3 * i + 0] = a * p[0];
    u[3 * i + 1] = b * p[1];
    u[3 * i + 2] = c * p[2];
  }
  ops.calculate_divergence(u.data(), div_u.data());
  for (int i = 0; i < n_sites; i++) {
    UT_ASSERT(div_u[i] < 1e20);
    vec3d x(vertices[i]);
    // LOG << div_u[i];
    double du =
        (1 - x[0] * x[0]) * a + (1 - x[1] * x[1]) * b + (1 - x[2] * x[2]) * c;
    // du = a + b + c;
    double ei = std::fabs(du - div_u[i]);
    if (ei < emin) emin = ei;
    if (ei > emax) emax = ei;
  }
  LOG << fmt::format("emin = {}, emax = {}", emin, emax);

  LOG << fmt::format("writing {} polygons", voronoi.polygons().n());
  if (voronoi.polygons().n() > 0) meshb::write(voronoi, "sphere.meshb");
}
UT_TEST_CASE_END(test_sphere)

UT_TEST_SUITE_END(voronoi_operators_test_suite)