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
#include "particles.h"

#include <cmath>

#include "graphics.h"
#include "io.h"
#include "library.h"
#include "math/vec.hpp"
#include "tester.h"
#include "util.h"
#include "voronoi.h"

using namespace vortex;

UT_TEST_SUITE(particles_test_suite)

UT_TEST_CASE(test_square_uniform) {
  // randomly initialize points in a square
  static const int dim = 3;
  size_t n_sites = 1e4;
  std::vector<coord_t> sites(n_sites * dim, 0.0);
  for (size_t k = 0; k < n_sites; k++) {
    sites[k * dim + 0] = double(rand()) / double(RAND_MAX);
    sites[k * dim + 1] = double(rand()) / double(RAND_MAX);
    sites[k * dim + 2] = 0.0;
  }

  // construct a better ordering of the points
  std::vector<index_t> order(n_sites);
  sort_points_on_zcurve(sites.data(), n_sites, dim, order);
  Vertices vertices(dim);
  vertices.reserve(n_sites);
  coord_t x[dim];
  for (size_t i = 0; i < n_sites; i++) {
    for (int d = 0; d < dim; d++) x[d] = sites[dim * order[i] + d];
    vertices.add(x);
  }

  // initialize the sphere domain
  SquareDomain domain;

  // smooth the initial point distribution with Lloyd relaxation
  VoronoiDiagram smoother(dim, vertices[0], n_sites);
  VoronoiDiagramOptions options;
  options.n_neighbors = 50;
  options.parallel = true;
  options.store_facet_data = false;
  options.neighbor_algorithm = NearestNeighborAlgorithm::kKdtree;
  int n_iter = 10;
  for (int iter = 1; iter <= n_iter; ++iter) {
    options.store_mesh = false;
    options.verbose = false;
    smoother.compute(domain, options);  // calculate voronoi diagram
    smoother.smooth(vertices, false);   // move sites to centroids
  }

  // set up the target mass for each cell (uniform)
  double v_target = 1.0 / n_sites;
  ParticleSimulation particles(n_sites, vertices[0], vertices.dim());
  particles.voronoi().weights().resize(n_sites, 0);
  std::vector<double> target_volume(n_sites, v_target);
  SimulationOptions sim_opts;
  sim_opts.volume_grad_tol = 1e-14;
  auto result = particles.optimize_volumes(domain, sim_opts, target_volume);
  UT_ASSERT(result.converged);

  // recompute the voronoi diagram to save it
  auto& voronoi = particles.voronoi();
  options.store_mesh = true;
  voronoi.compute(domain, options);

  // check the mass
  const auto& props = voronoi.properties();
  double total_volume = 0.0;
  for (auto& p : props) {
    UT_ASSERT_NEAR(p.volume, v_target, 1e-8);
    total_volume += p.volume;
  }
  LOG << props.size();
  LOG << fmt::format("total volume = {}", total_volume);

  // save the mesh
  LOG << fmt::format("writing {} polygons", voronoi.polygons().n());
  if (voronoi.polygons().n() > 0) meshb::write(voronoi, "particles0.meshb");
}
UT_TEST_CASE_END(test_square_uniform)

UT_TEST_CASE(test_sphere_uniform) {
  // randomly initialize points on the sphere
  auto irand = [](int min, int max) {
    return min + double(rand()) / (double(RAND_MAX) + 1.0) * (max - min);
  };
  static const int dim = 3;
  size_t n_sites = 1e4;
  std::vector<coord_t> sites(n_sites * dim, 0.0);
  for (size_t k = 0; k < n_sites; k++) {
    coord_t theta = 2.0 * M_PI * irand(0, 1);
    coord_t phi = acos(2.0 * irand(0, 1) - 1.0);
    sites[k * dim + 0] = cos(theta) * sin(phi);
    sites[k * dim + 1] = sin(theta) * sin(phi);
    sites[k * dim + 2] = cos(phi);
  }

  // construct a better ordering of the points
  std::vector<index_t> order(n_sites);
  sort_points_on_zcurve(sites.data(), n_sites, dim, order);
  Vertices vertices(dim);
  vertices.reserve(n_sites);
  coord_t x[dim];
  for (size_t i = 0; i < n_sites; i++) {
    for (int d = 0; d < dim; d++) x[d] = sites[dim * order[i] + d];
    vertices.add(x);
  }

  // initialize the sphere domain
  SphereDomain domain;
  domain.set_initialization_fraction(0.7);

  // smooth the initial point distribution with Lloyd relaxation
  VoronoiDiagram smoother(dim, vertices[0], n_sites);
  VoronoiDiagramOptions options;
  options.n_neighbors = 100;
  options.parallel = true;
  options.store_facet_data = false;
  options.neighbor_algorithm = NearestNeighborAlgorithm::kSphereQuadtree;
  int n_iter = 10;
  for (int iter = 1; iter <= n_iter; ++iter) {
    // vtk::write(vertices, fmt::format("particles/points{}.vtk", iter));
    options.store_mesh = false;
    options.verbose = false;
    smoother.compute(domain, options);  // calculate voronoi diagram
    smoother.smooth(vertices, true);    // move sites to centroids
  }

  // set up the target mass for each cell (uniform)
  double v_target = 4.0 * M_PI / n_sites;
  ParticleSimulation particles(n_sites, vertices[0], vertices.dim());
  particles.voronoi().weights().resize(n_sites, 0);
  std::vector<double> target_volume(n_sites, v_target);
  SimulationOptions sim_opts;
  sim_opts.volume_grad_tol = 1e-14;
  auto result = particles.optimize_volumes(domain, sim_opts, target_volume);
  UT_ASSERT(result.converged);

  // recompute the voronoi diagram to save it
  auto& voronoi = particles.voronoi();
  options.store_mesh = true;
  voronoi.compute(domain, options);

  // check the mass
  const auto& props = voronoi.properties();
  double total_volume = 0.0;
  for (auto& p : props) {
    UT_ASSERT_NEAR(p.volume, v_target, 1e-8);
    total_volume += p.volume;
  }
  LOG << fmt::format("total volume = {}", total_volume);

  // save the mesh
  LOG << fmt::format("writing {} polygons", voronoi.polygons().n());
  if (voronoi.polygons().n() > 0) meshb::write(voronoi, "particles1.meshb");
}
UT_TEST_CASE_END(test_sphere_uniform)

UT_TEST_CASE(test_sphere_nonuniform) {
  // randomly initialize points on the sphere
  auto irand = [](int min, int max) {
    return min + double(rand()) / (double(RAND_MAX) + 1.0) * (max - min);
  };
  static const int dim = 3;
  size_t n_sites = 1e4;
#if VORTEX_FULL_UNIT_TEST != 0
  n_sites = 1e4;
#endif
  std::vector<coord_t> sites(n_sites * dim, 0.0);
  for (size_t k = 0; k < n_sites; k++) {
    coord_t theta = 2.0 * M_PI * irand(0, 1);
    coord_t phi = acos(2.0 * irand(0, 1) - 1.0);
    sites[k * dim + 0] = cos(theta) * sin(phi);
    sites[k * dim + 1] = sin(theta) * sin(phi);
    sites[k * dim + 2] = cos(phi);
  }

  // construct a better ordering of the points
  std::vector<index_t> order(n_sites);
  sort_points_on_zcurve(sites.data(), n_sites, dim, order);
  Vertices vertices(dim);
  vertices.reserve(n_sites);
  coord_t x[dim];
  for (size_t i = 0; i < n_sites; i++) {
    for (int d = 0; d < dim; d++) x[d] = sites[dim * order[i] + d];
    vertices.add(x);
  }

  // initialize the sphere domain
  SphereDomain domain;
  domain.set_initialization_fraction(0.7);

  // smooth the initial point distribution with Lloyd relaxation
  VoronoiDiagram smoother(dim, vertices[0], n_sites);
  VoronoiDiagramOptions options;
  options.n_neighbors = 100;
  options.parallel = true;
  options.store_facet_data = false;
  int n_iter = 10;
  std::vector<double> target_volume(n_sites);
  for (int iter = 1; iter <= n_iter; ++iter) {
    options.store_mesh = false;
    options.verbose = false;
    smoother.compute(domain, options);    // calculate voronoi diagram
    for (size_t k = 0; k < n_sites; k++)  // save target mass
      target_volume[k] = smoother.properties()[k].volume;
    smoother.smooth(vertices, true);  // move sites to centroids
    // vtk::write(vertices, fmt::format("particles/points{}.vtk", iter));
  }

  // set up the target mass for each cell from the previous voronoi diagram
  ParticleSimulation particles(n_sites, vertices[0], vertices.dim());
  particles.voronoi().weights().resize(n_sites, 0);
  SimulationOptions sim_opts;
  sim_opts.volume_grad_tol = 1e-14;
  sim_opts.max_iter = 30;
  sim_opts.backtrack = true;
  auto result = particles.optimize_volumes(domain, sim_opts, target_volume);
  UT_ASSERT(result.converged);

  // recompute the voronoi diagram to save it
  auto& voronoi = particles.voronoi();
  options.store_mesh = true;
  voronoi.compute(domain, options);

  // check the mass
  const auto& props = voronoi.properties();
  double total_volume = 0.0;
  for (size_t k = 0; k < n_sites; k++) {
    UT_ASSERT_NEAR(props[k].volume, target_volume[k], 1e-8);
    total_volume += props[k].volume;
  }
  LOG << fmt::format("total volume = {}", total_volume);

  // save the mesh
  LOG << fmt::format("writing {} polygons", voronoi.polygons().n());
  if (voronoi.polygons().n() > 0) meshb::write(voronoi, "particles2.meshb");
}
UT_TEST_CASE_END(test_sphere_nonuniform)

UT_TEST_SUITE_END(particles_test_suite)