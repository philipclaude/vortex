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
#include "neighbors.h"

#include "io.h"
#include "log.h"
#include "mesh.h"
#include "tester.h"
#include "voronoi.h"

using namespace vortex;

UT_TEST_SUITE(neighbors_test_suite)

UT_TEST_CASE(test1) {
  SphereDomain domain;
  static const int dim = 4;
  size_t n_sites = 5e4;
  std::vector<coord_t> sites(n_sites * dim, 0.0);
  for (size_t k = 0; k < n_sites; k++) {
    auto point = domain.random_point();
    for (int d = 0; d < 3; d++) sites[k * dim + d] = point[d];
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

  VoronoiDiagramOptions options;
  options.n_neighbors = 50;
  options.parallel = true;
  options.store_mesh = false;
  options.verbose = true;

  // build a kdtree for comparison
  std::vector<index_t> knn(n_sites * options.n_neighbors);
  auto tree =
      get_nearest_neighbors<dim>(vertices[0], n_sites, vertices[0], n_sites,
                                 knn, options.n_neighbors, options);

  // calculate the voronoi diagram
  VoronoiDiagram voronoi(dim, vertices[0], n_sites);
  voronoi.compute(domain, options);

  Timer timer;
  timer.start();
  VoronoiNeighbors neighbors(voronoi, vertices[0], dim);
  neighbors.build();
  NearestNeighborsWorkspace search(options.n_neighbors);
  search.max_level = 10;

  for (size_t k = 0; k < n_sites; k++) {
    // check the nearest neighbors match those from the kdtree
    neighbors.knearest(k, search);
    const auto& result = search.sites;

    size_t m = result.size();
    if (result.size() > options.n_neighbors) m = options.n_neighbors;
    for (size_t j = 0; j < m; j++) {
      UT_ASSERT_EQUALS(knn[options.n_neighbors * k + j],
                       result.data()[j].first);
    }
  }
  timer.stop();
  size_t n_avg = search.total_neighbors / n_sites;
  LOG << fmt::format("avg = {}, time = {} sec.", n_avg, timer.seconds());
}
UT_TEST_CASE_END(test1)

UT_TEST_CASE(test2) {
  SphereDomain domain;
  static const int dim = 3;
  size_t n_sites = 1e5;
  std::vector<coord_t> sites(n_sites * dim, 0.0);
  for (size_t k = 0; k < n_sites; k++) {
    auto point = domain.random_point();
    for (int d = 0; d < 3; d++) sites[k * dim + d] = point[d];
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

  VoronoiDiagramOptions options;
  options.n_neighbors = 100;
  options.parallel = true;
  options.store_mesh = false;
  options.verbose = true;

  // build a kdtree for comparison
  std::vector<index_t> knn;
  knn.resize(n_sites * options.n_neighbors);
  auto tree =
      get_nearest_neighbors<dim>(vertices[0], n_sites, vertices[0], n_sites,
                                 knn, options.n_neighbors, options);

  int m = 10;
  int ns = std::log(n_sites / (Octahedron::n_faces * m)) / std::log(4);

  Timer timer;
  timer.start();
  SphereQuadtree neighbors(vertices[0], vertices.n(), dim, ns);
  timer.stop();
  LOG << fmt::format("setup time for ns = {}, nt = {}: {} s.", ns,
                     neighbors.n_triangles(), timer.seconds());

  timer.start();
  neighbors.build();
  timer.stop();
  LOG << fmt::format("build time = {} s.", timer.seconds());
  meshb::write(neighbors.mesh(), "neighbors.meshb");

  timer.start();
  size_t n_threads = std::thread::hardware_concurrency();
  std::vector<SphereQuadtreeWorkspace> searches(n_threads, options.n_neighbors);
  std::vector<size_t> nn(n_sites);
  std::parafor_i(0, n_sites, [&](int tid, size_t k) {
    auto& search = searches[tid];
    neighbors.knearest(k, search);
    UT_ASSERT(search.neighbors.size() > 0);
    const auto& result = search.neighbors;
    size_t r = search.n_neighbors < options.n_neighbors ? search.n_neighbors
                                                        : options.n_neighbors;
    double d = -1;
    for (size_t j = 0; j < r; j++) {
      UT_ASSERT(d <= result[j].second);
      if (r < 10)
        UT_ASSERT_EQUALS(result[j].first, knn[options.n_neighbors * k + j]);
      d = result[j].second;
    }
    nn[k] = search.n_neighbors;
  });
  timer.stop();
  LOG << fmt::format("computed nearest neighbors in {} s.", timer.seconds());

  double n_avg = 0;
  for (auto n : nn) n_avg += n;
  n_avg /= n_sites;
  LOG << fmt::format("average # neighbors = {}", n_avg);
}
UT_TEST_CASE_END(test2)

UT_TEST_SUITE_END(neighbors_test_suite)