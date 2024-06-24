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
#if VORTEX_FULL_UNIT_TEST != 0
  // skip this test when running the full suite of tests
  return;
#endif
  SphereDomain domain;
  static const int dim = 4;
  size_t n_sites = 1e7;
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

  VoronoiDiagram voronoi(dim, vertices[0], n_sites);
  VoronoiDiagramOptions options;
  options.n_neighbors = 100;
  options.parallel = true;
  int n_iter = 10;
  for (int iter = 1; iter <= n_iter; ++iter) {
    options.store_mesh = iter == n_iter;
    options.verbose = true;
    voronoi.vertices().clear();
    voronoi.vertices().set_dim(3);
    voronoi.polygons().clear();
    voronoi.triangles().clear();
    voronoi.compute(domain, options);
    // options.voronoi_neighbors = true;

    // move each site to the centroid of the corresponding cell
    voronoi.smooth(vertices, false);
    auto props = voronoi.analyze();
    LOG << fmt::format("iter = {}, area = {}", iter, props.area);

#if 0
    Timer timer;
    timer.start();
    VoronoiNeighbors neighbors(voronoi, vertices[0]);
    size_t n_threads = std::thread::hardware_concurrency();
    std::vector<NearestNeighborsWorkspace2> searches(n_threads,
                                                     options.n_neighbors);

    std::parafor_i(0, n_sites, [&neighbors, &searches](size_t tid, size_t k) {
      neighbors.knearest2(k, searches[tid]);
    });
    timer.stop();
    size_t n_avg = 0;
    for (int i = 0; i < n_threads; i++) n_avg += searches[i].n_avg;
    n_avg /= n_sites;
    LOG << fmt::format("avg # neighbors = {}, neighbors time = {} sec.", n_avg, timer.seconds());
#endif
  }
  meshb::write(voronoi, "voronoi.meshb");
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(neighbors_test_suite)