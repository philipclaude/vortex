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
#include "shallow_water.h"

#include <cmath>
#include <filesystem>

#include "io.h"
#include "library.h"
#include "particles.h"
#include "tester.h"
#include "util.h"
#include "voronoi.h"

using namespace vortex;

UT_TEST_SUITE(shallow_water_test_suite)

UT_TEST_CASE(test1) {
#if VORTEX_FULL_UNIT_TEST != 0
  // skip this test when running the full suite of tests
  return;
#endif

  static const int dim = 3;
  typedef SphereDomain Domain_t;
  Domain_t domain;

#if 0
  size_t n_sites = 1e4;
  int n_smooth = 100;
  // initialize random points
  std::vector<coord_t> sites(n_sites * dim, 0.0);
  for (size_t k = 0; k < n_sites; k++) {
    vec3d point = domain.random_point();
    for (int d = 0; d < 3; d++) sites[k * dim + d] = point[d];
  }
#else
  int n_smooth = 0;
  SubdividedSphere<Icosahedron> mesh(5);
  size_t n_sites = mesh.vertices().n();
  const auto& sites = mesh.vertices().data();
  LOG << fmt::format("n_sites = {}", n_sites);
#endif

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

  // set up the fluid simulator
  WilliamsonCase6 test_case;
  test_case.conserve_mass = true;
  test_case.add_artificial_viscosity = true;
  test_case.spring_stiffness = 0;
  ShallowWaterSimulation<Domain_t> solver(domain, n_sites, vertices[0],
                                          vertices.dim(), test_case);

  // set up simulation
  SimulationOptions solver_opts;
  solver_opts.save_initial_mesh = true;
  solver_opts.n_smoothing_iterations = n_smooth;
  solver_opts.advect_from_centroid = true;
  solver.initialize(domain, solver_opts);
  solver.setup();

  std::string output_dir = "swe";
  std::string prefix = output_dir + "/particles";
  size_t n_removed = std::filesystem::remove_all(output_dir);
  LOG << fmt::format("removed {} files with prefix {}", n_removed, prefix);
  std::filesystem::create_directories(output_dir);

  // step in time
  int days = 15;
  solver_opts.time_step = 60;
  solver_opts.verbose = false;
  solver_opts.backtrack = false;
  solver_opts.restart_zero_weights = true;
  solver_opts.max_iter = 5;
  solver_opts.skip_initial_calculation = true;
  double seconds = 0;
  double hour = 0;
  solver.save(prefix + "0.vtk");
  while (seconds <= days_to_seconds(days)) {
    double dt = solver.forward_euler_step(solver_opts);
    solver_opts.iteration++;
    seconds += dt;
    solver_opts.time = seconds;
    int current_hour = seconds / 3600;
    if (current_hour == hour + 1) {
      std::string filename = fmt::format("{}{}.vtk", prefix, ++hour);
      solver.save(filename);
    }
  }
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(shallow_water_test_suite)