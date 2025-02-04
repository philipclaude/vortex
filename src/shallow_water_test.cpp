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
  size_t n_sites = 1e5;
  std::vector<coord_t> sites(n_sites * dim, 0.0);

  auto irand = [](int min, int max) {
    return min + double(rand()) / (double(RAND_MAX) + 1.0) * (max - min);
  };
  for (size_t k = 0; k < n_sites; k++) {
    coord_t theta = 2.0 * M_PI * irand(0, 1);
    coord_t phi = acos(2.0 * irand(0, 1) - 1.0);
    sites[k * dim + 0] = cos(theta) * sin(phi);
    sites[k * dim + 1] = sin(theta) * sin(phi);
    sites[k * dim + 2] = cos(phi);
  }
  typedef SphereDomain Domain_t;
  EarthProperties earth;
  const double omega = earth.angular_velocity;
  const double a = earth.radius;
  const double g = earth.gravity;

#if 1  // case 5
  const double u0 = 20.0;
  auto initial_velocity = [u0](const double* x) -> vec3d {
    return {-u0 * x[1], u0 * x[0], 0.0};
  };

  const double h0 = 5960;
  const double hs0 = 2000;
  const double R = M_PI / 9.0;
  const double lambda_c = -M_PI / 2;
  const double theta_c = M_PI / 6;

  auto surface_height = [&](const double* x) -> double {
    double lambda = atan2(x[1], x[0]);  // atan2 returns in [-pi, pi]
    double z = x[2];
    ASSERT(std::fabs(x[2]) <= 1);
    double theta = asin(z);  // asin returns in [-pi/2, pi/2]
    double d_lambda = lambda - lambda_c;
    double d_theta = theta - theta_c;
    double d = std::min(R * R, d_lambda * d_lambda + d_theta * d_theta);
    double hs = hs0 * (1 - std::pow(d, 0.5) / R);
    // return hs0 * std::exp(-2.8 * 2.8 * d / (R * R));
    //  ASSERT(hs >= 0) << hs;
    //   if (hs > 0) LOG << hs;
    return hs;
  };
  auto initial_height = [&](const double* x) -> double {
    double h = h0 - (omega * a * u0 + 0.5 * u0 * u0) * x[2] * x[2] / g;
    return h;
  };

  auto coriolis_parameter = [&](const double* x) -> double {
    return 2.0 * omega * x[2];
  };
  double hmin = h0 - (omega * a * u0 + 0.5 * u0 * u0) / g;
  LOG << fmt::format("hmin = {}, hmax = {}", hmin, h0);

#elif 1  // case 2
  double h0 = 2.94e4 / g;
  double u0 = 2 * M_PI * a / (12 * 24 * 3600);
  auto surface_height = [](const double* x) -> double { return 0.0; };
  auto initial_height = [&](const double* x) -> double {
    return h0 - (omega * a * u0 + 0.5 * u0 * u0) * x[0] * x[0] / g;
  };
  auto initial_velocity = [&](const double* x) -> vec3d {
    double lambda = atan2(x[1], x[0]);  // atan2 returns in [-pi, pi]
    double z = x[2];
    ASSERT(std::fabs(x[2]) <= 1);
    double theta = asin(z);  // asin returns in [-pi/2, pi/2]
    double us = u0 * sin(theta) * cos(lambda);
    double vs = -u0 * sin(lambda);
    double ux = -us * sin(lambda) - vs * sin(theta) * cos(lambda);
    double uy = us * cos(lambda) - vs * sin(theta) * sin(lambda);
    double uz = vs * cos(theta);
    return {ux, uy, uz};
  };

  auto coriolis_parameter = [&](const double* x) -> double {
    return -2.0 * omega * x[0];
  };
#else    // case 1
  auto surface_height = [](const double* x) -> double { return 0.0; };
  auto initial_height = [&](const double* x) -> double {
    double lambda = atan2(x[1], x[0]);  // atan2 returns in [-pi, pi]
    double z = x[2];
    ASSERT(std::fabs(x[2]) <= 1);
    double theta = asin(z);  // asin returns in [-pi/2, pi/2]
    double R = earth.radius / 3;
    double h0 = 1000;
    double lc = 3 * M_PI / 2;
    double tc = 0;
    double dl = lambda - lc;
    double dt = theta - tc;
    double r = earth.radius *
               std::acos(sin(tc) * sin(theta) + cos(tc) * cos(theta) * cos(dl));
    if (r >= R) return 0;
    return 0.5 * h0 * (1 + cos(M_PI * r / R));
    return 0;
  };
  auto initial_velocity = [&](const double* x) -> vec3d {
    double lambda = atan2(x[1], x[0]);  //+ M_PI;  // atan2 returns in [-pi, pi]
    double z = x[2];
    ASSERT(std::fabs(x[2]) <= 1);
    double theta = asin(z);  // asin returns in [-pi/2, pi/2]
    double u0 = 2 * M_PI * earth.radius / (12 * 24 * 3600);
    double us = u0 * sin(theta) * cos(lambda);
    double vs = -u0 * sin(lambda);
    double ux = -us * sin(lambda) - vs * sin(theta) * cos(lambda);
    double uy = us * cos(lambda) - vs * sin(theta) * sin(lambda);
    double uz = vs * cos(theta);
    return {ux, uy, uz};
  };
  auto coriolis_parameter = [&](const double* x) -> double {
    return 2.0 * omega * x[2];
  };

#endif

  Domain_t domain;

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
  ShallowWaterSimulation<Domain_t> solver(domain, n_sites, vertices[0],
                                          vertices.dim());
  auto& swe_opts = solver.options();
  swe_opts.surface_height = surface_height;
  swe_opts.initial_height = initial_height;
  swe_opts.initial_velocity = initial_velocity;
  swe_opts.analytic_velocity = initial_velocity;
  swe_opts.use_analytic_velocity = false;
  swe_opts.project_points = false;
  swe_opts.project_velocity = false;
  swe_opts.conserve_mass = true;
  swe_opts.smoothing_iterations = 0;
  swe_opts.coriolis_parameter = coriolis_parameter;

  // set up simulation
  SimulationOptions solver_opts;
  solver_opts.save_initial_mesh = true;
  solver_opts.n_smoothing_iterations = 100;
  solver_opts.advect_from_centroid = true;
  solver.initialize(domain, solver_opts);
  meshb::write(solver.voronoi(), "voronoi.meshb");
  solver.setup();
  solver.start();

  // step in time
  int nt = 25000;
  solver_opts.time_step = 60;
  solver_opts.verbose = false;
  solver_opts.backtrack = false;
  solver_opts.time = 0;
  double hour = 0;
  solver.save("swe/particles0.vtk");
  // return;
  for (int t = 0; t < nt; t++) {
    int current_hour = solver_opts.time / 3600;
    if (current_hour == hour + 1) {
      std::string filename = fmt::format("swe/particles{}.vtk", ++hour);
      LOG << "saving: " << filename;
      solver.save(filename);
    }
    solver.forward_euler_step(solver_opts);
    solver_opts.iteration++;
  }
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(shallow_water_test_suite)