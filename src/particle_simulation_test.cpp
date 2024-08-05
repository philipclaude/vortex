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
#include <cmath>

#include "graphics.h"
#include "io.h"
#include "library.h"
#include "math/vec.hpp"
#include "particles.h"
#include "tester.h"
#include "util.h"
#include "voronoi.h"

using namespace vortex;

UT_TEST_SUITE(power_particles_test_suite)

UT_TEST_CASE(test1) {
#if VORTEX_FULL_UNIT_TEST != 0
  // skip this test when running the full suite of tests
  return;
#endif

  static const int dim = 3;
  size_t n_sites = 1e5;
  std::vector<coord_t> sites(n_sites * dim, 0.0);
#if 0  // square
  for (size_t k = 0; k < n_sites; k++) {
    sites[k * dim + 0] = double(rand()) / double(RAND_MAX);
    sites[k * dim + 1] = double(rand()) / double(RAND_MAX);
    sites[k * dim + 2] = 0.0;
  }
  typedef SquareDomain Domain_t;
  auto velocity = [](const double* x) -> vec3 { return {0, 0, 0}; };
  auto fext = [](const Particle& p) -> vec3d {
    vec3d f;
    f[1] -= 9.81 * p.mass;
    return f;
  };
  auto density = [](const double* x) -> double {
    double f = 0.1 * sin(3 * M_PI * x[0]);
    return x[1] - 0.5 > f ? 10 : 1;
  };
#else  // sphere
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
  vec3d omega{0., 0.01, 0.};
  auto velocity = [omega](const double* x) -> vec3d {
    vec3d p(x);
    vec3d v = cross(omega, p);
    return v;
  };
  auto fext = [omega](const Particle& particle) -> vec3d {
    vec3d f;
    double m = particle.mass;
    const auto& v = particle.velocity;
    const auto& p = particle.position;
    // f = -2 * m * cross(omega, v);
    f = -2 * m * omega[1] * p[1] * cross(p, v);
    return f;
  };
  auto density = [](const double* x) -> double {
    return std::fabs(x[1]) > 0.25 ? 1 : 10;
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
  SpringParticles<Domain_t> solver(domain, n_sites, vertices[0],
                                   vertices.dim());

  // assign initial velocities
  solver.particles().set_velocity(velocity);
  solver.particles().set_density(density);

  // set up fluid properties
  FluidProperties props;
  // props.viscosity = 0;

  SimulationOptions solver_opts;
  solver_opts.save_initial_mesh = true;
  solver_opts.n_smoothing_iterations = 100;
  solver.initialize(domain, solver_opts);
  // std::vector<double> target_vol(n_sites, 4.0 * M_PI / n_sites);
  // solver.optimize_volumes(domain, solver_opts, target_vol);
  meshb::write(solver.voronoi(), "voronoi.meshb");

  int nt = 50000;
  double hn = solver.voronoi().max_radius();
  // solver_opts.epsilon = 4e-2;
  solver_opts.epsilon = 10 * hn;
  // solver_opts.time_step = 2e-4;
  solver_opts.time_step = 0.15 * std::pow(solver_opts.epsilon, 2);
  LOG << fmt::format("hn = {:1.3e}, eps = {:1.3e}, dt = {:1.3e}", hn,
                     solver_opts.epsilon, solver_opts.time_step);
  solver_opts.verbose = false;
  solver_opts.backtrack = false;
  int s = 0;
  for (int t = 0; t < nt; t++) {
    if (t % 50 == 0)
      solver.particles().save(fmt::format("sphere100k/particles{}.vtk", s++));
    solver.step(fext, props, solver_opts);
    solver_opts.iteration++;
  }
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(power_particles_test_suite)