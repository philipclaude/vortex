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
  return;
#endif

  static const int dim = 3;
  size_t n_sites = 5e4;
  std::vector<coord_t> sites(n_sites * dim, 0.0);
#if 1
  for (size_t k = 0; k < n_sites; k++) {
    sites[k * dim + 0] = double(rand()) / double(RAND_MAX);
    sites[k * dim + 1] = double(rand()) / double(RAND_MAX);
    sites[k * dim + 2] = 0.0;
  }
  typedef SquareDomain Domain_t;
#else
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

  // smooth the initial point distribution with Lloyd relaxation
  VoronoiDiagram smoother(dim, vertices[0], n_sites);
  VoronoiDiagramOptions options;
  options.n_neighbors = 100;
  options.allow_reattempt = false;
  options.parallel = true;
  options.store_facet_data = true;
  int n_iter = 100;

  for (int iter = 1; iter <= n_iter; ++iter) {
    options.store_mesh = false;
    options.verbose = false;
    smoother.compute(domain, options);  // calculate voronoi diagram
    smoother.smooth(vertices, false);   // move sites to centroids
  }

  // set up the fluid simulator
  SpringParticles<Domain_t> solver(domain, n_sites, vertices[0],
                                   vertices.dim());

  // assign initial velocities
  auto velocity = [](const double* x) -> vec3 { return {0, 0, 0}; };
  solver.particles().set_velocity(velocity);

  auto density = [](const double* x) -> double {
    double f = 0.1 * sin(3 * M_PI * x[0]);
    return x[1] - 0.5 > f ? 10 : 1;
  };
  solver.particles().set_density(density);

  // set up the forcing function
  auto fext = [](const Particle& p) -> vec3d {
    vec3d omega{0., 1., 0.};
    vec3d f = -2 * p.mass * cross(omega, p.velocity);
    // f[1] -= 9.81 * p.mass;
    return f;
  };

  // set up fluid properties
  FluidProperties props;
  // props.viscosity = 0;

  SimulationOptions solver_opts;
  solver.initialize(domain, solver_opts);
  int nt = 100000;
  solver_opts.time_step = 3e-5;
  solver_opts.verbose = false;
  int s = 0;
  for (int t = 0; t < nt; t++) {
    if (t % 50 == 0)
      solver.particles().save(fmt::format("particles10k/particles{}.vtk", s++));
    solver.step(fext, props, solver_opts);
    solver_opts.time += solver_opts.time_step;
    solver_opts.iteration++;
  }
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(power_particles_test_suite)