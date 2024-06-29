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
#include "fem.h"

#include "io.h"
#include "library.h"
#include "math/vec.hpp"
#include "quadrature.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(fem_test_suite)

UT_TEST_CASE(square_poisson_solver_test) {
  std::vector<int> n = {10, 100, 250, 500};
  std::vector<double> h(n.size());
  std::vector<double> e(n.size());

  for (size_t k = 0; k < n.size(); k++) {
    Grid<Triangle> mesh({n[k], n[k]});

    GeneralPoissonSolver<Triangle> solver(mesh.vertices(), mesh.triangles());
    solver.setup();
    solver.set_force([](const vec3d& p) {
      return 2.0 * (M_PI * M_PI) * sin(M_PI * p[0]) * sin(M_PI * p[1]);
    });

    auto u_exact = [](const vec3d& x) {
      return sin(M_PI * x[0]) * sin(M_PI * x[1]);
    };
    solver.bcs().set_dirichlet_bcs(u_exact);

    PoissonSolverOptions options;
    solver.solve(options);
    if (n[k] <= 100) {
      solver.write("square");
      meshb::write(mesh, "square.meshb");
    }

    double error = solver.calculate_error(u_exact);
    LOG << fmt::format("nt = {} error = {}", mesh.triangles().n(), error);
    e[k] = error;
    h[k] = std::sqrt(mesh.vertices().n());
  }

  auto m = n.size() - 1;
  double slope = -std::log(e[m] / e[m - 1]) / std::log(h[m] / h[m - 1]);
  LOG << fmt::format("slope = {}", slope);
  UT_ASSERT_NEAR(slope, 2.0, 0.05);
}
UT_TEST_CASE_END(square_poisson_solver_test)

UT_TEST_CASE(potential_flow_test) {
  const std::string name = "airfoil";  //"circle-square";
  bool external = name == "airfoil" || name == "circle";
  double R = 0.01;
  int n = 40;
  Squircle mesh(R, 2 * n, n, true);

  if (name == "airfoil" || name == "circle") {
    mesh.vertices().clear();
    mesh.triangles().clear();
    mesh.lines().clear();
    meshb::read(name + ".meshb", mesh);
    R = 0.4;
  }

  double uinf = 2.0;
  PotentialFlowSolver<Triangle> solver(mesh.vertices(), mesh.triangles(), uinf);
  if (false && external) {
    solver.bcs().read(name + ".bc");
    mesh.lines().clear();
    for (const auto& [edge, bnd] : solver.bcs().bc_map()) {
      uint32_t line[2] = {edge.first, edge.second};
      size_t ne = mesh.lines().n();
      mesh.lines().add(line);
      mesh.lines().set_group(ne, bnd);
    }
    solver.bcs().bc_map().clear();
    solver.bcs().import(mesh.lines());

  } else {
    solver.bcs().import(mesh.lines());
  }

  PoissonSolverOptions options;
  options.max_linear_solver_iterations = 1000;
  solver.solve(options);
  solver.write("potential");
  meshb::write(mesh, "potential.meshb");

  auto u_exact = [&](const vec3d& x) {
    double r = length(x);
    double theta = atan2(x[1], x[0]);
    return uinf * r * (1 + R * R / (r * r)) * cos(theta);
  };

  // find a point on the circle and translate the solution
  double offset = std::numeric_limits<double>::max();
  for (size_t k = 0; k < mesh.vertices().n(); k++) {
    vec3d p(mesh.vertices()[k]);
    if (std::fabs(length(p) - R) < 1e-10) {
      offset = solver.solution()[k] - u_exact(p);
    }
  }
  LOG << fmt::format("offset = {}", offset);
  for (size_t k = 0; k < mesh.vertices().n(); k++) {
    solver.solution()[k] -= offset;
  }

  double error = solver.calculate_error(u_exact);
  LOG << fmt::format("error = {}", error);
}
UT_TEST_CASE_END(potential_flow_test)

UT_TEST_CASE(sphere_test) {
  SubdividedIcosahedron mesh(6);
  LOG << fmt::format("# triangles = {}", mesh.triangles().n());

  GeneralPoissonSolver<Triangle> solver(mesh.vertices(), mesh.triangles());
  solver.setup();
  double t = 1;
  double a = 6;
  solver.set_force(
      [t, a](const vec3d& p) { return a * std::exp(-a * t) * p[0] * p[1]; });
  auto u_exact = [t, a](const vec3d& p) {
    return std::exp(-a * t) * p[0] * p[1];
  };

  PoissonSolverOptions options;
  options.need_gradient = false;
  options.tol = 1e-10;
  Timer timer;
  timer.start();
  solver.solve(options);
  timer.stop();
  if (mesh.vertices().n() < 1e7) {
    solver.write("sphere");
    meshb::write(mesh, "sphere.meshb");
  }

  double error = solver.calculate_error(u_exact);
  LOG << fmt::format("error = {}, time = {} s.", error, timer.seconds());
}
UT_TEST_CASE_END(sphere_test)

UT_TEST_SUITE_END(fem_test_suite)