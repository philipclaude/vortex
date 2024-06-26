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

#include <fstream>

#include "io.h"
#include "math/vec.hpp"
#include "mesh.h"

namespace vortex {

void BoundaryConditions::read(const std::string& filename) {
  std::ifstream file(filename);
  std::string line;
  while (!file.eof()) {
    std::getline(file, line);
    int e0, e1, bnd;
    int n_items = sscanf(line.c_str(), "%d %d %d", &e0, &e1, &bnd);
    if (n_items <= 0) continue;
    bc_map_.insert({{e0, e1}, bnd});
  }
  LOG << fmt::format("read {} boundary edges", bc_map_.size());
}

template <typename Element_t>
void PoissonSolver::build() {
  laplacian_.clear();
  for (size_t k = 0; k < triangles_.n(); k++) {
    const auto* t = triangles_[k];
    const auto* pa = vertices_[t[0]];
    const auto* pb = vertices_[t[1]];
    const auto* pc = vertices_[t[2]];

    // TODO(philip) use quadrature if the basis functions are not linear
    const coord_t u = Element_t::center[0];
    const coord_t v = Element_t::center[1];
    vec3d ga, gb, gc;
    Element_t::get_basis_gradient(u, v, pa, pb, pc, ga, gb, gc);
    double area = Element_t::area(pa, pb, pc);

    mats<3, 3, double> M;
    M(0, 0) = dot(ga, ga);
    M(0, 1) = dot(ga, gb);
    M(0, 2) = dot(ga, gc);
    M(1, 0) = M(0, 1);
    M(1, 1) = dot(gb, gb);
    M(1, 2) = dot(gb, gc);
    M(2, 0) = M(0, 2);
    M(2, 1) = M(1, 2);
    M(2, 2) = dot(gc, gc);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) laplacian_(t[i], t[j]) += M(i, j) * area;
  }
}

template <typename Element_t>
void PoissonSolver::calculate_solution_gradient() {
  for (size_t k = 0; k < triangles_.n(); k++) {
    const auto* t = triangles_[k];
    const auto* pa = vertices_[t[0]];
    const auto* pb = vertices_[t[1]];
    const auto* pc = vertices_[t[2]];

    const coord_t u = Element_t::center[0];
    const coord_t v = Element_t::center[1];
    vec3d ga, gb, gc;
    Element_t::get_basis_gradient(u, v, pa, pb, pc, ga, gb, gc);
    grad_sol_[k] = ga * sol_[t[0]] + gb * sol_[t[1]] + gc * sol_[t[2]];
  }
}

void PoissonSolver::write(const std::string& prefix) const {
  std::vector<std::array<double, 1>> s(sol_.m());
  for (size_t k = 0; k < s.size(); k++) s[k] = {sol_[k]};
  meshb::write_sol<1>(s, true, prefix + ".sol");
  std::vector<std::array<double, 1>> speed(triangles_.n());
  for (size_t k = 0; k < speed.size(); k++) speed[k] = {length(grad_sol_[k])};
  meshb::write_sol<1>(speed, false, prefix + "-speed.sol");

  std::vector<std::array<double, 3>> velocity(sol_.m(), {0, 0, 0});
  std::vector<int> count(sol_.m(), 0);
  for (size_t k = 0; k < triangles_.n(); k++) {
    for (int i = 0; i < 3; i++) {
      count[triangles_[k][i]]++;
      for (int d = 0; d < 3; d++)
        velocity[triangles_[k][i]][d] += grad_sol_[k][d];
    }
  }
  for (size_t k = 0; k < velocity.size(); k++) {
    for (int d = 0; d < 3; d++) velocity[k][d] /= count[k];
  }
  meshb::write_sol<3>(velocity, true, prefix + "-velocity.sol");
}

PotentialFlowSolver::PotentialFlowSolver(const Vertices& vertices,
                                         const Triangles_t& triangles,
                                         double uinf)
    : PoissonSolver(vertices, triangles), uinf_(uinf) {}

void PotentialFlowSolver::apply_boundary_conditions() {
  for (const auto& [edge, bnd] : bcs_.bc_map()) {
    auto e0 = edge.first;
    auto e1 = edge.second;

    vec3d x0(vertices_[e0], 2);
    vec3d x1(vertices_[e1], 2);
    double ds = length(x1 - x0);
    double s = 0;
    if (bnd == 3)
      s = -1;
    else if (bnd == 1)
      s = 1;
    rhs_[e0] += s * 0.5 * uinf_ * ds;
    rhs_[e1] += s * 0.5 * uinf_ * ds;
  }
}

template void PoissonSolver::build<Triangle>();
template void PoissonSolver::calculate_solution_gradient<Triangle>();

}  // namespace vortex