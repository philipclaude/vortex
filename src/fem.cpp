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
#include "math/mat.hpp"
#include "math/vec.hpp"
#include "mesh.h"
#include "quadrature.hpp"

namespace vortex {

void BoundaryConditions::import(const Topology<Line>& lines) {
  for (size_t k = 0; k < lines.n(); k++) {
    int e0 = lines[k][0];
    int e1 = lines[k][1];
    int bnd = lines.group(k);
    bc_map_.insert({{e0, e1}, bnd});
  }
  LOG << fmt::format("imported {} boundary edges", bc_map_.size());
}

template <typename Element_t>
void PoissonSolver<Element_t>::build() {
  laplacian_.clear();
  const bool project_gradients_{false};
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

    // project the gradients
    if (project_gradients_) {
      mats<3, 3, double> P;
      P.eye();
      vec3d a(pa);
      vec3d b(pb);
      vec3d c(pc);
      vec3d n = cross(b - a, c - a);
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) P(i, j) -= n[i] * n[j];
      ga = P * ga;
      gb = P * gb;
      gc = P * gc;
    }

    // evaluate \int_{triangle} dv . dw
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
void PoissonSolver<Element_t>::calculate_solution_gradient() {
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

template <typename Element_t>
double PoissonSolver<Element_t>::calculate_error(
    const ExactSolution& u_exact) const {
  TriangleQuadrature<Triangle> quad(4);
  double error = 0;
  for (size_t k = 0; k < triangles_.n(); k++) {
    const auto* t = triangles_[k];
    const auto* pa = vertices_[t[0]];
    const auto* pb = vertices_[t[1]];
    const auto* pc = vertices_[t[2]];
    double ua = sol_[t[0]];
    double ub = sol_[t[1]];
    double uc = sol_[t[2]];

    double integral = 0;
    for (size_t k = 0; k < quad.points().size(); k++) {
      const vec2d& qk = quad.points()[k];
      vec3d point = Triangle::get_physical_coordinates(pa, pb, pc, &qk[0]);
      double ue = u_exact(point);

      vec3d basis;
      Triangle::get_basis(qk[0], qk[1], &basis[0]);
      double uh = basis[0] * ua + basis[1] * ub + basis[2] * uc;

      integral += quad.weights()[k] * (ue - uh) * (ue - uh);
    }
    double dj = 2.0 * Element_t::area(pa, pb, pc);

    error += integral * dj;
  }
  return std::sqrt(error);
}

template <typename Element_t>
double PoissonSolver<Element_t>::calculate_error_rms(
    const ExactSolution& u_exact) const {
  double error = 0;
  for (size_t k = 0; k < sol_.m(); k++) {
    vec3d p(vertices_[k]);
    double ue = u_exact(p);
    error += (sol_[k] - ue) * (sol_[k] - ue);
  }
  return std::sqrt(error / sol_.m());
}

void PoissonSolverBase::write(const std::string& prefix) const {
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

template <typename Element_t>
PotentialFlowSolver<Element_t>::PotentialFlowSolver(
    const Vertices& vertices, const Triangles_t& triangles, double uinf)
    : PoissonSolver<Element_t>(vertices, triangles), uinf_(uinf) {}

template <typename Element_t>
void PotentialFlowSolver<Element_t>::set_rhs() {
  for (const auto& [edge, bnd] : bcs_.bc_map()) {
    auto e0 = edge.first;
    auto e1 = edge.second;

    vec3d x0(vertices_[e0]);
    vec3d x1(vertices_[e1]);
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

template <typename Element_t>
GeneralPoissonSolver<Element_t>::GeneralPoissonSolver(
    const Vertices& vertices, const Triangles_t& triangles)
    : PoissonSolver<Element_t>(vertices, triangles) {
  force_ = [](const vec3d& x) { return 0.0; };
}

template <typename Element_t>
void GeneralPoissonSolver<Element_t>::setup() {
  std::unordered_set<std::pair<uint32_t, uint32_t>> edges;
  for (size_t k = 0; k < triangles_.n(); k++) {
    const auto* t = triangles_[k];
    for (int j = 0; j < 3; j++) {
      uint32_t p = t[j];
      uint32_t q = (j == 2) ? t[0] : t[j + 1];
      if (p > q) std::swap(p, q);
      auto it = edges.find({p, q});
      if (it == edges.end()) {
        edges.insert({p, q});
      } else
        edges.erase(it);
    }
  }

  const int bc_id = 1;
  for (auto& e : edges) bcs_.bc_map().insert({e, bc_id});
}

template <typename Element_t>
void GeneralPoissonSolver<Element_t>::set_rhs() {
  // integrate the forcing term over the mesh: \int_{mesh} basis * f
  TriangleQuadrature<Triangle> quad(4);  // TODO user-option for quad order
  for (size_t k = 0; k < triangles_.n(); k++) {
    const auto* t = triangles_[k];
    const auto* pa = vertices_[t[0]];
    const auto* pb = vertices_[t[1]];
    const auto* pc = vertices_[t[2]];

    // integrate over the triangle
    vec3d integral{0, 0, 0};
    vec3d basis;
    for (size_t k = 0; k < quad.points().size(); k++) {
      const auto& qk = quad.points()[k];
      vec3d point = Element_t::get_physical_coordinates(pa, pb, pc, &qk[0]);
      Triangle::get_basis(qk[0], qk[1], &basis[0]);
      integral = integral + quad.weights()[k] * basis * force_(point);
    }
    double dj = 2.0 * Element_t::area(pa, pb, pc);  // jacobian determinant
    integral = integral * dj;

    // add the contribution to the rhs
    for (int j = 0; j < 3; j++) rhs_[t[j]] += integral[j];
  }

  // set the dirichlet boundary conditions
  for (const auto& [edge, _] : bcs_.bc_map()) {
    auto e0 = edge.first;
    auto e1 = edge.second;

    // zero out the row for this DOF, then set diagonal term to 1
    laplacian_.rows()[e0].clear();
    laplacian_(e0, e0) = 1;
    laplacian_.rows()[e1].clear();
    laplacian_(e1, e1) = 1;

    // evaluate and set the dirichlet BC
    vec3d p0(vertices_[e0]);
    vec3d p1(vertices_[e1]);
    rhs_[e0] = bcs_.dirichlet()(p0);
    rhs_[e1] = bcs_.dirichlet()(p1);
  }
}

template class PoissonSolver<Triangle>;
template class PotentialFlowSolver<Triangle>;
template class GeneralPoissonSolver<Triangle>;

}  // namespace vortex