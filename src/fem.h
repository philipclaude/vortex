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
#pragma once

#include <unordered_map>

#include "defs.h"
#include "elements.h"
#include "math/spmat.h"
#include "mesh.h"

namespace vortex {

class BoundaryConditions {
 public:
  using DirichletFunction = std::function<double(const vec3d& x)>;
  BoundaryConditions() {
    dirichlet_ = [](const vec3d& x) { return 0.0; };
  }
  void read(const std::string& filename);
  void import(const Topology<Line>& lines);
  auto& bc_map() { return bc_map_; }
  const auto& bc_map() const { return bc_map_; }

  void set_dirichlet_bcs(const DirichletFunction& f) { dirichlet_ = f; }
  const auto& dirichlet() const { return dirichlet_; }

 private:
  std::unordered_map<std::pair<uint32_t, uint32_t>, uint32_t> bc_map_;
  DirichletFunction dirichlet_;
};

struct PoissonSolverOptions {
  double tol{1e-10};
  bool need_gradient{true};
  bool has_rhs{true};
  int max_linear_solver_iterations{100};
};

class PoissonSolverBase {
 public:
  using Triangles_t = Topology<Triangle>;
  PoissonSolverBase(const Vertices& vertices, const Triangles_t& triangles)
      : vertices_(vertices),
        triangles_(triangles),
        laplacian_(vertices.n(), vertices.n()),
        rhs_(vertices.n()),
        sol_(vertices.n()),
        grad_sol_(triangles.n()) {}

  const auto& laplacian() const { return laplacian_; }
  auto& solution() { return sol_; }
  const auto& solution() const { return sol_; }
  void write(const std::string& prefix) const;

 protected:
  const Vertices& vertices_;
  const Triangles_t& triangles_;
  spmat<double> laplacian_;
  vecd<double> rhs_;
  vecd<double> sol_;
  vecd<vec3d> grad_sol_;
};

template <typename Element_t>
class PoissonSolver : public PoissonSolverBase {
 public:
  using Triangles_t = Topology<Triangle>;

  PoissonSolver(const Vertices& vertices, const Triangles_t& triangles)
      : PoissonSolverBase(vertices, triangles) {}

  virtual ~PoissonSolver() {}
  virtual void set_rhs() = 0;

  void build();

  void solve(PoissonSolverOptions opts) {
    build();
    if (opts.has_rhs) set_rhs();
    SparseSolverOptions linear_opts;
    linear_opts.tol = opts.tol;
    linear_opts.symmetric = true;
    linear_opts.max_iterations = opts.max_linear_solver_iterations;
    laplacian_.solve_nl(rhs_, sol_, linear_opts);
    if (opts.need_gradient) calculate_solution_gradient();
  }

  void calculate_solution_gradient();

  // calculates L2 error between the computed and exact solution
  using ExactSolution = std::function<double(const vec3d& x)>;
  double calculate_error(const ExactSolution& u_exact) const;
  double calculate_error_rms(const ExactSolution& u_exact) const;
};

template <typename Element_t>
class PotentialFlowSolver : public PoissonSolver<Element_t> {
  using Base_t = PoissonSolver<Element_t>;
  using Triangles_t = typename Base_t::Triangles_t;
  using Base_t::laplacian_;
  using Base_t::rhs_;
  using Base_t::sol_;
  using Base_t::triangles_;
  using Base_t::vertices_;

 public:
  PotentialFlowSolver(const Vertices& vertices, const Triangles_t& triangles,
                      double uinf);
  void set_rhs();
  auto& bcs() { return bcs_; }

 private:
  double uinf_{2.0};
  BoundaryConditions bcs_;
};

template <typename Element_t>
class GeneralPoissonSolver : public PoissonSolver<Element_t> {
  using Base_t = PoissonSolver<Element_t>;
  using Triangles_t = typename Base_t::Triangles_t;
  using Base_t::laplacian_;
  using Base_t::rhs_;
  using Base_t::sol_;
  using Base_t::triangles_;
  using Base_t::vertices_;

 public:
  using ForcingFunction = std::function<double(const vec3d& x)>;
  GeneralPoissonSolver(const Vertices& vertices, const Triangles_t& triangles);
  void setup();
  void set_rhs();
  void set_force(const ForcingFunction& f) { force_ = f; }
  auto& bcs() { return bcs_; }

 private:
  ForcingFunction force_;
  BoundaryConditions bcs_;
};

}  // namespace vortex