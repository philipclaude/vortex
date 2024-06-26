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
  void read(const std::string& filename);
  auto& bc_map() { return bc_map_; }
  const auto& bc_map() const { return bc_map_; }

 private:
  std::unordered_map<std::pair<uint32_t, uint32_t>, uint32_t> bc_map_;
};

class PoissonSolver {
 public:
  using Triangles_t = Topology<Triangle>;
  PoissonSolver(const Vertices& vertices, const Triangles_t& triangles)
      : vertices_(vertices),
        triangles_(triangles),
        laplacian_(vertices.n(), vertices.n()),
        rhs_(vertices.n()),
        sol_(vertices.n()),
        grad_sol_(triangles.n()) {}

  virtual ~PoissonSolver() {}

  virtual void apply_boundary_conditions() = 0;

  template <typename Element_t>
  void build();

  template <typename Element_t>
  void solve() {
    build<Element_t>();
    apply_boundary_conditions();
    laplacian_.solve_nl(rhs_, sol_, 1e-10, true);
    calculate_solution_gradient<Element_t>();
  }

  template <typename Element_t>
  void calculate_solution_gradient();

  void write(const std::string& prefix) const;

 protected:
  const Vertices& vertices_;
  const Triangles_t& triangles_;
  spmat<double> laplacian_;
  vecd<double> rhs_;
  vecd<double> sol_;
  vecd<vec3d> grad_sol_;
};

class PotentialFlowSolver : public PoissonSolver {
 public:
  PotentialFlowSolver(const Vertices& vertices, const Triangles_t& triangles,
                      double uinf);
  void apply_boundary_conditions();

  auto& bcs() { return bcs_; }

 private:
  double uinf_{2.0};
  BoundaryConditions bcs_;
};

}  // namespace vortex