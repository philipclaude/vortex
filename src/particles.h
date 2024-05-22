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

#include <vector>

#include "math/spmat.h"
#include "voronoi.h"

namespace vortex {

struct SimulationOptions {
  double volume_grad_tol{1e-8};
  int max_iter{25};
  int n_neighbors{80};
  bool verbose{true};
};

struct SimulationConvergence {
  double error{1};
  bool converged{false};
  int n_iterations{0};
};

class Particles : public Vertices {
 public:
  Particles(size_t np, const coord_t* xp, int dim) : Vertices(4) {
    create(np, xp, dim);
  }

  void create(size_t np, const coord_t* xp, int dim);

  const auto& mass() const { return mass_; }
  auto& mass() { return mass_; }

 private:
  std::vector<double> volume_;
  std::vector<double> mass_;
};

class ParticleSimulation {
 public:
  ParticleSimulation(size_t np, const coord_t* xp, int dim)
      : particles_(np, xp, dim),
        voronoi_(4, particles_[0], np),
        hessian_(np, np),
        dw_(np),
        gradient_(np) {}

  template <typename Domain_t>
  SimulationConvergence conserve_mass(const Domain_t& domain,
                                      SimulationOptions sim_opts) {
    // calculate the voronoi diagram
    VoronoiDiagramOptions voro_opts;
    voro_opts.store_mesh = false;
    voro_opts.store_facet_data = true;
    voro_opts.allow_reattempt = false;
    voro_opts.n_neighbors = sim_opts.n_neighbors;
    voro_opts.verbose = false;

    // in case this is the first iteration
    if (voronoi_.weights().empty())
      voronoi_.weights().resize(particles_.n(), 0.0);

    // iterate until converged to the desired mass
    double error = 1;
    int iter = 0;
    for (iter = 0; iter < sim_opts.max_iter; iter++) {
      lift_sites(particles_, voronoi_.weights());
      voronoi_.compute(domain, voro_opts);
      compute_search_direction();

      // TODO backtracking to make sure cells always have positive volume
      error = length(gradient_);
      if (sim_opts.verbose)
        LOG << fmt::format("iter[{:3}]: error = {:.3e}", iter, error);
      if (error < sim_opts.volume_grad_tol)
        return {.error = error, .converged = true, .n_iterations = iter};

      for (size_t k = 0; k < voronoi_.weights().size(); k++)
        voronoi_.weights()[k] -= dw_[k];
    }
    return {.error = error, .converged = false, .n_iterations = iter};
  }

  void compute_search_direction();

  auto& particles() { return particles_; }
  const auto& voronoi() const { return voronoi_; }
  auto& voronoi() { return voronoi_; }

 private:
  Particles particles_;
  VoronoiDiagram voronoi_;
  spmat<double> hessian_;
  vecd<double> dw_;
  vecd<double> gradient_;
};

}  // namespace vortex