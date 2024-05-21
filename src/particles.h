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
  void conserve_mass(const Domain_t& domain) {
    // calculate the voronoi diagram
    VoronoiDiagramOptions options;
    options.store_mesh = false;
    options.store_facet_data = true;
    options.allow_reattempt = false;
    options.n_neighbors = 100;
    options.verbose = false;

    // in case this is the first iteration
    if (voronoi_.weights().empty())
      voronoi_.weights().resize(particles_.n(), 0.0);

    // iterate until converged to the desired mass
    for (int iter = 0; iter < 100; iter++) {
      lift_sites(particles_, voronoi_.weights());
      voronoi_.compute(domain, options);
      compute_search_direction();

      // TODO backtracking to make sure cells always have positive volume
      double error = length(gradient_);
      LOG << fmt::format("iter. {}: error = {}", iter, error);
      if (error < 1e-8) {
        LOG << fmt::format("converged in {} iterations", iter);
        break;
      }

      for (size_t k = 0; k < voronoi_.weights().size(); k++)
        voronoi_.weights()[k] -= dw_[k];
    }
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