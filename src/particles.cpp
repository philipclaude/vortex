#include "particles.h"

#include "math/vec.hpp"

namespace vortex {

void Particles::create(size_t np, const coord_t* xp, int dim) {
  Vertices::reserve(np);
  for (size_t k = 0; k < np; k++) Vertices::add(xp + dim * k);
  mass_.resize(np, 0.0);
}

void ParticleSimulation::compute_search_direction() {
  // calculate the gradient
  for (size_t k = 0; k < particles_.n(); k++)
    gradient_[k] = particles_.mass()[k] - voronoi_.properties()[k].mass;

  // set up the sparse matrix
  hessian_.clear();
  for (const auto& [b, f] : voronoi_.facets()) {
    size_t site_i = b[0];
    size_t site_j = b[1];
    vec3d pi(particles_[site_i]);
    vec3d pj(particles_[site_j]);
    double delta_ij = 0.5 * f.mass / length(pi - pj);
    hessian_(b[0], b[1]) = delta_ij;
    hessian_(b[1], b[0]) = delta_ij;
    hessian_(b[0], b[0]) -= delta_ij;
    hessian_(b[1], b[1]) -= delta_ij;
  }

  // solve for the search direction
  hessian_.solve_nl(gradient_, dw_, 1e-3, true);
}

}  // namespace vortex