//
//  vortex: Voronoi mesher and fluid simulator for the Earth's oceans and
//  atmosphere.
//
//  Copyright 2023 - 2025 Philip Claude Caplan
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
#include "particles.h"

#include "io.h"
#include "math/mat.hpp"
#include "math/vec.hpp"

namespace vortex {

template <>
bool has_boundary<SquareDomain>() {
  return true;
}

template <>
bool has_boundary<SphereDomain>() {
  return false;
}

template <>
void project_point<SquareDomain>(double* x) {
  x[2] = 0;
}

template <>
void project_point<SphereDomain>(double* x) {
  vec3d p(x);
  p = normalize(p);
  for (int d = 0; d < 3; d++) x[d] = p[d];
}

template <>
void project_velocity<SquareDomain>(double* x, double* v) {}

template <>
void project_velocity<SphereDomain>(double* x, double* v) {
  vec3d n(x);
  vec3d u(v);
  n = normalize(n);
  vec3d ur = dot(u, n) * n;  // radial component of velocity
  vec3d ut = u - ur;         // tangential component of velocity
  for (int d = 0; d < 3; d++) v[d] = ut[d];
}

void Particles::create(size_t np, const coord_t* xp, int dim) {
  Vertices::reserve(np);
  vec3d uvw;
  centroids_.reserve(np);
  velocity_.reserve(np);
  for (size_t k = 0; k < np; k++) {
    vec4d xk(xp + dim * k, dim);
    Vertices::add(&xk[0]);
    centroids_.add(&xk[0]);
    velocity_.add(&uvw[0]);
  }
  mass_.resize(np, 1);
  volume_.resize(np, 1);
}

void Particles::save(const std::string& filename) const {
  size_t n_data = n() * 3;
  std::vector<float> data(n_data);
  size_t m = 0;
  for (size_t k = 0; k < n(); k++)
    for (int d = 0; d < 3; d++) data[m++] = (*this)[k][d];
  FILE* fid = fopen(filename.c_str(), "wb");
  fprintf(fid, "# vtk DataFile Version 2.0\nvortex vertices\n");
  fprintf(fid, "BINARY\nDATASET UNSTRUCTURED_GRID\nPOINTS %zu float\n", n());
  std::parafor_i(0, n_data,
                 [&data](int tid, size_t k) { io::swap_end(data[k]); });
  fwrite(&data[0], 4, n_data, fid);
  fprintf(fid, "\nCELLS 0 0\nCELL_TYPES 0\n");
  fprintf(fid, "POINT_DATA %zu\n", n());
  fprintf(fid, "FIELD FieldData 1\n");
  fprintf(fid, "density 1 %zu float\n", n());
  std::vector<float> density_data(n());
  for (size_t k = 0; k < n(); k++) density_data[k] = density_[k] > 1 ? 1000 : 0;
  std::parafor_i(0, n(), [&density_data](int tid, size_t k) {
    io::swap_end(density_data[k]);
  });
  fwrite(&density_data[0], 4, n(), fid);
  fprintf(fid, "\n");
  fclose(fid);
}

void ParticleSimulation::compute_search_direction(
    const std::vector<double>& target_volumes) {
  hessian_.clear();

  // add regularization?
  // for (size_t k = 0; k < particles_.n(); k++)
  //  hessian_(k, k) = 1e-3 * target_volumes[k];

  // set up the sparse matrix
  ASSERT(voronoi_.facets().size() > 0);
  for (const auto& facet : voronoi_.facets()) {
    if (facet.bj < 0) continue;
    size_t site_i = facet.bi;
    size_t site_j = facet.bj;
    if (site_i >= particles_.n()) continue;
    if (site_j >= particles_.n()) continue;
    vec3d pi(particles_[site_i]);
    vec3d pj(particles_[site_j]);
    double delta_ij = 0.5 * facet.length / length(pi - pj);
    hessian_(site_i, site_j) = delta_ij;
    hessian_(site_j, site_i) = delta_ij;
    hessian_(site_i, site_i) -= delta_ij;
    hessian_(site_j, site_j) -= delta_ij;
  }

  // solve for the search direction
  // a tolerance of 1e-3 should give a good enough direction
  SparseSolverOptions opts;
  opts.tol = 1e-3;
  opts.symmetric = true;
  opts.max_iterations = 100;
  hessian_.solve_nl(gradient_, dw_, opts);
}

void ParticleSimulation::calculate_properties() {
  for (size_t k = 0; k < particles_.n(); k++) {
    ASSERT(voronoi_.properties()[k].site == k);
    particles_.volume()[k] = voronoi_.properties()[k].volume;
    ASSERT(particles_.volume()[k] > 0);
    for (int d = 0; d < 3; d++)
      particles_.centroids()(k, d) =
          voronoi_.properties()[k].moment[d] / particles_.volume()[k];
    max_displacement_[k] = std::numeric_limits<double>::max();
  }

  // determine the maximum displacement as the min (bisector distance)/2
  for (const auto& facet : voronoi_.facets()) {
    if (facet.bj < 0) continue;
    size_t site_i = facet.bi;
    size_t site_j = facet.bj;
    if (site_i >= particles_.n()) continue;
    if (site_j >= particles_.n()) continue;
    ASSERT(site_i < particles_.n());
    ASSERT(site_j < particles_.n());
    vec3d pi(particles_[site_i]);
    vec3d pj(particles_[site_j]);
    double wi = voronoi_.weights()[site_i];
    double wj = voronoi_.weights()[site_j];
    coord_t l = length(pj - pi);
    coord_t d1 = (l * l + wi - wj) / (2.0 * l);
    vec3d m = pi + d1 * (pj - pi) / l;

    double d = 0.5 * length(pi - m);
    if (d < max_displacement_[site_i]) max_displacement_[site_i] = d;
    if (d < max_displacement_[site_j]) max_displacement_[site_j] = d;
  }
}

}  // namespace vortex