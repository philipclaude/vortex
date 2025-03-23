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
#include "operators.h"

#include "math/linalg.h"
#include "math/vec.hpp"
#include "voronoi.h"

namespace vortex {

template <typename Domain_t>
VoronoiOperators<Domain_t>::VoronoiOperators(const VoronoiDiagram& voronoi)
    : voronoi_(voronoi) {}

template <typename Domain_t>
void VoronoiOperators<Domain_t>::calculate_gradient(const coord_t* f,
                                                    coord_t* grad_f) {
  const size_t n_sites = voronoi_.n_sites();
  const coord_t* sites = voronoi_.sites();
  const int dim = voronoi_.dim();  // could be 3 or 4

  const auto& facets = voronoi_.facets();
  ASSERT(facets.size() > 0);

  // zero gradient
  for (size_t i = 0; i < n_sites * 3; i++) grad_f[i] = 0.0;

  // add contribution to sites on both sides of each facet
  for (const auto& facet : facets) {
    const auto i = facet.bi;
    const auto j = facet.bj;
    if (j < 0 || size_t(j) >= n_sites) {
      for (int d = 0; d < 3; d++) grad_f[i * 3 + d] = boundary_value_;
      continue;
    }
    const vec3d xi(sites + dim * i, 3);
    const vec3d xj(sites + dim * j, 3);
    double rij = length(xi - xj);
    // rij = Domain_t::length(xi, xj);
    const double fi = f[i];
    const double fj = f[j];
    const vec3d mij(&facet.midpoint.x);
    const double lij = facet.length;
    const double wi = voronoi_.properties()[i].volume;
    const double wj = voronoi_.properties()[j].volume;

    ASSERT(wi > 0 || wj > 0) << fmt::format("wi = {}, wj = {}", wi, wj);
    ASSERT(rij > 0);
    ASSERT(!std::isnan(lij));
    for (int d = 0; d < 3; d++) {
#if 1
      grad_f[3 * i + d] += lij * (xi[d] - mij[d]) * (fi - fj) / (rij * wi);
      grad_f[3 * j + d] += lij * (xj[d] - mij[d]) * (fj - fi) / (rij * wj);
#else
      const coord_t cij_d = mij[d] - 0.5 * (xi[d] + xj[d]);
      grad_f[3 * i + d] += lij * cij_d * (fj - fi) / (rij * wi);
      grad_f[3 * j + d] += lij * cij_d * (fi - fj) / (rij * wj);
      grad_f[3 * i + d] -= 0.5 * lij * (fi + fj) * (xi[d] - xj[d]) / (rij * wi);
      grad_f[3 * j + d] -= 0.5 * lij * (fi + fj) * (xj[d] - xi[d]) / (rij * wj);
#endif
    }
  }
}

template <typename Domain_t>
void VoronoiOperators<Domain_t>::calculate_divergence(const coord_t* u,
                                                      coord_t* div_u) {
  const size_t n_sites = voronoi_.n_sites();
  const coord_t* sites = voronoi_.sites();
  const int dim = voronoi_.dim();  // could be 3 or 4

  const auto& facets = voronoi_.facets();
  ASSERT(facets.size() > 0);

  // zero divergence
  for (size_t i = 0; i < n_sites; i++) div_u[i] = 0.0;

  // add contribution to sites on both sides of each facet
  for (const auto& facet : facets) {
    const auto i = facet.bi;
    const auto j = facet.bj;
    if (j < 0 || size_t(j) >= n_sites) {
      for (int d = 0; d < 3; d++) div_u[i] = boundary_value_;
      continue;
    }
    const vec3d xi(sites + dim * i, 3);
    const vec3d xj(sites + dim * j, 3);
    double rij = length(xi - xj);
    // rij = Domain_t::length(xi, xj);
    const vec3d ui(u + 3 * i);
    const vec3d uj(u + 3 * j);
    const vec3d mij(&facet.midpoint.x);
    const double lij = facet.length;
    const double wi = voronoi_.properties()[i].volume;
    const double wj = voronoi_.properties()[j].volume;
    const vec3d uij = ui - uj;

    div_u[i] += lij * dot(xi - mij, uij) / (rij * wi);
    div_u[j] -= lij * dot(xj - mij, uij) / (rij * wj);
  }
}

template class VoronoiOperators<SquareDomain>;
template class VoronoiOperators<SphereDomain>;

}  // namespace vortex