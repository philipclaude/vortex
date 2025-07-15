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

    if (method_ == OperatorMethod::kKincl) {
      for (int d = 0; d < 3; d++) {
        grad_f[3 * i + d] += lij * (xi[d] - mij[d]) * (fi - fj) / (rij * wi);
        grad_f[3 * j + d] += lij * (xj[d] - mij[d]) * (fj - fi) / (rij * wj);
      }
    } else if (method_ == OperatorMethod::kSpringel) {
      vec3d cij = mij - 0.5 * (xi + xj);
      vec3d grad_ij =
          (lij / rij) * ((fj - fi) * cij - 0.5 * (fi + fj) * (xi - xj)) / wi;
      vec3d grad_ji =
          (lij / rij) * ((fi - fj) * cij - 0.5 * (fi + fj) * (xj - xi)) / wj;
      for (int d = 0; d < 3; d++) {
        grad_f[3 * i + d] += grad_ij[d];
        grad_f[3 * j + d] += grad_ji[d];
      }
    } else
      ASSERT(false) << "unknown method";
  }

  if (project_) {
    for (size_t i = 0; i < n_sites; i++) {
      vec3d gip;
      Domain_t::project(sites + dim * i, grad_f + 3 * i, &gip[0]);
      for (size_t j = 0; j < 3; j++) {
        grad_f[3 * i + j] = gip[j];
      }
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

    if (method_ == OperatorMethod::kKincl) {
      const vec3d uij = ui - uj;
      vec3d uip = uij, ujp = uij;
      if (project_) {
        // We want to project (xi - mij) to calculate uij . P(xi - mij).
        // This is = uij^T P(xi - mij) = (P(xi - mij))^T uij = (xi - mij)^TP^T
        // uij, which is = (xi - mij)^T P uij = (xi - mij) . (P uij) since the
        // projection P = I - nn^T is symmetric. So projecting uij works too.
        Domain_t::project(sites + dim * i, &uij[0], &uip[0]);
        Domain_t::project(sites + dim * j, &uij[0], &ujp[0]);
      }

      div_u[i] += lij * dot(xi - mij, uip) / (rij * wi);
      div_u[j] -= lij * dot(xj - mij, ujp) / (rij * wj);
    } else if (method_ == OperatorMethod::kSpringel) {
      vec3d cij = mij - 0.5 * (xi + xj);
      vec3d cijp = cij, cjip = cij;
      vec3d xij = xi - xj, xji = xj - xi;
      vec3d rijp = xij, rjip = xji;
      vec3d uij_avg = 0.5 * (ui + uj);
      if (project_) {
        Domain_t::project(sites + dim * i, &cij[0], &cijp[0]);
        Domain_t::project(sites + dim * j, &cij[0], &cjip[0]);
        Domain_t::project(sites + dim * i, &xij[0], &rijp[0]);
        Domain_t::project(sites + dim * j, &xji[0], &rjip[0]);
      }

      double div_i =
          (lij / rij) * (dot(cijp, uj - ui) - dot(uij_avg, rijp)) / wi;
      double div_j =
          (lij / rij) * (dot(cjip, ui - uj) - dot(uij_avg, rjip)) / wj;
      div_u[i] += div_i;
      div_u[j] += div_j;
    } else
      ASSERT(false) << "unknown method";
  }
}

template <typename Domain_t>
void VoronoiOperators<Domain_t>::calculate_relative_vorticity(const coord_t* u,
                                                              coord_t* w) {
  const size_t n_sites = voronoi_.n_sites();
  const coord_t* sites = voronoi_.sites();
  const int dim = voronoi_.dim();  // could be 3 or 4

  const auto& facets = voronoi_.facets();
  ASSERT(facets.size() > 0);

  // zero vorticity
  for (size_t i = 0; i < n_sites; i++) w[i] = 0.0;

  // add contribution to sites on both sides of each facet
  for (const auto& facet : facets) {
    const auto i = facet.bi;
    const auto j = facet.bj;
    if (j < 0 || size_t(j) >= n_sites) {
      for (int d = 0; d < 3; d++) w[i] = boundary_value_;
      continue;
    }
    const vec3d xi(sites + dim * i, 3);
    const vec3d xj(sites + dim * j, 3);
    double rij = length(xi - xj);
    // rij = Domain_t::length(xi, xj);
    const vec3d ui(u + 3 * i);
    const vec3d uj(u + 3 * j);
    const vec3d mij(&facet.midpoint.x);
    const vec3d dl = normalize(cross(0.5 * (xi + xj), xj - xi));
    const double lij = facet.length;
    const double wi = voronoi_.properties()[i].volume;
    const double wj = voronoi_.properties()[j].volume;

    ASSERT(wi > 0 || wj > 0) << fmt::format("wi = {}, wj = {}", wi, wj);
    ASSERT(rij > 0);
    ASSERT(!std::isnan(lij));

    const vec3d uij = 0.5 * (ui + uj);
    vec3d uijf = uij;
    if (project_) {
      Domain_t::project(&mij[0], &uij[0], &uijf[0]);
    }

    w[i] += dot(uijf, dl) * lij / wi;
    w[j] -= dot(uijf, dl) * lij / wj;
  }
}

template class VoronoiOperators<SquareDomain>;
template class VoronoiOperators<SphereDomain>;

}  // namespace vortex