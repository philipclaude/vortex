#include "operators.h"

#include "math/linalg.h"
#include "math/vec.hpp"
#include "voronoi.h"

namespace vortex {

#define USE_PROJECTION 0

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
  for (int i = 0; i < n_sites * 3; i++) grad_f[i] = 0.0;

  // add contribution to sites on both sides of each facet
  for (const auto& facet : facets) {
    const auto i = facet.bi;
    const auto j = facet.bj;
    if (j < 0 || j >= n_sites) {
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

    vec3d gi, gj;
    for (int d = 0; d < 3; d++) {
      gi[d] = lij * (xi[d] - mij[d]) * (fi - fj) / (rij * wi);
      gj[d] = lij * (xj[d] - mij[d]) * (fj - fi) / (rij * wj);
    }

    vec3d grad_i = gi, grad_j = gj;
#if USE_PROJECTION
    Domain_t::project(&xi[0], &gi[0], &grad_i[0]);
    Domain_t::project(&xj[0], &gj[0], &grad_j[0]);
#endif

    for (int d = 0; d < 3; d++) {
      grad_f[i * 3 + d] += grad_i[d];
      grad_f[j * 3 + d] += grad_j[d];
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
  for (int i = 0; i < n_sites; i++) div_u[i] = 0.0;

  // add contribution to sites on both sides of each facet
  for (const auto& facet : facets) {
    const auto i = facet.bi;
    const auto j = facet.bj;
    if (j < 0 || j >= n_sites) {
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

    vec3d uij = ui - uj;
    vec3d gi = uij, gj = uij;
#if USE_PROJECTION
    Domain_t::project(&xi[0], &uij[0], &gi[0]);
    Domain_t::project(&xj[0], &uij[0], &gj[0]);
#endif
    div_u[i] += lij * dot(xi - mij, gi) / (rij * wi);
    div_u[j] -= lij * dot(xj - mij, gj) / (rij * wj);
  }
}

template class VoronoiOperators<SquareDomain>;
template class VoronoiOperators<SphereDomain>;

}  // namespace vortex