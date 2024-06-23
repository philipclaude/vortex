#include "neighbors.h"

#include "voronoi.h"

namespace vortex {

VoronoiNeighbors::VoronoiNeighbors(const VoronoiDiagram& voronoi,
                                   const coord_t* vertices)
    : voronoi_(voronoi), points_(vertices), ring_(-1) {
  build();
}

void VoronoiNeighbors::build() {
  // count how many vertices are in each vertex one-ring
  const auto& facets = voronoi_.facets();
  std::vector<uint32_t> count(voronoi_.n_sites(), 0);
  for (const auto& [f, _] : facets) {
    uint32_t p = f.first;
    uint32_t q = f.second;
    count[p]++;
    count[q]++;
  }

  // allocate the rings
  ring_.reserve(voronoi_.n_sites(), -2 * facets.size());
  std::vector<uint32_t> ring;
  for (size_t k = 0; k < voronoi_.n_sites(); k++) {
    ring.resize(count[k]);
    ring_.add(ring.data(), ring.size());
  }

  // save ring data
  std::vector<uint32_t> idx(voronoi_.n_sites(), 0);
  for (const auto& [f, _] : facets) {
    uint32_t p = f.first;
    uint32_t q = f.second;
    ring_(p, idx[p]++) = q;
    ring_(q, idx[q]++) = p;
  }
}

template <int dim>
double distance_squared(const double* p, const double* q) {
  double ds = 0;
  for (int i = 0; i < dim; ++i) ds += (p[i] - q[i]) * (p[i] - q[i]);
  return ds;
}

void VoronoiNeighbors::knearest(uint32_t p,
                                NearestNeighborsWorkspace& search) const {
  const int max_level = 2;  // TODO input option
  static constexpr int dim = 3;

  search.reset();
  search.add(p, 0, 0);
  while (!search.sites.empty()) {
    auto site = search.sites.front().first;
    search.sites.pop();
    uint8_t level = search.neighbors.at(site);
    if (level >= max_level) continue;

    for (int j = 0; j < ring_.length(site); j++) {
      auto n = ring_[site][j];
      if (n == p) continue;
      if (search.neighbors.find(n) != search.neighbors.end()) continue;
      double d = distance_squared<dim>(&points_[dim * p], &points_[dim * n]);
      if (level > 1 && d > search.max_distance) continue;
      search.add(n, level + 1, d);
    }
  }
  search.sort();
  search.n_avg += search.sites.size();
}

}  // namespace vortex