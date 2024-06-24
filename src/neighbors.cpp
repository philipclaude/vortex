#include "neighbors.h"

#include "voronoi.h"

namespace vortex {

VoronoiNeighbors::VoronoiNeighbors(const VoronoiDiagram& voronoi,
                                   const coord_t* points, int dim)
    : dim_(dim), voronoi_(voronoi), points_(points), ring_(-1) {}

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

  std::vector<index_t> first(voronoi_.n_sites(), 0);
  size_t m = count[0];
  for (size_t k = 1; k < voronoi_.n_sites(); k++) {
    first[k] = m;
    m += count[k];
  }

  // allocate the rings
  ring_.clear();
  ring_.set(first, count, m);

  // save ring data
  std::vector<uint32_t> idx(voronoi_.n_sites(), 0);
  for (const auto& [f, _] : facets) {
    uint32_t p = f.first;
    uint32_t q = f.second;
    ring_(p, idx[p]++) = q;
    ring_(q, idx[q]++) = p;
  }
}

double distance_squared(const double* p, const double* q, int dim) {
  double ds = 0;
  for (int i = 0; i < dim; ++i) {
    const double dx = p[i] - q[i];
    ds += dx * dx;
  }
  return ds;
}

void VoronoiNeighbors::knearest(uint32_t p,
                                NearestNeighborsWorkspace& search) const {
  search.reset();
  search.add(p, 0, 0);
  while (!search.sites.empty()) {
    uint32_t site = search.next();
    uint8_t level = search.neighbors.at(site);
    if (level >= search.max_level) continue;

    for (int j = 0; j < ring_.length(site); j++) {
      auto n = ring_[site][j];
      if (n == p) continue;
      if (search.neighbors.find(n) != search.neighbors.end()) continue;
      double d = distance_squared(&points_[dim_ * p], &points_[dim_ * n], dim_);
      if (level > 1 && d > search.max_distance) continue;
      search.add(n, level + 1, d);
    }
  }
  search.sort();
  search.total_neighbors += search.sites.size();
}

}  // namespace vortex