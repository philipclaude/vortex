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
#include "neighbors.h"

#include <stlext.h>

#include "library.h"
#include "math/vec.hpp"
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
      if (search.limit == NeighborSearchLimit::kDistanceBased && level > 1 &&
          d > search.max_distance)
        continue;
      if (search.limit == NeighborSearchLimit::kCapacityBased &&
          search.neighbors.size() > search.n_neighbors)
        continue;
      search.add(n, level + 1, d);
    }
  }
  search.sort();
  search.total_neighbors += search.sites.size();
}

SphereQuadtree::SphereQuadtree(const coord_t* points, size_t n_points, int dim,
                               int ns)
    : points_(points), n_points_(n_points), dim_(dim), mesh_(n_points, ns) {
  setup();
}

SphereQuadtree::Subdivision::Subdivision(int np, int ns)
    : Mesh(3), children(4), n_levels(ns + 1) {
  using T = Octahedron;
  if (ns < 0) {
    int k = MAX_NEIGHBOR_CAPACITY;
    int n = np;
    int m = double(k) / kOneRingSize;
    ns = std::floor(std::log(n / (T ::n_faces * m)) / std::log(4.0));
  }
  n_levels = ns + 1;
  n_points_per_triangle = np / (Octahedron::n_faces * std::pow(4, ns));
  LOG << fmt::format("# levels = {}, # triangles = {}, # pts/tri ~ {}",
                     n_levels, t_last(ns) - t_first(ns), n_points_per_triangle);

  // initial octahedron
  for (int i = 0; i < T::n_vertices; i++) vertices_.add(T::coordinates[i]);
  for (int i = 0; i < T::n_faces; i++) {
    triangles_.add(T::faces[i]);
    triangles_.set_group(i, 0);
  }

  triangles_.reserve(T::n_faces * std::pow(4, ns));
  children.reserve(T::n_faces * std::pow(4, ns));
  std::unordered_map<Edge, index_t> edges;
  int n_triangles0 = 0;
  for (int j = 0; j < ns; j++) {
    // perform the subdivision
    edges.clear();
    int n_triangles = triangles_.n();
    for (size_t k = n_triangles0; k < n_triangles; k++) {
      int edge_indices[3];
      for (int j = 0; j < 3; j++) {
        index_t e0 = triangles_(k, j);
        index_t e1 = triangles_(k, j == 2 ? 0 : j + 1);

        index_t idx;
        auto it = edges.find({e1, e0});
        if (it == edges.end()) {
          // create a new point on this edge
          vec3d p0(vertices_[e0]);
          vec3d p1(vertices_[e1]);
          vec3d q = normalize(0.5 * (p0 + p1));

          idx = vertices_.n();
          vertices_.add(q.data());
          edges.insert({{e0, e1}, idx});
        } else {
          idx = it->second;
          edges.erase(it);
        }
        edge_indices[j] = idx;
      }

      // create the four new triangles from the subdivision
      index_t t0 = triangles_(k, 0);
      index_t t1 = triangles_(k, 1);
      index_t t2 = triangles_(k, 2);
      index_t e0 = edge_indices[0];
      index_t e1 = edge_indices[1];
      index_t e2 = edge_indices[2];

      index_t new_triangles[4][3] = {
          {t0, e0, e2}, {e0, t1, e1}, {e2, e1, t2}, {e0, e1, e2}};
      int c[4];
      size_t nt = triangles_.n();
      for (int i = 0; i < 4; i++) {
        triangles_.add(new_triangles[i]);
        triangles_.set_group(nt + i, j + 1);  // subdivision level
        c[i] = nt + i;
      }
      children.add(c);  // save the children of the current triangle
    }
    n_triangles0 = n_triangles;
  }
}

void SphereQuadtree::setup() {
  int last_level = mesh_.n_levels - 1;
  size_t t0 = mesh_.t_first(last_level);
  size_t t1 = mesh_.t_last(last_level);
  size_t n_triangles = t1 - t0;
  ASSERT(t1 == mesh_.triangles().n()) << t1 << ", " << mesh_.triangles().n();

  std::vector<std::vector<uint32_t>> v2t(mesh_.vertices().n());
  for (size_t k = t0; k < t1; k++) {
    for (int j = 0; j < 3; j++) {
      auto v = mesh_.triangles()(k, j);
      v2t[v].push_back(k);
    }
  }

  std::unordered_set<uint32_t> stri;
  std::array<int32_t, kOneRingSize> vtri;
  one_ring_tris_.resize(n_triangles);
  for (size_t k = t0; k < t1; k++) {
    stri.clear();
    for (int j = 0; j < 3; j++) {
      auto v = mesh_.triangles()[k][j];
      for (auto t : v2t[v]) stri.insert(t);
    }

    std::fill(vtri.begin(), vtri.end(), -1);
    int i = 0;
    for (auto t : stri) vtri[i++] = t - t0;
    one_ring_tris_[k - t0] = vtri;
  }

  if (mesh_.n_points_per_triangle < 10) {
    LOG << "activating two-ring";
    use_two_ring_ = true;
    struct TwoRingWorkspace {
      std::unordered_set<uint32_t> stri;
      std::array<int32_t, kTwoRingSize> vtri;
    };
    std::vector<TwoRingWorkspace> ws(std::thread::hardware_concurrency());
    two_ring_tris_.resize(n_triangles);
    std::parafor_i(0, t1 - t0, [&](int tid, size_t k) {
      auto& stri2 = ws[tid].stri;
      stri2.clear();
      const auto& tris = one_ring_tris_[k];
      for (auto t : tris) {
        for (int j = 0; j < 3; j++) {
          auto v = mesh_.triangles()[t + t0][j];
          for (auto tt : v2t[v]) stri2.insert(tt);
        }
      }

      std::fill(ws[tid].vtri.begin(), ws[tid].vtri.end(), -1);
      int i = 0;
      for (auto t : stri2) ws[tid].vtri[i++] = t - t0;
      two_ring_tris_[k] = ws[tid].vtri;
    });
  }
}

void SphereQuadtree::build() {
  // utility to determine if a point is inside a spherical triangle
  auto intriangle = [](const vec3d& v0, const vec3d& v1, const vec3d& v2,
                       const vec3d& v) {
    return dot(cross(v0, v1), v) >= 0 && dot(cross(v1, v2), v) >= 0 &&
           dot(cross(v2, v0), v) >= 0;
  };

  int last_level = mesh_.n_levels - 1;
  int t0 = mesh_.t_first(last_level);
  int t1 = mesh_.t_last(last_level);
  size_t n_triangles = t1 - t0;

  // build the point2triangle info
  point2triangle_.resize(n_points_);
  std::parafor_i(0, n_points_, [&](int tid, size_t k) {
    // first determine if we are in the first four or last four triangles
    int c[4] = {0, 1, 2, 3};
    vec3d p(points_ + dim_ * k);
    if (p[2] < 0) {
      for (int i = 0; i < 4; i++) c[i] += 4;
    }

    int level = -1;
    while (true) {
      int child = -1;
      for (int i = 0; i < 4; i++) {
        const auto* t = mesh_.triangles()[c[i]];
        vec3d v0(mesh_.vertices()[t[0]]);
        vec3d v1(mesh_.vertices()[t[1]]);
        vec3d v2(mesh_.vertices()[t[2]]);
        if (intriangle(v0, v1, v2, p)) {
          child = c[i];
          break;
        }
      }
      ASSERT(child >= 0);
      ASSERT(mesh_.triangles().group(child) == level + 1);

      level = mesh_.triangles().group(child);
      if (level == last_level) {
        point2triangle_[k] = child - t0;
        break;
      }

      // save the triangles for the next level
      for (int i = 0; i < 4; i++) c[i] = mesh_.children[child][i];
    }
  });

  // convert to triangle2point
  triangle2points_.clear();
  triangle2points_.resize(n_triangles);
  for (size_t k = 0; k < point2triangle_.size(); k++) {
    triangle2points_[point2triangle_[k]].push_back(k);
  }

  // gather statistics
  min_leaf_size_ = std::numeric_limits<int>::max();
  max_leaf_size_ = 0;
  for (const auto& pts : triangle2points_) {
    size_t n = pts.size();
    if (n < min_leaf_size_) min_leaf_size_ = n;
    if (n > max_leaf_size_) max_leaf_size_ = n;
  }
}

void SphereQuadtree::knearest(uint32_t p,
                              SphereQuadtreeWorkspace& search) const {
  search.reset();

  // get the triangle containing this point
  auto t = point2triangle_[p];

  // loop through all the triangles in the first layer around this triangle
  const auto* triangles =
      use_two_ring_ ? two_ring_tris_[t].data() : one_ring_tris_[t].data();
  size_t nt = use_two_ring_ ? kTwoRingSize : kOneRingSize;
  for (size_t k = 0; k < nt; k++) {
    // loop through all the points in this triangle
    if (triangles[k] < 0) break;
    const auto& points = triangle2points_[triangles[k]];
    for (size_t j = 0; j < points.size(); j++) {
      double d = distance_squared(points_ + dim_ * p,
                                  points_ + dim_ * points[j], dim_);
      search.add(points[j], d);
    }
  }
  search.sort();
}

}  // namespace vortex