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
#include <unordered_set>
#include <vector>

#include "device_util.h"
#include "elements.h"
#include "voronoi.h"

namespace vortex {

template <typename Domain_t>
class VoronoiPolygon {
 public:
  using Cell_t = typename Domain_t::Cell_t;
  using Vertex_t = typename Cell_t::Vertex_t;
  typedef Polygon Element_t;

  // allocate space for the cell
  VoronoiPolygon(VoronoiCellMemoryPool& pool, size_t tid)
      : polygon_(pool.polygon.block(tid), pool.polygon.CAPACITY),
        plane_(pool.planes.block(tid), pool.planes.CAPACITY),
        bisector_to_site_(pool.bisector_to_site.block(tid),
                          pool.bisector_to_site.CAPACITY),
        neighbor_data_(512),
        distance_data_(512),
        site_stack_(128) {
    clear();
  }

  void set_sites(const coord_t* sites, size_t n_sites) {
    sites_ = sites;
    n_sites_ = n_sites;
  }
  void set_weights(const coord_t* weights) { weights_ = weights; }
  void set_neighbors(const index_t* neighbors, uint32_t n_neighbors) {
    neighbors_ = neighbors;
    n_neighbors_ = n_neighbors;
  }
  void set_max_neighbors(const uint16_t* max_neighbors) {
    max_neighbors_ = max_neighbors;
  }

  void set_kdtree(trees::KdTreeNd<coord_t, index_t>* tree) { tree_ = tree; }

  void clear() {
    polygon_.clear();
    plane_.clear();
    bisector_to_site_.clear(-1);
    max_radius_ = 0;
  }

  uint8_t new_plane(const vec4& p, int64_t n) {
    size_t b = plane_.size();
    plane_.push_back(p);
    bisector_to_site_.push_back(n);
    return b;
  }

  uint8_t compute_side(uint8_t bl, uint8_t br, const vec4& eqn) {
    return cell_.side(plane_[bl], plane_[br], eqn);
  }

  // calculate the Voronoi cell for site i clipped to the domain k_elem
  VoronoiStatusCode compute(const Domain_t& domain, int dim, size_t site,
                            VoronoiMesh& mesh) {
    assert(sites_);
    assert(neighbors_);
    const coord_t* zi = sites_ + site * dim;
    const coord_t wi = (weights_) ? weights_[site] : 0.0;
    const vec4 ui(zi, dim);

    // store initial neighbors
    size_t n_neighbors = max_neighbors_[site];
    neighbor_data_.resize(n_neighbors);
    for (size_t j = 0; j < n_neighbors; j++)
      neighbor_data_[j] = neighbors_[site * n_neighbors_ + j];

    bool security_radius_reached = false;
    VoronoiStatusCode status = VoronoiStatusCode::kIncomplete;
    for (int attempt = 1; !security_radius_reached; ++attempt) {
      // initialize the cell
      clear();
      domain.initialize(ui, *this);
      bisector_to_site_.set_size(plane_.size());
      for (size_t j = 1; j < n_neighbors; j++) {
        // get the next point and weight
        const index_t n = neighbor_data_[j];
        ASSERT(n < std::numeric_limits<index_t>::max())
            << "points may have duplicate coordinates, attempt = " << attempt
            << "# neighbors = " << n_neighbors << ", d = " << distance_data_[j];
        const coord_t* zj = sites_ + n * dim;
        const coord_t wj = (weights_) ? weights_[n] : 0.0;
        const vec4 uj(zj, dim);

        // create the plane and clip
        const vec4 eqn = plane_equation(ui, uj, wi, wj);
        const uint8_t b = new_plane(eqn, n);
        clip_by_plane(b);

        // check if no more bisectors contribute to the cell
        double sr = squared_radius(ui);
        max_radius_ = std::sqrt(sr);
        // security_radius_reached = 4.01 * sr < distance_squared(ui, uj);
        security_radius_reached =
            4.01 * sr < pow(Domain_t::length(ui.xyz(), uj.xyz()), 2) - wi + wj;
        if (security_radius_reached) break;
      }

      if (security_radius_reached) {
        status = VoronoiStatusCode::kSuccess;
        break;
      }
      if (!tree_ || attempt == kMaxClippingAttempts) {
        status = VoronoiStatusCode::kRadiusNotReached;
        break;
      }

      // get more neighbors
      n_neighbors *= 2;
      if (n_sites_ < n_neighbors) n_neighbors = n_sites_;
      if (n_neighbors < neighbor_data_.size()) {
        status = VoronoiStatusCode::kIncomplete;
        break;
      }
      neighbor_data_.resize(n_neighbors);
      distance_data_.resize(n_neighbors);
      trees::NearestNeighborSearch<index_t, coord_t> search(
          n_neighbors, neighbor_data_.data(), distance_data_.data());
      tree_->knearest(zi, search);
      ASSERT(neighbor_data_[0] == site);
    }

    // append to the mesh if necessary
    if (mesh.save_mesh()) append_to_mesh(mesh, site);
    if (mesh.save_facets()) save_facets(mesh, site);
    if (mesh.save_delaunay()) save_delaunay(mesh, site);

    return status;
  }

  VoronoiStatusCode compute(const TriangulationDomain& domain, int dim,
                            uint64_t triangle, uint64_t site,
                            VoronoiCellProperties* properties,
                            VoronoiMesh& mesh) {
    NOT_POSSIBLE;
    return VoronoiStatusCode::kSuccess;
  }

  void clip_by_plane(uint8_t b) {
    // retrieve the plane equation
    const vec4& eqn = plane_[b];

    const uint8_t b0 = polygon_[0];
    const uint8_t b1 = polygon_[1];
    int si = compute_side(b0, b1, eqn);

    uint8_t bi = b0;
    uint8_t bj = b1;
    uint8_t bk = polygon_[2];
    const size_t m = polygon_.size();
    size_t n = 0;
    for (size_t i = 0; i < m; i++) {
      const int sj = compute_side(bj, bk, eqn);

      if (si != sj) {
        // intersection
        if (si == INSIDE) {
          polygon_[n++] = bi;
          polygon_[n++] = bj;
        } else {
          polygon_[n++] = b;
        }
      } else if (si == INSIDE) {
        // both vertices are in this cell
        polygon_[n++] = bi;
      } else {
        // both vertices are outside the cell
        // no vertices get added
      }

      bi = bj;
      if (i + 3 == m) {
        bj = bk;
        bk = b0;
      } else if (i + 2 == m) {
        bj = b0;
        bk = b1;
      } else {
        bj = bk;
        bk = polygon_[i + 3];
      }

      si = sj;
    }
    polygon_.set_size(n);
  }

  double squared_radius(const vec4& c) const {
    if (polygon_.size() == 0) return 0.0;
    double r = 0.0;
    uint8_t bi = polygon_.back();
    for (size_t k = 0; k < polygon_.size(); k++) {
      vec4 p = cell_.compute(plane_[bi], plane_[polygon_[k]]);
      double rk = distance_squared(c, {p.x / p.w, p.y / p.w, p.z / p.w, 0});
      if (rk > r) r = rk;
      bi = polygon_[k];
    }
    return r;
  }

  void get_properties(VoronoiCellProperties& props, bool reset) const {
    if (reset) props.reset();
    cell_.get_properties(polygon_, plane_, props);
  }

  bool append_to_mesh(VoronoiMesh& mesh, index_t site) const {
    if (polygon_.size() < 3) return false;

    // initialize arrays for Voronoi polygon and Delaunay triangle
    std::vector<index_t> polygon(polygon_.size());
    std::array<index_t, 3> t;
    t[0] = site;

    index_t v_offset = mesh.vertices().n();
    for (size_t k = 0; k < polygon_.size(); k++) {
      uint8_t bi = polygon_[k];
      uint8_t bj = polygon_[(k + 1) % polygon_.size()];
      vec4 p = cell_.compute(plane_[bi], plane_[bj]);
      if (p.w == 0.0) p.w = 1.0;
      p = p / p.w;
      mesh.vertices().add(&p.x);
      polygon[k] = k + v_offset;

      // add Delaunay triangle associated with this Voronoi vertex
      auto tj = bisector_to_site_[bi];
      auto tk = bisector_to_site_[bj];
      if (tj < 0) tj = kMaxSite;
      if (tk < 0) tk = kMaxSite;
      t[1] = tj;
      t[2] = tk;
      mesh.triangles().add(t.data());
    }
    size_t id = mesh.polygons().n();
    mesh.polygons().add(polygon.data(), polygon.size());
    mesh.polygons().set_group(id, site);
    return true;
  }

  void save_facets(VoronoiMesh& mesh, index_t site_i) const {
    if (polygon_.size() < 3) return;

    // last point in the polygon
    size_t m = polygon_.size() - 1;
    vec4 p = cell_.compute(plane_[polygon_[m]], plane_[polygon_[0]]);
    if (p.w == 0.0) p.w = 1.0;
    p = p / p.w;
    for (size_t k = 0; k < polygon_.size(); k++) {
      uint8_t bi = polygon_[k];
      uint8_t bj = polygon_[(k + 1) % polygon_.size()];
      vec4 q = cell_.compute(plane_[bi], plane_[bj]);
      if (q.w == 0.0) q.w = 1.0;
      q = q / q.w;

      // TODO(philip) use the specialized cell_ to compute geometric quantities
      double l = Domain_t::length(p.xyz(), q.xyz());
      ASSERT(!std::isnan(l));
      // double l = length(p.xyz() - q.xyz());
      // double l = std::acos(dot(p.xyz(), q.xyz()));
      // vec3 c = 0.5 * (p.xyz() + q.xyz());
      vec3 c = Domain_t::project(0.5 * (p.xyz() + q.xyz()));

      // save the vertex for the next iteration
      p = q;

      // add Delaunay triangle associated with this Voronoi vertex
      auto site_j = bisector_to_site_[bi];
      m = k;
      // ASSERT(site_j >= 0);
      if (site_j < 0) {
        // mesh.n_boundary_facets()++;
        // mesh.boundary_area() += l;
        mesh.add_facet(site_i, -1, 0.0, c);
        continue;
      }
      mesh.add_facet(site_i, site_j, l, c);
    }
  }

  void save_delaunay(VoronoiMesh& mesh, index_t site) const {
    std::array<uint32_t, 3> t;
    t[0] = site;
    for (size_t k = 0; k < polygon_.size(); k++) {
      // add Delaunay triangle associated with this Voronoi vertex
      uint8_t bi = polygon_[k];
      uint8_t bj = polygon_[(k + 1) % polygon_.size()];
      auto tj = bisector_to_site_[bi];
      auto tk = bisector_to_site_[bj];
      if (tj < 0) continue;
      if (tk < 0) continue;
      t[1] = tj;
      t[2] = tk;
      mesh.add_triangle(t);
    }
  }

  Cell_t& cell() { return cell_; }
  auto& vertices() { return polygon_; }
  auto& planes() { return plane_; }
  const auto* sites() const { return sites_; }
  const auto& bisector_to_site(uint8_t b) { return bisector_to_site_[b]; }
  double max_radius() const { return max_radius_; }

 private:
  Cell_t cell_;
  device_vector<Vertex_t> polygon_;
  device_vector<vec4> plane_;
  device_vector<int64_t> bisector_to_site_;

  const index_t* neighbors_;
  const uint16_t* max_neighbors_;
  uint32_t n_neighbors_{0};
  const coord_t* sites_;
  size_t n_sites_{0};
  const coord_t* weights_{nullptr};
  device_vector<index_t> neighbor_data_;
  device_vector<coord_t> distance_data_;

  // for element-based clipping
  device_vector<index_t> site_stack_;
  device_hash_set<index_t> site_visited_;

  double max_radius_{0};

  // only CPU
  trees::KdTreeNd<coord_t, index_t>* tree_{nullptr};
};

template <>
VoronoiStatusCode VoronoiPolygon<TriangulationDomain>::compute(
    const TriangulationDomain& domain, int dim, uint64_t triangle,
    uint64_t site, VoronoiCellProperties* properties, VoronoiMesh& mesh) {
  assert(sites_);
  assert(neighbors_);

  site_stack_.clear();
  site_visited_.clear();

  // get the first site
  size_t i = site;
  site_stack_.push_back(i);
  site_visited_.insert(i);
  int iter = 0;
  while (!site_stack_.empty()) {
    iter++;
    clear();
    i = site_stack_.back();
    site_stack_.pop_back();
    const coord_t* zi = sites_ + i * dim;
    const vec4 ui(zi, dim);
    const coord_t wi = 0.0;

    size_t n_neighbors = n_neighbors_;
    neighbor_data_.resize(n_neighbors);
    for (int j = 0; j < n_neighbors; j++) {
      neighbor_data_[j] = neighbors_[i * n_neighbors_ + j];
    }

    bool security_radius_reached = false;
    for (int attempt = 1; attempt < kMaxClippingAttempts; ++attempt) {
      clear();
      domain.initialize(triangle, *this);
      bisector_to_site_.set_size(plane_.size());
      for (size_t j = 1; j < n_neighbors; j++) {
        // get the next point, TODO and weight
        const index_t n = neighbor_data_[j];
        if (!site_visited_.contains(n)) {
          site_stack_.push_back(n);
          site_visited_.insert(n);
        }
        const coord_t* zj = sites_ + n * dim;
        const vec4 uj(zj, dim);
        const coord_t wj = 0.0;

        // create the plane and clip
        const vec4 eqn = plane_equation(ui, uj, wi, wj);
        const uint8_t b = new_plane(eqn, n);
        clip_by_plane(b);

        // check if no more bisectors contribute to the cell
        double sr = squared_radius(ui);
        max_radius_ = std::sqrt(sr);
        security_radius_reached = 4.01 * sr < distance_squared(ui, uj);
        if (security_radius_reached) break;
      }
      if (security_radius_reached) break;
      // retrieve more neighbors and keep clipping
      if (!tree_) break;
      n_neighbors *= 2;
      if (n_sites_ < n_neighbors) n_neighbors = n_sites_;
      if (n_neighbors < neighbor_data_.size()) break;
      neighbor_data_.resize(n_neighbors);
      distance_data_.resize(n_neighbors);
      trees::NearestNeighborSearch<index_t, coord_t> search(
          n_neighbors, neighbor_data_.data(), distance_data_.data());
      tree_->knearest(zi, search);
    }
    // append to the mesh and properties if necessary
    if (mesh.save_mesh() && security_radius_reached) append_to_mesh(mesh, i);
    if (properties) get_properties(properties[i], false);
  }
  return VoronoiStatusCode::kSuccess;
}

}  // namespace vortex