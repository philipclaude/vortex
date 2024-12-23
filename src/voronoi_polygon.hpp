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
#include <unordered_set>
#include <vector>

#include "device_util.h"
#include "elements.h"
#include "voronoi.h"

namespace vortex {

struct ElementVoronoiWorkspace {
  std::vector<index_t> site_stack;
  std::unordered_set<index_t> site_visited;
  void clear() {
    site_stack.clear();
    site_visited.clear();
  }
};

template <typename Domain_t>
class VoronoiPolygon {
 public:
  using Cell_t = typename Domain_t::Cell_t;
  using Vertex_t = typename Cell_t::Vertex_t;
  typedef Polygon Element_t;

  // allocate space for the cell
  VoronoiPolygon(size_t capacity)
      : p_(capacity),
        q_(capacity),
        plane_(capacity),
        bisector_to_site_(capacity),
        neighbor_data_(capacity),
        distance_data_(capacity) {
    clear();
  }

  void set_sites(const coord_t* sites) { sites_ = sites; }
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
    p_.clear();
    q_.clear();
    plane_.clear();
    bisector_to_site_.clear(-1);
    max_radius_ = 0;
  }

  uint8_t new_plane(const vec4& p, int64_t n) {
    size_t b = plane_.size();
    plane_.push_back(p);  // TODO use emplace back
    bisector_to_site_.push_back(n);
    return b;
  }

  uint8_t compute_side(uint8_t k, const vec4& eqn) {
    return cell_.side(plane_[p_[k].bl], plane_[p_[k].br], eqn);
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
    for (int attempt = 0; !security_radius_reached; ++attempt) {
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
        security_radius_reached = 4.01 * sr < distance_squared(ui, uj);
        if (security_radius_reached) break;
      }

      // append to the mesh if necessary
      if (mesh.save_mesh()) append_to_mesh(mesh, site);
      if (mesh.save_facets()) save_facets(mesh, site);
      if (mesh.save_delaunay()) save_delaunay(mesh, site);

      if (security_radius_reached) break;
      if (!tree_) return VoronoiStatusCode::kRadiusNotReached;

      // get more neighbors
      n_neighbors *= 2;
      neighbor_data_.resize(n_neighbors);
      distance_data_.resize(n_neighbors);
      trees::NearestNeighborSearch<index_t, coord_t> search(
          n_neighbors, neighbor_data_.data(), distance_data_.data());
      tree_->knearest(zi, search);
      ASSERT(neighbor_data_[0] == site);
    }

    return VoronoiStatusCode::kSuccess;
  }

  VoronoiStatusCode compute(const TriangulationDomain& domain, int dim,
                            uint64_t triangle, uint64_t site,
                            ElementVoronoiWorkspace& workspace,
                            VoronoiCellProperties* properties,
                            VoronoiMesh& mesh) {
    NOT_POSSIBLE;
    return VoronoiStatusCode::kSuccess;
  }

  void clip_by_plane(uint8_t b) {
    // retrieve the plane equation
    const vec4& eqn = plane_[b];

    // size_t nq = 0;
    q_.clear();
    int si = compute_side(0, eqn);
    for (size_t i = 0; i < p_.size(); i++) {
      int j = (i + 1 == p_.size()) ? 0 : i + 1;
      int sj = compute_side(j, eqn);

      if (si != sj) {
        // intersection
        if (si == INSIDE) {
          // q_[nq].bl = p_[i].bl;
          // q_[nq].br = p_[i].br;
          auto& v0 = q_.emplace_back();
          v0.bl = p_[i].bl;
          v0.br = p_[i].br;
          //++nq;
          // q_[nq].bl = p_[i].br;
          // q_[nq].br = b;
          auto& v1 = q_.emplace_back();
          v1.bl = p_[i].br;
          v1.br = b;
          //++nq;
        } else {
          // q_[nq].bl = b;
          // q_[nq].br = p_[j].bl;
          auto& v = q_.emplace_back();
          v.bl = b;
          v.br = p_[j].bl;
          //++nq;
        }
      } else if (si == INSIDE) {
        // both vertices are in this cell
        // q_[nq].bl = p_[i].bl;
        // q_[nq].br = p_[i].br;
        auto& v = q_.emplace_back();
        v.bl = p_[i].bl;
        v.br = p_[i].br;
        //++nq;
      } else {
        // both vertices are outside the cell
        // no vertices get added
      }
      si = sj;
    }

    p_.swap(q_);
  }

  double squared_radius(const vec4& c) const {
    double r = 0.0;
    for (size_t k = 0; k < p_.size(); k++) {
      const vec4 p = cell_.compute(plane_[p_[k].bl], plane_[p_[k].br]);
      double rk = distance_squared(c, {p.x / p.w, p.y / p.w, p.z / p.w, 0});
      if (rk > r) r = rk;
    }
    return r;
  }

  void get_properties(VoronoiCellProperties& props, bool reset) const {
    if (reset) props.reset();
    for (size_t k = 0; k < p_.size(); k++) ASSERT(p_[k].bl != p_[k].br);
    cell_.get_properties(p_, plane_, props);
    props.rmax = max_radius_;
  }

  bool append_to_mesh(VoronoiMesh& mesh, index_t site) const {
    if (p_.size() < 3) return false;

    // initialize arrays for Voronoi polygon and Delaunay triangle
    std::vector<index_t> polygon(p_.size());
    std::array<index_t, 3> t;
    t[0] = site;

    index_t v_offset = mesh.vertices().n();
    for (size_t k = 0; k < p_.size(); k++) {
      vec4 p = cell_.compute(plane_[p_[k].bl], plane_[p_[k].br]);
      if (p.w == 0.0) p.w = 1.0;
      p = p / p.w;
      mesh.vertices().add(&p.x);
      polygon[k] = k + v_offset;

      // add Delaunay triangle associated with this Voronoi vertex
      auto tj = bisector_to_site_[p_[k].bl];
      auto tk = bisector_to_site_[p_[k].br];
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
    if (p_.size() < 3) return;

    // last point in the polygon
    size_t m = p_.size() - 1;
    vec4 p = cell_.compute(plane_[p_[m].bl], plane_[p_[m].br]);
    if (p.w == 0.0) p.w = 1.0;
    p = p / p.w;
    for (size_t k = 0; k < p_.size(); k++) {
      vec4 q = cell_.compute(plane_[p_[k].bl], plane_[p_[k].br]);
      if (q.w == 0.0) q.w = 1.0;
      q = q / q.w;

      // TODO(philip) use the specialized cell_ to compute geometric quantities
      double l = length(p.xyz() - q.xyz());
      // double l = std::acos(dot(p.xyz(), q.xyz()));

      // save the vertex for the next iteration
      p = q;

      // add Delaunay triangle associated with this Voronoi vertex
      auto site_j = bisector_to_site_[p_[k].bl];
      ASSERT(p_[k].bl == p_[m].br);
      m = k;
      // ASSERT(site_j >= 0);
      if (site_j < 0) {
        // mesh.n_boundary_facets()++;
        // mesh.boundary_area() += l;
        continue;
      }
      mesh.add(site_i, site_j, l);
    }
  }

  void save_delaunay(VoronoiMesh& mesh, index_t site) const {
    std::array<uint32_t, 3> t;
    t[0] = site;
    for (size_t k = 0; k < p_.size(); k++) {
      // add Delaunay triangle associated with this Voronoi vertex
      auto tj = bisector_to_site_[p_[k].bl];
      auto tk = bisector_to_site_[p_[k].br];
      if (tj < 0) continue;
      if (tk < 0) continue;
      t[1] = tj;
      t[2] = tk;
      mesh.add_triangle(t);
    }
  }

  Cell_t& cell() { return cell_; }
  auto& vertices() { return p_; }
  auto& planes() { return plane_; }
  const auto* sites() const { return sites_; }
  const auto& bisector_to_site(uint8_t b) { return bisector_to_site_[b]; }

 private:
  Cell_t cell_;
  device_vector<Vertex_t> p_;
  device_vector<Vertex_t> q_;
  device_vector<vec4> plane_;
  const index_t* neighbors_;
  const uint16_t* max_neighbors_;
  uint32_t n_neighbors_{0};
  const coord_t* sites_;
  const coord_t* weights_{nullptr};
  device_vector<int64_t> bisector_to_site_;
  device_vector<index_t> neighbor_data_;
  device_vector<coord_t> distance_data_;

  double max_radius_;

  // only CPU
  trees::KdTreeNd<coord_t, index_t>* tree_{nullptr};
};

template <>
VoronoiStatusCode VoronoiPolygon<TriangulationDomain>::compute(
    const TriangulationDomain& domain, int dim, uint64_t triangle,
    uint64_t site, ElementVoronoiWorkspace& workspace,
    VoronoiCellProperties* properties, VoronoiMesh& mesh) {
  assert(sites_);
  assert(neighbors_);

  // get the first site
  size_t i = site;
  workspace.site_stack.push_back(i);
  workspace.site_visited.insert(i);
  int iter = 0;
  while (!workspace.site_stack.empty()) {
    iter++;
    clear();
    i = workspace.site_stack.back();
    workspace.site_stack.pop_back();
    const coord_t* zi = sites_ + i * dim;
    const vec4 ui(zi, dim);
    const coord_t wi = 0.0;

    size_t n_neighbors = n_neighbors_;
    neighbor_data_.resize(n_neighbors);
    for (int j = 0; j < n_neighbors; j++) {
      neighbor_data_[j] = neighbors_[i * n_neighbors_ + j];
    }

    bool security_radius_reached = false;
    while (!security_radius_reached) {
      domain.initialize(triangle, *this);
      bisector_to_site_.set_size(plane_.size());
      for (size_t j = 1; j < n_neighbors; j++) {
        // get the next point, TODO and weight
        ASSERT(j < neighbor_data_.size()) << j;
        const index_t n = neighbor_data_[j];
        if (workspace.site_visited.find(n) == workspace.site_visited.end()) {
          workspace.site_stack.push_back(n);
          workspace.site_visited.insert(n);
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
        security_radius_reached = 4.01 * sr < distance_squared(ui, uj);
        if (security_radius_reached) break;
      }
      if (security_radius_reached) break;
      // retrieve more neighbors and keep clipping
      if (tree_) {
        n_neighbors *= 2;
        neighbor_data_.resize(n_neighbors);
        device_vector<coord_t> distance(n_neighbors);
        trees::NearestNeighborSearch<index_t, coord_t> search(
            n_neighbors, neighbor_data_.data(), distance.data());
        tree_->knearest(zi, search);
        clear();
      } else {
        break;
      }
    }
    // append to the mesh and properties if necessary
    if (mesh.save_mesh() && security_radius_reached) append_to_mesh(mesh, i);
    if (properties) get_properties(properties[i], false);
  }
  return VoronoiStatusCode::kSuccess;
}

}  // namespace vortex