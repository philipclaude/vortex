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

  // initialize the vertex and plane pools
  VoronoiPolygon(Vertex_t* v, size_t nv, vec4* p, size_t np)
      : p_(v, nv / 2), q_(v + nv / 2, nv / 2), plane_(p, np) {
    assert(nv % 2 == 0);
  }

  void set_sites(const coord_t* sites) { sites_ = sites; }
  void set_weights(const coord_t* weights) { weights_ = weights; }
  void set_neighbors(const index_t* neighbors, uint32_t n_neighbors) {
    neighbors_ = neighbors;
    n_neighbors_ = n_neighbors;
  }

  void set_kdtree(trees::KdTreeNd<coord_t, index_t>* tree) { tree_ = tree; }

  void clear() {
    p_.clear();
    q_.clear();
    plane_.clear();
    std::fill(bisector_to_site_, bisector_to_site_ + 256, -1);
  }

  uint8_t new_plane(const vec4& p, int64_t n) {
    size_t b = plane_.size();
    ASSERT(b < plane_.capacity());
    plane_.push_back(p);
    bisector_to_site_[b] = n;
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
    clear();
    const coord_t* zi = sites_ + site * dim;
    const coord_t wi = (weights_) ? weights_[site] : 0.0;
    const vec4 ui(zi, dim);
    domain.initialize(ui, *this);
    bool security_radius_reached = false;

    for (size_t j = 1; j < n_neighbors_; j++) {
      // get the next point and weight
      const index_t n = neighbors_[site * n_neighbors_ + j];
      const coord_t* zj = sites_ + n * dim;
      const coord_t wj = (weights_) ? weights_[n] : 0.0;
      const vec4 uj(zj, dim);

      // create the plane and clip
      const vec4 eqn = plane_equation(ui, uj, wi, wj);
      const uint8_t b = new_plane(eqn, n);
      clip_by_plane(b);

      // check if no more bisectors contribute to the cell
      double r = squared_radius(ui);
      security_radius_reached = 4.01 * r < distance_squared(ui, uj);
      if (security_radius_reached) break;
    }

    // append to the mesh if necessary
    if (mesh.save_mesh()) append_to_mesh(mesh, site);
    if (mesh.save_facets()) save_facets(mesh, site);

    // TODO(philip) retrieve more neighbors and keep clipping
    if (!security_radius_reached) return VoronoiStatusCode::kRadiusNotReached;

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

    size_t nq = 0;
    q_.clear();
    int si = compute_side(0, eqn);
    for (size_t i = 0; i < p_.size(); i++) {
      int j = (i + 1 == p_.size()) ? 0 : i + 1;
      int sj = compute_side(j, eqn);

      if (si != sj) {
        // intersection
        if (si == INSIDE) {
          q_[nq].bl = p_[i].bl;
          q_[nq].br = p_[i].br;
          ++nq;
          q_[nq].bl = p_[i].br;
          q_[nq].br = b;
          ++nq;
        } else {
          q_[nq].bl = b;
          q_[nq].br = p_[j].bl;
          ++nq;
        }
      } else if (si == INSIDE) {
        // both vertices are in this cell
        q_[nq].bl = p_[i].bl;
        q_[nq].br = p_[i].br;
        ++nq;
      } else {
        // both vertices are outside the cell
        // no vertices get added
      }
      si = sj;
    }

    q_.set_size(nq);
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
    cell_.get_properties(p_, plane_, props);
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
      // double l = length(p.xyz() - q.xyz());
      double l = std::acos(dot(p.xyz(), q.xyz()));
      vec3 c = 0.5 * (p.xyz() + q.xyz());

      // save the vertex for the next iteration
      p = q;

      // add Delaunay triangle associated with this Voronoi vertex
      auto site_j = bisector_to_site_[p_[k].bl];
      // ASSERT(site_j >= 0);
      if (site_j < 0) continue;
      mesh.add(site_i, site_j, l, c);
    }
  }

  Cell_t& cell() { return cell_; }
  auto& vertices() { return p_; }
  auto& planes() { return plane_; }
  const auto* sites() const { return sites_; }
  const auto& bisector_to_site(uint8_t b) { return bisector_to_site_[b]; }

 private:
  Cell_t cell_;
  pool<Vertex_t> p_;
  pool<Vertex_t> q_;
  pool<vec4> plane_;
  const index_t* neighbors_;
  uint32_t n_neighbors_{0};
  const coord_t* sites_;
  const coord_t* weights_{nullptr};
  int64_t bisector_to_site_[256];  // 256 since bisectors are uint8_t

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
    domain.initialize(triangle, *this);
    for (size_t j = 1; j < n_neighbors_; j++) {
      // get the next point, TODO and weight
      const index_t n = neighbors_[i * n_neighbors_ + j];
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
      double r = squared_radius(ui);
      if (4.01 * r < distance_squared(ui, uj)) break;
    }
    // append to the mesh and properties if necessary
    if (mesh.save_mesh()) append_to_mesh(mesh, i);
    if (properties) get_properties(properties[i], false);
  }
  return VoronoiStatusCode::kSuccess;
}

}  // namespace vortex