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

template <typename Domain_t> class VoronoiPolygon {
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
  void set_neighbors(const index_t* neighbors, uint32_t n_neighbors) {
    neighbors_ = neighbors;
    n_neighbors_ = n_neighbors;
  }

  void set_kdtree(maple::KdTreeNd<coord_t, index_t>* tree) { tree_ = tree; }

  void clear() {
    p_.clear();
    q_.clear();
    plane_.clear();
  }

  uint8_t new_plane(const vec4& p) {
    size_t b = plane_.size();
    ASSERT(b < plane_.capacity());
    plane_.push_back(p);
    return b;
  }

  uint8_t compute_side(uint8_t k, const vec4& eqn) {
    return cell_.side(plane_[p_[k].bl], plane_[p_[k].br], eqn);
  }

  // calculate the Voronoi cell for site i clipped to the domain k_elem
  VoronoiStatusCode compute(const Domain_t& domain, int dim, size_t site,
                            Mesh* mesh = nullptr) {
    assert(sites_);
    assert(neighbors_);
    clear();
    const coord_t* zi = sites_ + site * dim;
    const vec4 ui(zi, dim);  // TODO use weights
    domain.initialize(ui, *this);
    bool security_radius_reached = false;

    for (size_t j = 1; j < n_neighbors_; j++) {
      // get the next point, TODO and weight
      const index_t n = neighbors_[site * n_neighbors_ + j];
      const coord_t* zj = sites_ + n * dim;
      const vec4 uj(zj, dim);

      // create the plane and clip
      const vec4 eqn = cell_.plane_equation(ui, uj);
      const uint8_t b = new_plane(eqn);
      clip_by_plane(b);

      // check if no more bisectors contribute to the cell
      double r = squared_radius(ui.xyz());
      security_radius_reached = 4.01 * r < distance_squared(ui.xyz(), uj.xyz());
      if (security_radius_reached) break;
    }
    if (!security_radius_reached) return VoronoiStatusCode::kRadiusNotReached;

    // append to the mesh if necessary
    if (mesh) append_to_mesh(*mesh, site);

    return VoronoiStatusCode::kSuccess;
  }

  VoronoiStatusCode compute(const TriangulationDomain& domain, int dim,
                            uint64_t triangle, uint64_t site,
                            ElementVoronoiWorkspace& workspace,
                            Mesh* mesh = nullptr) {
    NOT_POSSIBLE;
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

  double squared_radius(const vec3& c) const {
    double r = 0.0;
    for (size_t k = 0; k < p_.size(); k++) {
      const vec4 p = cell_.compute(plane_[p_[k].bl], plane_[p_[k].br]);
      double rk = distance_squared(c, {p.x / p.w, p.y / p.w, p.z / p.w});
      if (rk > r) r = rk;
    }
    return r;
  }

  void get_properties(VoronoiCellProperties& props) const {
    cell_.get_properties(p_, plane_, props);
  }

  bool append_to_mesh(Mesh& mesh, int site) const {
    if (p_.size() < 3) return false;
    index_t v_offset = mesh.vertices().n();

    std::vector<index_t> polygon(p_.size());
    for (size_t k = 0; k < p_.size(); k++) {
      vec4 p = cell_.compute(plane_[p_[k].bl], plane_[p_[k].br]);
      if (p.w == 0.0) p.w = 1.0;
      ASSERT(p.w != 0.0) << fmt::format("p = {}, {}, {}, {}", p.x, p.y, p.z,
                                        p.w);
      p = p / p.w;
      mesh.vertices().add(&p.x);
      polygon[k] = k + v_offset;
    }
    size_t id = mesh.polygons().n();
    mesh.polygons().add(polygon.data(), polygon.size());
    mesh.polygons().set_group(id, site);
    return true;
  }

  Cell_t& cell() { return cell_; }
  auto& vertices() { return p_; }
  auto& planes() { return plane_; }

 private:
  Cell_t cell_;
  pool<Vertex_t> p_;
  pool<Vertex_t> q_;
  pool<vec4> plane_;
  const index_t* neighbors_;
  uint32_t n_neighbors_{0};
  const coord_t* sites_;

  // only CPU
  maple::KdTreeNd<coord_t, index_t>* tree_{nullptr};
};

template <>
VoronoiStatusCode VoronoiPolygon<TriangulationDomain>::compute(
    const TriangulationDomain& domain, int dim, uint64_t triangle,
    uint64_t site, ElementVoronoiWorkspace& workspace, Mesh* mesh) {
  assert(sites_);
  assert(neighbors_);

  // get the first site
  size_t i = site;
  workspace.site_stack.push_back(i);
  int iter = 0;
  while (!workspace.site_stack.empty()) {
    iter++;
    clear();
    i = workspace.site_stack.back();
    workspace.site_stack.pop_back();
    workspace.site_visited.insert(i);
    const coord_t* zi = sites_ + i * dim;
    const vec4 ui(zi, dim);  // TODO use weights
    domain.initialize(triangle, *this);
    for (size_t j = 1; j < n_neighbors_; j++) {
      // get the next point, TODO and weight
      const index_t n = neighbors_[i * n_neighbors_ + j];
      if (workspace.site_visited.find(n) == workspace.site_visited.end())
        workspace.site_stack.push_back(n);
      const coord_t* zj = sites_ + n * dim;
      const vec4 uj(zj, dim);

      // create the plane and clip
      const vec4 eqn = cell_.plane_equation(ui, uj);
      const uint8_t b = new_plane(eqn);
      clip_by_plane(b);

      // check if no more bisectors contribute to the cell
      double r = squared_radius(ui.xyz());
      if (4.01 * r < distance_squared(ui.xyz(), uj.xyz())) break;
    }
    // append to the mesh if necessary
    if (mesh) append_to_mesh(*mesh, i);
  }
  return VoronoiStatusCode::kSuccess;
}

}  // namespace vortex