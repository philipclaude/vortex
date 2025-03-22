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
#include "voronoi.h"

#include <Predicates_psm.h>

#include <memory>
#include <nlohmann/json.hpp>
#include <set>
#include <thread>
#include <unordered_set>

#include "elements.h"
#include "log.h"
#include "mesh.h"
#include "neighbors.h"
#include "stlext.h"
#include "trees/kdtree.h"
#include "voronoi_polygon.hpp"

namespace vortex {

void SquareDomain::initialize(vec4,
                              VoronoiPolygon<SquareDomain>& polygon) const {
  auto& vertices = polygon.vertices();
  auto& planes = polygon.planes();
  polygon.cell().initialize(&points_[0], 4, vertices, planes);
}

void SphereDomain::initialize(vec4 site,
                              VoronoiPolygon<SphereDomain>& polygon) const {
  auto& vertices = polygon.vertices();
  auto& planes = polygon.planes();
  vec3 center = unit_vector(site.xyz());
  const double r = initialization_fraction * radius;
  polygon.cell().center = center;

  // compute the normal and tangent to the sphere
  vec3 n = unit_vector(center);
  vec3 t = unit_vector(vec3{-n.y, n.x, 0.0});
  if (n.x * n.x + n.y * n.y == 0.0) {
    t = {1, 0, 0};
  }
  vec3 u = cross(t, n);

  vec3 p0 = center - r * u;
  vec3 p1 = center - r * t;
  vec3 p2 = center + r * u;
  vec3 p3 = center + r * t;

  // the equations of the four planes of the bounding square
  planes.resize(4);
  planes[0] = {u.x, u.y, u.z, -dot(p0, u)};
  planes[1] = {t.x, t.y, t.z, -dot(p1, t)};
  planes[2] = {-u.x, -u.y, -u.z, dot(p2, u)};
  planes[3] = {-t.x, -t.y, -t.z, dot(p3, t)};

  // create the 4 vertices for the bounding square
  vertices.resize(4);
  vertices[0] = 0;
  vertices[1] = 1;
  vertices[2] = 2;
  vertices[3] = 3;
}

uint8_t SphericalVoronoiPolygon::side(const vec4& pi, const vec4& pj,
                                      const vec4& p) const {
  return plane_side(compute(pi, pj), p);
}

vec4 SphericalVoronoiPolygon::compute(const vec4& pi, const vec4& pj) const {
  // 1. compute the line of intersection of two planes: x(t) = p + r * t
  // robust version described here:
  // https://stackoverflow.com/questions/6408670/line-of-intersection-between-two-planes/18092154#18092154
  vec3 n1 = pi.xyz();
  vec3 n2 = pj.xyz();
  coord_t h1 = pi.w;
  coord_t h2 = pj.w;

  vec3 r = cross(n1, n2);
  coord_t det = dot(r, r);
  ASSERT(det != 0.0) << "planes are parallel";
  vec3 p = (1.0 / det) * (h1 * cross(r, n2) + h2 * cross(n1, r));
  r = (1.0 / sqrt(det)) * r;

  // 2. compute the intersection between the line and a sphere
  coord_t b = dot(r, p);
  coord_t c = dot(p, p) - 1.0;

  coord_t discriminant = b * b - c;
  if (discriminant < 0) return {1e10, 1e10, 1e10, 1};
  if (std::isnan(discriminant)) {
    LOG << fmt::format("n1 = {}, {}, {}", n1.x, n1.y, n1.z);
    LOG << fmt::format("n2 = {}, {}, {}", n2.x, n2.y, n2.z);
    LOG << fmt::format("h1 = {}, h2 = {}", h1, h2);
  }
  ASSERT(discriminant >= 0.0) << fmt::format("sqrt({})", discriminant);

  coord_t t1 = -b + sqrt(discriminant);
  coord_t t2 = -b - sqrt(discriminant);

  vec3 pa = p + t1 * r;
  vec3 pb = p + t2 * r;

  coord_t distance1 = distance_squared(pa, center);
  coord_t distance2 = distance_squared(pb, center);

  // 3. keep the one that is closest to the center (site)
  if (distance1 < distance2) {
    return {pa.x, pa.y, pa.z, 1.0};
  }
  return {pb.x, pb.y, pb.z, 1.0};
}

void SphericalVoronoiPolygon::get_properties(
    const device_vector<Vertex_t>& p, const device_vector<vec4>& planes,
    VoronoiCellProperties& props) const {
  if (p.size() < 3) return;

  // https://math.stackexchange.com/questions/1143354/numerically-stable-method-for-angle-between-3d-vectors
  auto get_angle = [](const vec3& u, const vec3& v) {
    return 2 * atan2(length(u - v), length(u + v));
  };

  vec4 ah = compute(planes[p[0]], planes[p[1]]);
  vec4 bh = compute(planes[p[1]], planes[p[2]]);
  vec3 a = (1.0 / ah.w) * ah.xyz();
  vec3 b = (1.0 / bh.w) * bh.xyz();

  // https://stackoverflow.com/questions/19897187/locating-the-centroid-center-of-mass-of-spherical-polygons#answer-38201499
  coord_t tab = get_angle(a, b);
  props.moment = props.moment + 0.5 * tab * unit_vector(cross(a, b));

  for (size_t k = 2; k < p.size(); k++) {
    vec4 ch = compute(planes[p[k]], planes[p[(k + 1) % p.size()]]);
    vec3 c = (1.0 / ch.w) * ch.xyz();

    // https://www.johndcook.com/blog/2021/11/29/area-of-spherical-triangle/
    coord_t num = std::fabs(dot(a, cross(b, c)));
    coord_t den = 1.0 + dot(a, b) + dot(b, c) + dot(a, c);
    coord_t ak = 2.0 * std::atan2(num, den);
    ASSERT(ak == ak);

    const coord_t tbc = get_angle(b, c);
    props.moment = props.moment + 0.5 * tbc * unit_vector(cross(b, c));
    props.volume += ak;
    b = c;
  }

  tab = get_angle(b, a);
  props.moment = props.moment + 0.5 * tab * unit_vector(cross(b, a));
}

void PlanarVoronoiPolygon::initialize(const vec3* points, const size_t n_points,
                                      device_vector<Vertex_t>& vertices,
                                      device_vector<vec4>& planes) {
  vertices.resize(n_points);
  planes.resize(n_points);
  const vec3 u = points[1] - points[0];
  const vec3 v = points[n_points - 1] - points[0];
  const vec3 n = unit_vector(cross(u, v));
  base = {n.x, n.y, n.z, -dot(n, points[0])};

  for (size_t i = 0; i < n_points; ++i) {
    const uint16_t j = (i + 1) % n_points;
    vertices[i] = i;
    const vec3 w = unit_vector(cross(n, points[j] - points[i]));
    planes[i] = {w.x, w.y, w.z, -dot(w, points[i])};
  }
}

uint8_t PlanarVoronoiPolygon::side(const vec4& pi, const vec4& pj,
                                   const vec4& p) const {
  // The implementation for this function was inspired by voronoi.cu in the
  // Supplemental Material for "Meshless voronoi on the GPU":
  // (https://dl.acm.org/doi/10.1145/3272127.3275092)
  const vec4& pk = base;
  bool exact = true;
  double det = det4x4(pi.x, pj.x, pk.z, p.x, pi.y, pj.y, pk.y, p.y, pi.z, pj.z,
                      pk.z, p.z, pi.w, pj.w, pk.w, p.w);
  double maxx = std::max({fabs(pi.x), fabs(pj.x), fabs(pk.x), fabs(p.x)});
  double maxy = std::max({fabs(pi.y), fabs(pj.y), fabs(pk.y), fabs(p.y)});
  double maxz = std::max({fabs(pi.z), fabs(pj.z), fabs(pk.z), fabs(p.z)});
  double eps = 1.2466136531027298e-13 * maxx * maxy * maxz;
  // double min_max = std::min({maxx, maxy, maxz});
  double max_max = std::max({maxx, maxy, maxz});

  eps *= (max_max * max_max);
  if (fabs(det) < eps) exact = true;

  if (exact) {
    // see Partial Optimal Transport for a Constant-Volume Lagrangian Mesh
    // with Free Boundaries [Levy, 2021], page 19 ("The first phase")
    GEO::Sign det1 = GEO::PCK::det_3d(&pi.x, &pj.x, &pk.x);
    GEO::Sign det2 = GEO::PCK::det_4d(&pi.x, &pj.x, &pk.x, &p.x);
    return (det1 * det2 <= 0) ? OUTSIDE : INSIDE;
  }
  return plane_side(compute(pi, pj), p);
}

vec4 PlanarVoronoiPolygon::compute(const vec4& pi, const vec4& pj) const {
  // The implementation of this function was inspired by the
  // "ConvexCell::compute_triangle_point" function in
  // https://github.com/BrunoLevy/geogram/blob/main/src/lib/geogram/voronoi/convex_cell.cpp
  // Here, we need to find the intersection of the three planes using Cramer's
  // formula. The system of equations to be solved is:
  // pi.x * x + pi.y * y + pi.z * z = -pi.w
  // pj.x * x + pj.y * y + pj.z * z = -pj.w
  // pk.x * x + pk.y * y + pk.z * z = -pk.w
  //
  const vec4& pk = base;
  vec4 p;
  p.x = -det3x3(pi.w, pi.y, pi.z, pj.w, pj.y, pj.z, pk.w, pk.y, pk.z);
  p.y = -det3x3(pi.x, pi.w, pi.z, pj.x, pj.w, pj.z, pk.x, pk.w, pk.z);
  p.z = -det3x3(pi.x, pi.y, pi.w, pj.x, pj.y, pj.w, pk.x, pk.y, pk.w);
  p.w = +det3x3(pi.x, pi.y, pi.z, pj.x, pj.y, pj.z, pk.x, pk.y, pk.z);
  return p;
}

void PlanarVoronoiPolygon::get_properties(const device_vector<Vertex_t>& p,
                                          const device_vector<vec4>& planes,
                                          VoronoiCellProperties& props) const {
  if (p.size() < 3) return;
  // vec4 ah = compute(planes[p[0].bl], planes[p[0].br]);
  // vec4 bh = compute(planes[p[1].bl], planes[p[1].br]);
  vec4 ah = compute(planes[p[0]], planes[p[1]]);
  vec4 bh = compute(planes[p[1]], planes[p[2]]);
  vec3 a = (1.0 / ah.w) * ah.xyz();
  vec3 b = (1.0 / bh.w) * bh.xyz();

  for (size_t k = 2; k < p.size(); k++) {
    // vec4 ch = compute(planes[p[k].bl], planes[p[k].br]);
    vec4 ch = compute(planes[p[k]], planes[p[(k + 1) % p.size()]]);
    vec3 c = (1.0 / ch.w) * ch.xyz();
    coord_t ak = 0.5 * length(cross(b - a, c - a));
    vec3 ck = (1.0 / 3.0) * (a + b + c);

    props.moment = props.moment + ak * ck;
    props.volume += ak;
    b = c;
  }
}

void TriangulationDomain::initialize(
    int64_t elem, VoronoiPolygon<TriangulationDomain>& cell) const {
  vec3 vertices[3];
  for (int i = 0; i < 3; i++) {
    for (int d = 0; d < 3; d++)
      vertices[i][d] = points[3 * triangles[3 * elem + i] + d];
  }
  auto& planes = cell.planes();
  cell.cell().initialize(vertices, 3, cell.vertices(), planes);
}

namespace {

template <typename Domain_t>
class VoronoiThreadBlock : public VoronoiMesh {
  using Cell_t = VoronoiPolygon<Domain_t>;
  using Vertex_t = typename Cell_t::Vertex_t;
  using Element_t = typename Cell_t::Element_t;

 public:
  VoronoiThreadBlock(const Domain_t& domain, VoronoiCellMemoryPool& pool,
                     size_t tid)
      : VoronoiMesh(3),
        domain_(domain),  // copy the domain
        cell_(pool, tid),
        max_radius_(0) {
    allocate(20);
  }

  void append_to_mesh(VoronoiMesh& mesh) {
    // add vertices
    int n = mesh.vertices().n();
    mesh.vertices().reserve(mesh.vertices().n() + vertices_.n());
    mesh.triangles().reserve(mesh.triangles().n() + vertices_.n());
    for (size_t k = 0; k < vertices_.n(); k++) {
      mesh.vertices().add(vertices_[k]);
      bool on_boundary = false;
      for (int j = 0; j < 3; j++) {
        if (triangles_[k][j] == kMaxSite) on_boundary = true;
      }
      int group = 0;
      if (on_boundary) {
        for (int j = 0; j < 3; j++) triangles_[k][j] = 0;
        group = -1;
      }
      size_t id = mesh.triangles().n();
      mesh.triangles().add(triangles_[k]);
      mesh.triangles().set_group(id, group);
    }

    // add polygons
    mesh.polygons().reserve(mesh.polygons().n() + polygons_.n());
    std::vector<index_t> polygon(128);
    for (size_t k = 0; k < polygons_.n(); k++) {
      polygon.resize(polygons_.length(k));
      for (size_t j = 0; j < polygon.size(); j++)
        polygon[j] = polygons_[k][j] + n;
      mesh.polygons().add(polygon.data(), polygon.size());
      mesh.polygons().set_group(mesh.polygons().n() - 1, polygons_.group(k));
    }
  }

  void set_properties_ptr(VoronoiCellProperties* p) { properties_ = p; }
  void set_mesh_lock(std::mutex* lock) { append_mesh_lock_ = lock; }
  void set_status_ptr(VoronoiStatusCode* s) { status_ = s; }
  void set_weights_ptr(const coord_t* w) { weights_ = w; }
  Cell_t& cell() { return cell_; }
  double max_radius() const { return max_radius_; }

 protected:
  // only CPU
  std::mutex* append_mesh_lock_;

  // CPU or GPU
  Domain_t domain_;
  Cell_t cell_;
  VoronoiCellProperties* properties_{nullptr};
  VoronoiStatusCode* status_{nullptr};
  const coord_t* weights_;

  double max_radius_{0};
};

template <typename Domain_t>
class SiteThreadBlock : public VoronoiThreadBlock<Domain_t> {
  using Cell_t = VoronoiPolygon<Domain_t>;
  using Element_t = typename Cell_t::Element_t;
  using Base_t = VoronoiThreadBlock<Domain_t>;
  using Base_t::append_mesh_lock_;
  using Base_t::cell_;
  using Base_t::max_radius_;
  using Base_t::properties_;
  using Base_t::status_;
  using Base_t::vertices_;

 public:
  // only CPU
  SiteThreadBlock(const Domain_t& domain, VoronoiCellMemoryPool& pool,
                  size_t tid)
      : VoronoiThreadBlock<Domain_t>(domain, pool, tid) {}

  void compute(int dim, size_t m, size_t n, VoronoiMesh& mesh) {
    // compute all cells in range [m, n)
    if (mesh.save_mesh()) {
      vertices_.reserve((n - m) * 10);
      Base_t::template get<Element_t>().reserve(n - m);
    }
    for (size_t k = m; k < n; ++k) {
      // if (status_[k] == VoronoiStatusCode::kSuccess) continue;
      status_[k] = compute(dim, k, mesh);

      // save the maximum radius
      if (cell_.max_radius() > max_radius_) max_radius_ = cell_.max_radius();
    }
    if (mesh.save_mesh()) {
      append_mesh_lock_->lock();
      Base_t::append_to_mesh(mesh);
      append_mesh_lock_->unlock();
    }
  }

 private:
  VoronoiStatusCode compute(int dim, uint64_t site, VoronoiMesh& mesh) {
    auto status = Base_t::cell_.compute(Base_t::domain_, dim, site, *this);
    if (Base_t::properties_) {
      properties_[site].site = site;
      cell_.get_properties(properties_[site], true);
    }
    return status;
  }
};

template <typename Domain_t>
class ElementThreadBlock : public VoronoiThreadBlock<Domain_t> {
  using Cell_t = VoronoiPolygon<Domain_t>;
  using Element_t = typename Cell_t::Element_t;
  using Base_t = VoronoiThreadBlock<Domain_t>;
  using Base_t::append_mesh_lock_;
  using Base_t::cell_;
  using Base_t::domain_;
  using Base_t::max_radius_;
  using Base_t::properties_;
  using Base_t::status_;
  using Base_t::vertices_;

 public:
  ElementThreadBlock(const Domain_t& domain, VoronoiCellMemoryPool& pool,
                     size_t tid)
      : VoronoiThreadBlock<Domain_t>(domain, pool, tid) {}

  void compute(int dim, size_t m, size_t n, VoronoiMesh& mesh) {
    // compute all cells in range [m, n)
    if (mesh.save_mesh()) {
      vertices_.reserve((n - m) * 10);
      Base_t::template get<Element_t>().reserve(n - m);
    }
    for (size_t k = m; k < n; ++k) {
      // if (status_[k] == VoronoiStatusCode::kSuccess) continue;
      status_[k] = compute(dim, k, mesh);
    }
    if (mesh.save_mesh()) {
      append_mesh_lock_->lock();
      Base_t::append_to_mesh(mesh);
      append_mesh_lock_->unlock();
    }
  }

  void set_elem2site_ptr(const index_t* elem2site) { elem2site_ = elem2site; }

 private:
  VoronoiStatusCode compute(int dim, uint64_t elem, VoronoiMesh& mesh) {
    auto status =
        cell_.compute(domain_, dim, elem, elem2site_[elem], properties_, *this);
    return status;
  }

  const index_t* elem2site_{nullptr};
};

template <typename ThreadBlock_t>
void clip(ThreadBlock_t* block, int dim, size_t m, size_t n,
          VoronoiMesh* mesh) {
  block->compute(dim, m, n, *mesh);
}

template <int dim>
std::shared_ptr<trees::KdTreeNd<coord_t, index_t>> get_nearest_neighbor(
    const coord_t* p, uint64_t np, const coord_t* q, uint64_t nq,
    std::vector<index_t>& nn, const VoronoiDiagramOptions& options,
    std::shared_ptr<trees::KdTreeNd<coord_t, index_t>> ptree = nullptr) {
  Timer timer;
  trees::KdTreeOptions kdtree_opts;
  kdtree_opts.max_dim = options.max_kdtree_axis_dim;
  if (kdtree_opts.max_dim < 0) kdtree_opts.max_dim = dim;
  using kdtree_t = trees::KdTree<dim, coord_t, index_t>;
  if (!ptree) {
    timer.start();
    ptree = std::make_shared<kdtree_t>(p, np, kdtree_opts);
    timer.stop();
    if (options.verbose)
      LOG << "kdtree created in " << timer.seconds() << " s.";
  }

  auto* tree = static_cast<kdtree_t*>(ptree.get());
  timer.start();
  std::parafor_i(0, nq,
                 [&](int tid, int k) { nn[k] = tree->nearest(&q[k * dim]); });
  timer.stop();
  if (options.verbose)
    LOG << "nearest neighbors computed in " << timer.seconds() << " s.";
  return ptree;
}

std::shared_ptr<trees::KdTreeNd<coord_t, index_t>> get_kdtree(
    const coord_t* p, uint64_t np, int dim,
    const VoronoiDiagramOptions& options) {
  trees::KdTreeOptions kdtree_opts;
  kdtree_opts.max_dim = options.max_kdtree_axis_dim;
  if (kdtree_opts.max_dim < 0) kdtree_opts.max_dim = dim;
  if (dim == 2) {
    return std::make_shared<trees::KdTree<2, coord_t, index_t>>(p, np,
                                                                kdtree_opts);
  } else if (dim == 3) {
    return std::make_shared<trees::KdTree<3, coord_t, index_t>>(p, np,
                                                                kdtree_opts);
  } else if (dim == 4) {
    return std::make_shared<trees::KdTree<4, coord_t, index_t>>(p, np,
                                                                kdtree_opts);
  } else
    NOT_IMPLEMENTED;
  return nullptr;
}

}  // namespace

template <int dim>
std::shared_ptr<trees::KdTreeNd<coord_t, index_t>> get_nearest_neighbors(
    const coord_t* p, uint64_t np, const coord_t* q, uint64_t nq,
    std::vector<index_t>& knn, size_t n_neighbors,
    const VoronoiDiagramOptions& options, VoronoiStatistics& stats,
    std::shared_ptr<trees::KdTreeNd<coord_t, index_t>> ptree) {
  Timer timer;
  trees::KdTreeOptions kdtree_opts;
  kdtree_opts.max_dim = options.max_kdtree_axis_dim;
  if (kdtree_opts.max_dim < 0) kdtree_opts.max_dim = dim;
  using kdtree_t = trees::KdTree<dim, coord_t, index_t>;
  if (!ptree) {
    timer.start();
    ptree = std::make_shared<kdtree_t>(p, np, kdtree_opts);
    timer.stop();
    if (options.verbose)
      LOG << "kdtree created in " << timer.seconds() << " s.";
    stats.t_kdtree_build = timer.seconds();
  }

  auto* tree = static_cast<kdtree_t*>(ptree.get());
  timer.start();
  std::parafor_i(0, nq, [&](int tid, int k) {
    index_t* neighbors = (index_t*)alloca(n_neighbors * sizeof(index_t));
    coord_t* distances = (coord_t*)alloca(n_neighbors * sizeof(coord_t));
    trees::NearestNeighborSearch<index_t, coord_t> search(n_neighbors,
                                                          neighbors, distances);
    tree->knearest(&q[k * dim], search);
    for (size_t j = 0; j < n_neighbors; ++j)
      knn[k * n_neighbors + j] = neighbors[j];
  });
  timer.stop();
  if (options.verbose)
    LOG << "nearest neighbors computed in " << timer.seconds() << " s.";
  stats.t_kdtree_query = timer.seconds();
  return ptree;
}

VoronoiDiagram::VoronoiDiagram(int dim, const coord_t* sites, uint64_t n_sites)
    : VoronoiMesh(3),
      dim_(dim),
      sites_(sites),
      n_sites_(n_sites),
      neighbors_(*this, sites, dim) {}

void VoronoiDiagram::create_sqtree(int n_subdiv) {
  sqtree_ = std::make_unique<SphereQuadtree>(sites_, n_sites_, dim_, n_subdiv);
}

int check_closed_delaunay(
    const absl::flat_hash_set<std::array<uint32_t, 3>>& triangles) {
  absl::flat_hash_set<std::pair<uint32_t, uint32_t>> edges;
  for (const auto& t : triangles) {
    for (int j = 0; j < 3; j++) {
      uint32_t p = t[j];
      uint32_t q = t[j == 2 ? 0 : j + 1];
      if (p > q) std::swap(p, q);
      auto it = edges.find({p, q});
      if (it == edges.end())
        edges.insert({p, q});
      else
        edges.erase(it);
    }
  }
  return edges.size();
}

template <typename Domain_t>
void VoronoiDiagram::compute(const Domain_t& domain,
                             VoronoiDiagramOptions options) {
  ASSERT(sites_);
  using ThreadBlock_t = SiteThreadBlock<Domain_t>;
  Timer timer;
  size_t n_threads = std::thread::hardware_concurrency();

  // reset all the statistics
  statistics_.reset();
  statistics_.energy = -1;  // not computed yet
  statistics_.n_neighbors = options.n_neighbors;
  average_statistics_.n_sites = n_sites_;
  statistics_.n_sites = n_sites_;

  Timer global_timer;
  global_timer.start();

  auto alg = options.neighbor_algorithm;
  bool is_sphere = std::is_same<Domain_t, SphereDomain>::value;
  bool use_kdtree = alg == NearestNeighborAlgorithm::kKdtree;
  bool use_bfs = alg == NearestNeighborAlgorithm::kVoronoiBFS;
  bool use_sqtree =
      alg == NearestNeighborAlgorithm::kSphereQuadtree && is_sphere;
  if (use_bfs && facets_.empty()) use_kdtree = true;

  // calculate nearest neighbors
  size_t n_neighbors = options.n_neighbors;
  if (n_sites_ < n_neighbors) n_neighbors = n_sites_;
  std::vector<index_t> knn(n_sites_ * n_neighbors);
  std::vector<uint16_t> max_neighbors(n_sites_, options.n_neighbors);
  std::shared_ptr<trees::KdTreeNd<coord_t, index_t>> tree{nullptr};
  if (use_kdtree) {
    if (dim_ == 2)
      tree = get_nearest_neighbors<2>(sites_, n_sites_, sites_, n_sites_, knn,
                                      n_neighbors, options, statistics_);
    else if (dim_ == 3)
      tree = get_nearest_neighbors<3>(sites_, n_sites_, sites_, n_sites_, knn,
                                      n_neighbors, options, statistics_);
    else if (dim_ == 4)
      tree = get_nearest_neighbors<4>(sites_, n_sites_, sites_, n_sites_, knn,
                                      n_neighbors, options, statistics_);
    else
      NOT_IMPLEMENTED;
  } else if (use_sqtree) {
    if (!sqtree_) create_sqtree(-1);

    timer.start();
    sqtree_->build();
    timer.stop();
    double t_build = timer.seconds();
    statistics_.n_sqtree_maxleaf = sqtree_->min_leaf_size();
    statistics_.n_sqtree_maxleaf = sqtree_->max_leaf_size();
    statistics_.t_sqtree_build = t_build;

    std::vector<SphereQuadtreeWorkspace> searches(n_threads,
                                                  options.n_neighbors);
    timer.start();
    std::parafor_i(0, n_sites_,
                   [this, &searches, &knn, &options, &max_neighbors](size_t tid,
                                                                     size_t k) {
                     sqtree_->knearest(k, searches[tid]);
                     const auto& result = searches[tid].neighbors;
                     size_t m = searches[tid].size();
                     if (searches[tid].size() > options.n_neighbors)
                       m = options.n_neighbors;
                     ASSERT(result[0].first == k);
                     for (size_t j = 0; j < m; j++) {
                       knn[options.n_neighbors * k + j] = result[j].first;
                     }
                     max_neighbors[k] = m;
                   });
    timer.stop();
    statistics_.t_sqtree_query = timer.seconds();

    if (options.verbose)
      LOG << fmt::format(
          "neighbors computed in {} s. (build = {} s.), total = {} s.",
          timer.seconds(), t_build, timer.seconds() + t_build);

    // still build the kdtree to retrieve neighbors if the cell is incomplete
    timer.start();
    tree = get_kdtree(sites_, n_sites_, dim_, options);
    timer.stop();
    statistics_.t_kdtree_build = timer.seconds();
  } else if (use_bfs) {
    timer.start();
    neighbors_.build();
    timer.stop();
    double t_build = timer.seconds();
    statistics_.t_bfs_build = t_build;
    statistics_.n_bfs_level = options.bfs_max_level;

    std::vector<NearestNeighborsWorkspace> searches(n_threads,
                                                    options.n_neighbors);
    timer.start();
    std::parafor_i(0, n_sites_,
                   [this, &searches, &knn, &options, &max_neighbors](size_t tid,
                                                                     size_t k) {
                     searches[tid].max_level = options.bfs_max_level;
                     neighbors_.knearest(k, searches[tid]);
                     const auto& result = searches[tid].sites;
                     size_t m = result.size();
                     if (result.size() > options.n_neighbors)
                       m = options.n_neighbors;
                     ASSERT(result.data()[0].first == k);
                     for (size_t j = 0; j < m; j++) {
                       knn[options.n_neighbors * k + j] =
                           result.data()[j].first;
                     }
                     max_neighbors[k] = result.size();
                   });

    timer.stop();
    statistics_.t_bfs_query = timer.seconds();

    if (options.verbose)
      LOG << fmt::format(
          "neighbors computed in {} s. (build = {} s.), total = {} s.",
          timer.seconds(), t_build, timer.seconds() + t_build);

    // still build the kdtree to retrieve neighbors if the cell is incomplete
    timer.start();
    tree = get_kdtree(sites_, n_sites_, dim_, options);
    timer.stop();
    statistics_.t_kdtree_build = timer.seconds();
  } else {
    LOG << "unknown nearest neighbor algorithm";
    NOT_POSSIBLE;
  }

  // always add sites first for the Delaunay triangulation
  if (options.store_mesh) {
    vertices_.clear();
    for (size_t k = 0; k < n_sites_; k++) vertices_.add(sites_ + k * dim_);
  }

  // compute voronoi diagram
  timer.start();
  if (!options.parallel) n_threads = 1;
  // n_threads = 1;
  if (n_threads > VORTEX_NUM_CORES) n_threads = VORTEX_NUM_CORES;
  std::mutex append_mesh_lock;
  set_save_mesh(options.store_mesh);
  set_save_facets(options.store_facet_data);
  set_save_delaunay(options.store_delaunay_triangles);
  facets_.clear();
  delaunay_.clear();
  polygons_.clear();
  triangles_.clear();
  lines_.clear();

  std::vector<std::thread> threads;
  std::vector<std::shared_ptr<ThreadBlock_t>> blocks;
  properties_.resize(n_sites_);

  VoronoiCellMemoryPool pool;
  status_.resize(n_sites_, VoronoiStatusCode::kIncomplete);
  size_t n_sites_per_block = n_sites_ / n_threads;
  for (size_t k = 0; k < n_threads; k++) {
    blocks.push_back(std::make_shared<ThreadBlock_t>(domain, pool, k));
    blocks[k]->cell().set_sites(sites_, n_sites_);
    blocks[k]->cell().set_neighbors(knn.data(), n_neighbors);
    blocks[k]->cell().set_max_neighbors(max_neighbors.data());
    blocks[k]->cell().set_kdtree(tree.get());
    blocks[k]->set_mesh_lock(&append_mesh_lock);
    blocks[k]->set_properties_ptr(properties_.data());
    blocks[k]->set_status_ptr(status_.data());
    if (weights_.size() == n_sites_) {
      blocks[k]->cell().set_weights(weights_.data());
    }
    blocks[k]->set_save_mesh(options.store_mesh);
    blocks[k]->set_save_facets(options.store_facet_data);
    blocks[k]->set_save_delaunay(options.store_delaunay_triangles);
    size_t m = k * n_sites_per_block;
    size_t n = m + n_sites_per_block;
    if (k + 1 == n_threads) n = n_sites_;
    std::thread t(clip<ThreadBlock_t>, blocks[k].get(), dim_, m, n, this);
    threads.push_back(std::move(t));
  }

  for (auto& t : threads) t.join();
  timer.stop();
  if (options.verbose) {
    LOG << "voronoi computed in " << timer.seconds() << " s.";
  }
  statistics_.t_voronoi = timer.seconds();

  // save the max radius
  max_radius_ = 0;
  for (const auto& block : blocks) {
    if (block->max_radius() > max_radius_) max_radius_ = block->max_radius();
  }

  // merge the facets
  if (options.store_facet_data) {
    timer.start();
    size_t n_facets = 0;
    for (const auto& block : blocks) n_facets += block->facets().size();
    facets_.reserve(n_facets);
    for (auto& block : blocks) append(*block);
    timer.stop();
    if (options.verbose)
      LOG << fmt::format("# facets = {}, (done in {} sec.)", facets_.size(),
                         timer.seconds());
    statistics_.t_facets = timer.seconds();
  }

  // store the delaunay triangles
  if (options.store_delaunay_triangles) {
    timer.start();
    size_t n_triangles = 0;
    for (const auto& block : blocks) n_triangles += block->delaunay().size();
    delaunay_.reserve(n_triangles);
    for (auto& block : blocks) append_triangles(*block);
    timer.stop();
    if (options.verbose)
      LOG << fmt::format("# triangles = {}, (done in {} sec.)",
                         delaunay_.size(), timer.seconds());
    statistics_.t_delaunay = timer.seconds();

    if (options.check_closed) {
      timer.start();
      statistics_.n_bnd_delaunay_edges = check_closed_delaunay(delaunay_);
      timer.stop();
      if (options.verbose)
        LOG << fmt::format("delaunay check done in {} sec.", timer.seconds());
    }
  }

  uint64_t n_incomplete = 0;
  for (auto s : status_) {
    if (s == VoronoiStatusCode::kRadiusNotReached) n_incomplete++;
  }
  if (options.verbose && n_incomplete > 0)
    LOG << fmt::format("n_incomplete = {}", n_incomplete);
  statistics_.n_incomplete = n_incomplete;

  statistics_.area = analyze().area;
  statistics_.area_error = statistics_.area - domain.area();

  global_timer.stop();
  statistics_.t_total = global_timer.seconds();
  average_statistics_.append_to_average(statistics_);
  if (track_statistics_history_) statistics_history_.push_back(statistics_);
}

template <>
void VoronoiDiagram::compute(const TriangulationDomain& domain,
                             VoronoiDiagramOptions options) {
  ASSERT(sites_);
  using ThreadBlock_t = ElementThreadBlock<TriangulationDomain>;

  // calculate nearest neighbors
  size_t n_neighbors = options.n_neighbors;
  if (n_sites_ < n_neighbors) n_neighbors = n_sites_;
  std::vector<index_t> knn(n_sites_ * n_neighbors);
  std::shared_ptr<trees::KdTreeNd<coord_t, index_t>> tree{nullptr};
  if (dim_ == 2)
    tree = get_nearest_neighbors<2>(sites_, n_sites_, sites_, n_sites_, knn,
                                    n_neighbors, options, statistics_);
  else if (dim_ == 3)
    tree = get_nearest_neighbors<3>(sites_, n_sites_, sites_, n_sites_, knn,
                                    n_neighbors, options, statistics_);
  else if (dim_ == 4)
    tree = get_nearest_neighbors<4>(sites_, n_sites_, sites_, n_sites_, knn,
                                    n_neighbors, options, statistics_);
  else
    NOT_IMPLEMENTED;

  // compute voronoi diagram
  Timer timer;
  timer.start();
  size_t n_threads = std::thread::hardware_concurrency();
  if (!options.parallel) n_threads = 1;
  std::mutex append_mesh_lock;
  allocate(n_sites_);
  set_save_mesh(options.store_mesh);
  set_save_facets(options.store_facet_data);
  std::vector<std::thread> threads;
  std::vector<std::shared_ptr<ThreadBlock_t>> blocks;
  properties_.resize(n_sites_);

  // nearest neighbor to vertices of the mesh (recycle existing tree)
  std::vector<index_t> vnn(domain.n_points);
  if (dim_ == 2)
    get_nearest_neighbor<2>(sites_, n_sites_, domain.points, domain.n_points,
                            vnn, options, tree);
  else if (dim_ == 3)
    get_nearest_neighbor<3>(sites_, n_sites_, domain.points, domain.n_points,
                            vnn, options, tree);
  else if (dim_ == 4)
    get_nearest_neighbor<4>(sites_, n_sites_, domain.points, domain.n_points,
                            vnn, options, tree);

  // use triangle nearest neighbor as nearest neighbor of first vertex
  std::vector<index_t> tnn(domain.n_triangles);
  for (size_t k = 0; k < domain.n_triangles; k++)
    tnn[k] = vnn[domain.triangles[3 * k]];

  // reset properties
  for (auto& props : properties_) props.reset();

  // always add sites first for the Delaunay triangulation
  if (options.store_mesh) {
    for (size_t k = 0; k < n_sites_; k++) vertices_.add(sites_ + k * dim_);
  }

  // set up the thread blocks
  size_t n_elems = domain.n_elems();
  status_.resize(n_elems, VoronoiStatusCode::kIncomplete);
  if (n_threads >= n_elems) n_threads = 1;
  if (n_threads > VORTEX_NUM_CORES) n_threads = VORTEX_NUM_CORES;
  // n_threads = 1;
  VoronoiCellMemoryPool pool;
  size_t n_elems_per_block = n_elems / n_threads;
  for (size_t k = 0; k < n_threads; k++) {
    blocks.push_back(std::make_shared<ThreadBlock_t>(domain, pool, k));
    blocks[k]->cell().set_sites(sites_, n_sites_);
    blocks[k]->cell().set_neighbors(knn.data(), n_neighbors);
    blocks[k]->set_mesh_lock(&append_mesh_lock);
    blocks[k]->set_elem2site_ptr(tnn.data());
    blocks[k]->set_properties_ptr(properties_.data());
    blocks[k]->set_status_ptr(status_.data());
    blocks[k]->cell().set_kdtree(tree.get());
    blocks[k]->set_save_mesh(options.store_mesh);
    blocks[k]->set_save_facets(options.store_facet_data);
    size_t m = k * n_elems_per_block;
    size_t n = m + n_elems_per_block;
    if (k + 1 == n_threads) n = n_elems;
    std::thread t(clip<ThreadBlock_t>, blocks[k].get(), dim_, m, n, this);
    threads.push_back(std::move(t));
  }

  for (auto& t : threads) t.join();
  timer.stop();
  if (options.verbose)
    LOG << "voronoi computed in " << timer.seconds() << " s.";

  // save the max radius
  max_radius_ = 0;
  for (const auto& block : blocks) {
    if (block->max_radius() > max_radius_) max_radius_ = block->max_radius();
  }

  uint64_t n_incomplete = 0;
  for (auto s : status_) {
    if (s == VoronoiStatusCode::kRadiusNotReached) n_incomplete++;
  }
  if (options.verbose && n_incomplete > 0)
    LOG << fmt::format("n_incomplete = {}", n_incomplete);
}

void lift_sites(Vertices& sites, const std::vector<double>& weights) {
  ASSERT(sites.n() == weights.size() && !weights.empty());
  ASSERT(sites.dim() == 4);
  double wmax = *std::max_element(weights.begin(), weights.end());
  for (size_t k = 0; k < sites.n(); k++) {
    sites[k][3] = std::sqrt(std::max(0.0, wmax - weights[k]));
  }
}

void VoronoiDiagram::smooth(Vertices& sites, bool on_sphere) const {
  vec3 x;
  for (size_t k = 0; k < n_sites_; k++) {
    if (properties_[k].volume == 0) continue;
    ASSERT(properties_[k].volume > 0);
    x = static_cast<float>(1.0 / properties_[k].volume) * properties_[k].moment;
    if (on_sphere) x = unit_vector(x);
    for (int d = 0; d < 3; d++) sites[k][d] = x[d];
  }
}

VoronoiDiagramProperties VoronoiDiagram::analyze() const {
  // TODO(philip) calculate energy and gradient norms
  VoronoiDiagramProperties props;
  for (size_t k = 0; k < n_sites_; k++) {
    props.area += properties_[k].volume;
  }
  return props;
}

void VoronoiDiagram::merge() {
  Vertices vertices(vertices_.dim());
  vertices.reserve(vertices_.n());

  std::unordered_map<std::array<index_t, 3>, index_t> tmap;
  std::vector<index_t> vmap(vertices_.n());
  tmap.reserve(vertices_.n());

  // copy the sites
  for (size_t k = 0; k < n_sites_; k++) vertices.add(vertices_[k]);

  Topology<Triangle> triangles;
  ASSERT(vertices_.n() == triangles_.n() + n_sites_);
  for (size_t k = 0; k < triangles_.n(); k++) {
    if (triangles_.group(k) < 0) {
      // boundary vertex
      vmap[k + n_sites_] = vertices.n();
      vertices.add(vertices_[k + n_sites_]);
      continue;
    }
    const auto* t = triangles_[k];
    std::array<index_t, 3> triangle = {t[0], t[1], t[2]};
    std::sort(triangle.begin(), triangle.end());
    auto it = tmap.find(triangle);
    if (it == tmap.end()) {
      tmap.insert({triangle, vertices.n()});
      vertices.add(vertices_[k + n_sites_]);
      triangles.add(t);
    }
    vmap[k + n_sites_] = tmap.at(triangle);
  }
  LOG << fmt::format("found {} unique vertices", tmap.size());

  for (size_t k = 0; k < polygons_.n(); k++) {
    for (int j = 0; j < polygons_.length(k); j++) {
      polygons_[k][j] = vmap.at(polygons_[k][j]);
    }
  }
  vertices.copy(vertices_);
  triangles.copy(triangles_);
}

nlohmann::json VoronoiStatistics::to_json() const {
  nlohmann::json data;

  data["n_neighbors"] = n_neighbors;
  data["t_kdtree_build"] = t_kdtree_build;
  data["t_kdtree_query"] = t_kdtree_query;
  data["t_bfs_build"] = t_bfs_build;
  data["t_bfs_query"] = t_bfs_query;
  data["t_sqtree_build"] = t_sqtree_build;
  data["t_sqtree_query"] = t_sqtree_query;
  data["n_sqtree_minleaf"] = n_sqtree_minleaf;
  data["n_sqtree_maxleaf"] = n_sqtree_maxleaf;
  data["t_voronoi"] = t_voronoi;
  data["n_incomplete"] = n_incomplete;
  data["t_facets"] = t_facets;
  data["t_delaunay"] = t_delaunay;
  data["energy"] = energy;
  data["t_total"] = t_total;
  data["count"] = count;
  data["area_error"] = area_error;
  data["area"] = area;
  data["n_bfs_level"] = n_bfs_level;
  data["n_sites"] = n_sites;
  data["n_triangles"] = n_triangles;
  data["n_bnd_delaunay_edges"] = n_bnd_delaunay_edges;

  return data;
}

template void VoronoiDiagram::compute(const SphereDomain&,
                                      VoronoiDiagramOptions);
template void VoronoiDiagram::compute(const SquareDomain&,
                                      VoronoiDiagramOptions);

#define INSTANTIATE_NEAREST_NEIGHBORS(DIM)                            \
  template std::shared_ptr<trees::KdTreeNd<coord_t, index_t>>         \
  get_nearest_neighbors<DIM>(                                         \
      const coord_t* p, uint64_t np, const coord_t* q, uint64_t nq,   \
      std::vector<index_t>& knn, size_t n_neighbors,                  \
      const VoronoiDiagramOptions& options, VoronoiStatistics& stats, \
      std::shared_ptr<trees::KdTreeNd<coord_t, index_t>> ptree);

INSTANTIATE_NEAREST_NEIGHBORS(2)
INSTANTIATE_NEAREST_NEIGHBORS(3)
INSTANTIATE_NEAREST_NEIGHBORS(4)

#undef INSTANTIATE_NEAREST_NEIGHBORS

}  // namespace vortex