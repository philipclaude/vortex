#include "voronoi.h"

#include <Predicates_psm.h>

#include <memory>
#include <set>
#include <thread>
#include <unordered_set>

#include "elements.h"
// #define HAVE_NANOFLANN 1
#include "kdtree.h"
#include "mesh.h"
#include "stlext.h"
#include "voronoi_polygon.hpp"

namespace vortex {

void SphereDomain::initialize(vec4 site,
                              VoronoiPolygon<SphereDomain>& polygon) const {
  auto& vertices = polygon.vertices();
  auto& planes = polygon.planes();
  vec3 center = unit_vector(site.xyz());
  const double r = 0.7 * radius;  // sqrt(2)/2 would graze the sphere
  polygon.cell().center = center;

  // compute the normal and tangent to the sphere
  vec3 n = unit_vector(center);
  vec3 t = unit_vector(vec3{-n.y, n.x, 0.0});
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
  vertices[0].bl = 0;
  vertices[0].br = 1;
  vertices[1].bl = 1;
  vertices[1].br = 2;
  vertices[2].bl = 2;
  vertices[2].br = 3;
  vertices[3].bl = 3;
  vertices[3].br = 0;
}

vec4 SphericalVoronoiPolygon::plane_equation(const vec4& ui, const vec4& uj) {
  coord_t wi = ui.w;
  coord_t wj = uj.w;

  // https://en.wikipedia.org/wiki/Radical_axis
  coord_t d = length((uj - ui).xyz());
  coord_t d1 = (d * d + wi - wj) / (2.0 * d);
  vec4 m = ui + d1 * (uj - ui) / d;
  vec3 n = unit_vector((ui - uj).xyz());
  return {n.x, n.y, n.z, -dot(n, m.xyz())};
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
    const pool<Vertex_t>& p, const pool<vec4>& planes,
    VoronoiCellProperties& props) const {
  props.mass = 0.0;
  props.moment = {0, 0, 0};
  vec3 a = compute(planes[p[0].bl], planes[p[0].br]).xyz();
  for (size_t i = 0; i < p.size(); i++) {
    size_t j = (i + 1) == p.size() ? 0 : i + 1;
    vec3 b = compute(planes[p[j].bl], planes[p[j].br]).xyz();
    const vec3& c = center;

    // https://www.johndcook.com/blog/2021/11/29/area-of-spherical-triangle/
    coord_t num = std::fabs(dot(a, cross(b, c)));
    coord_t den = 1.0 + dot(a, b) + dot(b, c) + dot(a, c);
    coord_t ak = 2.0 * std::atan2(num, den);
    vec3 ck = unit_vector((1.0 / 3.0) * (a + b + c));

    props.moment = props.moment + ak * ck;
    props.mass += ak;
    a = b;
  }
}

namespace {

template <typename Domain_t>
class VoronoiThreadBlock : public Mesh {
  using Cell_t = VoronoiPolygon<Domain_t>;
  using Vertex_t = typename Cell_t::Vertex_t;
  using Element_t = typename Cell_t::Element_t;

 public:
  // only CPU
  VoronoiThreadBlock(const Domain_t& domain, size_t nv = 512, size_t np = 512)
      : Mesh(3),
        vertex_pool_(nv),
        plane_pool_(np),
        domain_(domain),  // copy the domain
        cell_(&vertex_pool_[0], nv, &plane_pool_[0], np) {}

  VoronoiStatusCode compute(int dim, uint64_t k, Mesh* mesh) {
    auto status = cell_.compute(domain_, dim, k);
    if (properties_) {
      cell_.get_properties(properties_[k]);
      properties_[k].site = k;
    }
    if (mesh) cell_.append_to_mesh(*this);
    return status;
  }

  void compute(int dim, size_t m, size_t n, Mesh* mesh) {
    // compute all cells in range [m, n)
    if (mesh) {
      vertices_.reserve((n - m) * 10);
      get<Element_t>().reserve(n - m);
    }
    for (size_t k = m; k < n; ++k) {
      // if (status_[k] == VoronoiStatusCode::kSuccess) continue;
      status_[k] = compute(dim, k, mesh);
    }
    if (mesh) {
      append_mesh_lock_->lock();
      append_to_mesh(*mesh);
      append_mesh_lock_->unlock();
    }
  }

  void append_to_mesh(Mesh& mesh) {
    // add vertices
    int n = mesh.vertices().n();
    mesh.vertices().reserve(mesh.vertices().n() + vertices_.n());
    for (int k = 0; k < vertices_.n(); k++) {
      mesh.vertices().add(vertices_[k]);
    }

    // add polygons
    mesh.polygons().reserve(mesh.polygons().n() + polygons_.n());
    std::vector<index_t> polygon(128);
    for (int k = 0; k < polygons_.n(); k++) {
      polygon.resize(polygons_.length(k));
      for (size_t j = 0; j < polygon.size(); j++)
        polygon[j] = polygons_[k][j] + n;
      mesh.polygons().add(polygon.data(), polygon.size());
    }
  }

  void set_properties_ptr(VoronoiCellProperties* p) { properties_ = p; }
  void set_mesh_lock(std::mutex* lock) { append_mesh_lock_ = lock; }
  void set_status_ptr(VoronoiStatusCode* s) { status_ = s; }
  Cell_t& cell() { return cell_; }

 private:
  // only CPU
  std::vector<Vertex_t> vertex_pool_;
  std::vector<vec4> plane_pool_;
  std::mutex* append_mesh_lock_;

  // CPU or GPU
  Domain_t domain_;
  Cell_t cell_;
  VoronoiCellProperties* properties_{nullptr};
  VoronoiStatusCode* status_{nullptr};
};

template <typename ThreadBlock_t>
void clip(ThreadBlock_t* block, int dim, size_t m, size_t n, Mesh* mesh) {
  block->compute(dim, m, n, mesh);
}

template <int dim>
std::shared_ptr<maple::KdTreeNd<coord_t, index_t>> get_nearest_neighbors(
    const coord_t* p, uint64_t np, const coord_t* q, uint64_t nq,
    std::vector<index_t>& knn, size_t n_neighbors,
    const VoronoiDiagramOptions& options,
    std::shared_ptr<maple::KdTreeNd<coord_t, index_t>> ptree = nullptr) {
  size_t n_threads = std::thread::hardware_concurrency();
  Timer timer;
  maple::KdTreeOptions kdtree_opts;
  kdtree_opts.max_dim = options.max_kdtree_axis_dim;
  if (kdtree_opts.max_dim < 0) kdtree_opts.max_dim = dim;
  using kdtree_t = maple::KdTree<dim, coord_t, index_t>;
  //  using kdtree_t = KdTree_nanoflann<dim, coord_t, index_t>;
  if (!ptree) {
    timer.start();
    ptree = std::make_shared<kdtree_t>(p, np, kdtree_opts);
    timer.stop();
    if (options.verbose)
      LOG << "kdtree created in " << timer.seconds() << " s.";
  }
  if (options.interleave_neighbors) return ptree;

  auto* tree = static_cast<kdtree_t*>(ptree.get());
  timer.start();
  std::parafor_i(0, nq, [&](int tid, int k) {
    index_t* neighbors = (index_t*)alloca(n_neighbors * sizeof(index_t));
    coord_t* distances = (coord_t*)alloca(n_neighbors * sizeof(coord_t));
    maple::NearestNeighborSearch<index_t, coord_t> search(n_neighbors,
                                                          neighbors, distances);
    tree->knearest(&q[k * dim], search);
    for (int j = 0; j < n_neighbors; ++j)
      knn[k * n_neighbors + j] = neighbors[j];
  });
  timer.stop();
  if (options.verbose)
    LOG << "nearest neighbors computed in " << timer.seconds() << " s.";
  return ptree;
}

}  // namespace

template <typename Domain_t>
void VoronoiDiagram::compute(const Domain_t& domain,
                             VoronoiDiagramOptions options) {
  ASSERT(sites_);
  using ThreadBlock_t = VoronoiThreadBlock<Domain_t>;

  // calculate nearest neighbors
  size_t n_neighbors = options.n_neighbors;
  if (n_sites_ < n_neighbors) n_neighbors = n_sites_;
  std::vector<index_t> knn(n_sites_ * n_neighbors);
  std::shared_ptr<maple::KdTreeNd<coord_t, index_t>> tree{nullptr};
  if (dim_ == 2)
    tree = get_nearest_neighbors<2>(sites_, n_sites_, sites_, n_sites_, knn,
                                    n_neighbors, options);
  else if (dim_ == 3)
    tree = get_nearest_neighbors<3>(sites_, n_sites_, sites_, n_sites_, knn,
                                    n_neighbors, options);
  else if (dim_ == 4)
    tree = get_nearest_neighbors<4>(sites_, n_sites_, sites_, n_sites_, knn,
                                    n_neighbors, options);
  else
    NOT_IMPLEMENTED;

  // compute voronoi diagram
  bool set_kdtree = options.interleave_neighbors;
  bool recomputing = false;
compute_voronoi:
  Timer timer;
  timer.start();
  size_t n_threads = std::thread::hardware_concurrency();
  if (!options.parallel) n_threads = 1;
  std::mutex append_mesh_lock;
  Mesh* mesh = options.store_mesh ? this : nullptr;
  // (options.store_mesh) ? ((options.mesh) ? options.mesh : this) : nullptr;
  std::vector<std::thread> threads;
  std::vector<std::shared_ptr<ThreadBlock_t>> blocks;
  properties_.resize(n_sites_);

  if (!recomputing) status_.resize(n_sites_, VoronoiStatusCode::kIncomplete);
  size_t n_sites_per_block = n_sites_ / n_threads;
  for (size_t k = 0; k < n_threads; k++) {
    blocks.push_back(std::make_shared<ThreadBlock_t>(domain));
    blocks[k]->cell().set_sites(sites_);
    blocks[k]->cell().set_neighbors(knn.data(), n_neighbors);
    blocks[k]->set_mesh_lock(&append_mesh_lock);
    blocks[k]->set_properties_ptr(properties_.data());
    blocks[k]->set_status_ptr(status_.data());
    if (set_kdtree) blocks[k]->cell().set_kdtree(tree.get());
    size_t m = k * n_sites_per_block;
    size_t n = m + n_sites_per_block;
    if (k + 1 == n_threads) n = n_sites_;
    std::thread t(clip<ThreadBlock_t>, blocks[k].get(), dim_, m, n, mesh);
    threads.push_back(std::move(t));
  }

  for (auto& t : threads) t.join();
  timer.stop();
  if (options.verbose)
    LOG << "voronoi computed in " << timer.seconds() << " s.";

  uint64_t n_incomplete = 0;
  for (auto s : status_) {
    if (s == VoronoiStatusCode::kRadiusNotReached) n_incomplete++;
  }
  if (options.verbose) LOG << fmt::format("n_incomplete = {}", n_incomplete);
  if (n_incomplete > 0 && options.allow_reattempt) {
    set_kdtree = true;
    recomputing = true;
    goto compute_voronoi;
  }
}

template void VoronoiDiagram::compute(const SphereDomain&,
                                      VoronoiDiagramOptions);

}  // namespace vortex