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
#pragma once

#include <absl/container/flat_hash_map.h>
#include <absl/container/flat_hash_set.h>
#include <stlext.h>

#include <cassert>
#include <cmath>
#include <cstdint>
#include <nlohmann/json_fwd.hpp>

#include "defs.h"
#include "device_util.h"
#include "log.h"
#include "mesh.h"
#include "neighbors.h"
#include "util.h"

#ifndef VORTEX_NUM_CORES
#define VORTEX_NUM_CORES 1
#endif

#define NEW_FACETS 1

namespace trees {
template <typename coord_t, typename index_t>
class KdTreeNd;
}  // namespace trees
namespace vortex {

struct VoronoiStatistics {
  VoronoiStatistics() { reset(); }
  int n_neighbors;
  int n_sites;
  int n_triangles;
  double t_kdtree_build;
  double t_kdtree_query;
  double t_bfs_build;
  double t_bfs_query;
  double t_sqtree_build;
  double t_sqtree_query;
  int n_sqtree_minleaf;
  int n_sqtree_maxleaf;
  int n_bfs_level;
  double t_voronoi;
  double n_incomplete;
  double t_facets;
  double t_delaunay;
  double energy;
  double t_total;
  int count = 0;
  double area;
  double area_error;
  int n_bnd_delaunay_edges;
  void reset() {
    n_neighbors = 0;
    t_kdtree_build = 0;
    t_kdtree_query = 0;
    t_bfs_build = 0;
    t_bfs_query = 0;
    t_sqtree_build = 0;
    t_sqtree_query = 0;
    n_sqtree_minleaf = 0;
    n_sqtree_maxleaf = 0;
    t_voronoi = 0;
    n_incomplete = 0;
    t_facets = 0;
    t_delaunay = 0;
    energy = 0;
    t_total = 0;
    count = 0;
    area = 0;
    area_error = -1;
    n_bfs_level = -1;
    n_sites = -1;
    n_triangles = -1;
    n_bnd_delaunay_edges = -1;
  }
  nlohmann::json to_json() const;
  void append_to_average(const VoronoiStatistics& stats) {
    auto avg = [this](auto& x, const auto& y) {
      x = double(x * count + y) / double(count + 1);
    };
    avg(t_kdtree_build, stats.t_kdtree_build);
    avg(t_kdtree_query, stats.t_kdtree_query);
    avg(t_bfs_build, stats.t_bfs_build);
    avg(t_bfs_query, stats.t_bfs_query);
    avg(t_sqtree_build, stats.t_sqtree_build);
    avg(t_sqtree_query, stats.t_sqtree_query);
    avg(t_voronoi, stats.t_voronoi);
    avg(t_facets, stats.t_facets);
    avg(t_delaunay, stats.t_delaunay);
    avg(t_total, stats.t_total);
    count++;
  }
};

static constexpr int kMaxClippingAttempts = 2;
static constexpr index_t kMaxSite = std::numeric_limits<index_t>::max();
enum { INSIDE = 0, OUTSIDE = 1 };
enum class VoronoiStatusCode : uint8_t {
  kIncomplete,
  kSuccess,
  kRadiusNotReached,
  kNeedPredicates
};

/// @brief Specialized vector (separate from math/vec.h) for computing Voronoi
/// diagrams.
/// @tparam T type of the elements in the vector (float or double).
/// @tparam dim dimension of the vector.
template <typename T, int dim>
struct vec;

/// @brief Specialized 4d vector for computing Voronoi diagrams.
/// @tparam T type of the elements in the vector (float or double).
template <typename T>
struct vec<T, 4> {
  T x{0};
  T y{0};
  T z{0};
  T w{0};
  vec<T, 3> xyz() const;
  template <typename R>
  vec(const R* coord, int dim) {
    x = coord[0];
    y = coord[1];
    if (dim > 2) z = coord[2];
    if (dim > 3) w = coord[3];
  }
  vec(T a, T b, T c, T d) {
    x = a;
    y = b;
    z = c;
    w = d;
  }
  vec() : x(0), y(0), z(0), w(0) {}
};

/// @brief  Specialized 3d vector for computing Voronoi diagrams.
/// @tparam T type of the elements in the vector (float or double).
template <typename T>
struct vec<T, 3> {
  T x{0};
  T y{0};
  T z{0};
  vec(T a, T b, T c) {
    x = a;
    y = b;
    z = c;
  }
  vec() : x(0), y(0), z(0) {}
  T& operator[](int d) { return *(&x + d); }
  const T& operator[](int d) const { return *(&x + d); }
};

template <typename T>
vec<T, 3> vec<T, 4>::xyz() const {
  return {x, y, z};
}

template <typename T>
vec<T, 3> cross(const vec<T, 3>& u, const vec<T, 3>& v) {
  return {u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
}

template <typename T>
vec<T, 3> operator-(const vec<T, 3>& u, const vec<T, 3>& v) {
  return {u.x - v.x, u.y - v.y, u.z - v.z};
}

template <typename T>
vec<T, 3> operator+(const vec<T, 3>& u, const vec<T, 3>& v) {
  return {u.x + v.x, u.y + v.y, u.z + v.z};
}

template <typename T, typename R>
vec<T, 3> operator*(const R& a, const vec<T, 3>& u) {
  return {a * u.x, a * u.y, a * u.z};
}

template <typename T>
T dot(const vec<T, 3>& u, const vec<T, 3>& v) {
  return u.x * v.x + u.y * v.y + u.z * v.z;
}

template <typename T>
vec<T, 4> operator-(const vec<T, 4>& u, const vec<T, 4>& v) {
  return {u.x - v.x, u.y - v.y, u.z - v.z, u.w - v.w};
}

template <typename T>
vec<T, 4> operator+(const vec<T, 4>& u, const vec<T, 4>& v) {
  return {u.x + v.x, u.y + v.y, u.z + v.z, u.w + v.w};
}

template <typename T, typename R>
vec<T, 4> operator*(const R& a, const vec<T, 4>& u) {
  return {a * u.x, a * u.y, a * u.z, a * u.w};
}

template <typename T, typename R>
vec<T, 4> operator/(const vec<T, 4>& u, const R& a) {
  return {u.x / a, u.y / a, u.z / a, u.w / a};
}

template <typename T>
T dot(const vec<T, 4>& u, const vec<T, 4>& v) {
  return u.x * v.x + u.y * v.y + u.z * v.z + u.w * v.w;
}

template <typename T>
T length(const vec<T, 3>& u) {
  return std::sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
}

template <typename T>
vec<T, 3> unit_vector(const vec<T, 3>& u) {
  const T len = length(u);
  return {u.x / len, u.y / len, u.z / len};
}

template <typename T>
T distance_squared(const vec<T, 3>& u, const vec<T, 3>& v) {
  return (u.x - v.x) * (u.x - v.x) + (u.y - v.y) * (u.y - v.y) +
         (u.z - v.z) * (u.z - v.z);
}

template <typename T>
T distance_squared(const vec<T, 4>& u, const vec<T, 4>& v) {
  return (u.x - v.x) * (u.x - v.x) + (u.y - v.y) * (u.y - v.y) +
         (u.z - v.z) * (u.z - v.z) + (u.w - v.w) * (u.w - v.w);
}

inline double det2x2(double a11, double a12, double a21, double a22) {
  return a11 * a22 - a12 * a21;
}

inline double det3x3(double a11, double a12, double a13, double a21, double a22,
                     double a23, double a31, double a32, double a33) {
  return a11 * det2x2(a22, a23, a32, a33) - a21 * det2x2(a12, a13, a32, a33) +
         a31 * det2x2(a12, a13, a22, a23);
}

inline double det4x4(double a11, double a12, double a13, double a14, double a21,
                     double a22, double a23, double a24, double a31, double a32,
                     double a33, double a34, double a41, double a42, double a43,
                     double a44) {
  double m12 = a21 * a12 - a11 * a22;
  double m13 = a31 * a12 - a11 * a32;
  double m14 = a41 * a12 - a11 * a42;
  double m23 = a31 * a22 - a21 * a32;
  double m24 = a41 * a22 - a21 * a42;
  double m34 = a41 * a32 - a31 * a42;

  double m123 = m23 * a13 - m13 * a23 + m12 * a33;
  double m124 = m24 * a13 - m14 * a23 + m12 * a43;
  double m134 = m34 * a13 - m14 * a33 + m13 * a43;
  double m234 = m34 * a23 - m24 * a33 + m23 * a43;

  return (m234 * a14 - m134 * a24 + m124 * a34 - m123 * a44);
}

using vec3 = vec<double, 3>;
using vec4 = vec<double, 4>;

/// @brief Returns the side of a point with respect to a plane.
/// @param p point represented in 4d in homogeneous coordinates.
/// @param eqn plane equation used to determine the side.
inline uint8_t plane_side(vec4 p, const vec4& eqn) {
  p = p / p.w;  // divide by homogeneous coordinate
  auto s = p.x * eqn.x + p.y * eqn.y + p.z * eqn.z + eqn.w;
  return (s < 0) ? OUTSIDE : INSIDE;
}

inline vec4 plane_equation(const vec4& ui, const vec4& uj, const coord_t& wi,
                           const coord_t& wj) {
  // https://en.wikipedia.org/wiki/Radical_axis
  // the weight is the radius squared, M1 = ui.xyz() and M2 = uj.xyz()
  coord_t d = length((uj - ui).xyz());
  if (d == 0.0) return {1e20, 1e20, 1e20, 1e20};
  ASSERT(d > 0) << d;
  coord_t d1 = (d * d + wi - wj) / (2.0 * d);
  vec4 m = ui + d1 * (uj - ui) / d;  // (uj - ui) / d is unit_vector(M2 - M1)
  vec3 n = unit_vector((ui - uj).xyz());  // normal points *into* cell i
  return {n.x, n.y, n.z, -dot(n, m.xyz())};
}

enum class NearestNeighborAlgorithm : uint8_t {
  kKdtree,
  kVoronoiBFS,
  kSphereQuadtree
};

struct VoronoiDiagramOptions {
  bool store_mesh{false};  // should the mesh be stored? (this makes the Voronoi
                           // diagram calculation slower).
  int max_kdtree_axis_dim{
      -1};  // maximum dimension used to alternate between dimensions when
            // building the kdtree (-1 means use all dimensions)
            // This is useful when building a kdtree of planar points in 3d
            // (i.e. each point has a z-coordinate of 0).
  int n_neighbors{50};  // number of nearest neighbors to precompute
  bool verbose{true};   // whether to print timing info during the calculation
  bool parallel{true};  // whether to parallelize the calculation
  Mesh* mesh{nullptr};  // destination of the mesh when store_mesh is true
  bool store_facet_data{true};
  bool store_delaunay_triangles{false};
  int bfs_max_level{3};
  bool check_closed{false};
  NearestNeighborAlgorithm neighbor_algorithm{
      NearestNeighborAlgorithm::kVoronoiBFS};
};

struct VoronoiCellProperties {
  VoronoiCellProperties() { reset(); }
  uint64_t site;         // which site/cell is this for?
  vec3 moment{0, 0, 0};  // centroid * volume
  double volume{0};      // volume
  void reset() {
    moment = {0, 0, 0};
    volume = 0;
  }
};

struct VoronoiFacetData {
  VoronoiFacetData(uint32_t a, uint32_t b, double l, vec3 c) {
    bi = a;
    bj = b;
    length = l;
    midpoint = c;
  }
  int32_t bi;
  int32_t bj;
  double length;
  vec3 midpoint;
};

struct VoronoiDiagramProperties {
  double area{0};    // total area of the Voronoi diagram
  double energy{0};  // \sum_{site \in sites} [ \int_{\vec{x} \in cell(site)}
                     // \rho \lVert \vec{x} - \vec{site}\rVert - w_{site}
                     // \mathrm{d}\vec{x} + \nu_{site} w_{site} ]
};

class VoronoiMesh : public Mesh {
 protected:
  using Mesh::Mesh;

 public:
  void allocate(size_t n_sites, int n_facets_per_site = 12) {
    size_t n_facets = n_sites * n_facets_per_site / 2;
    facets_.reserve(n_facets);
  }

  void set_save_mesh(bool x) { save_mesh_ = x; }
  void set_save_facets(bool x) { save_facets_ = x; }
  void set_save_delaunay(bool x) { save_delaunay_ = x; }

  bool save_mesh() const { return save_mesh_; }
  bool save_facets() const { return save_facets_; }
  bool save_delaunay() const { return save_delaunay_; }

  void add_facet(int32_t bi, int32_t bj, double volume, vec3 midpoint) {
#if NEW_FACETS
    if (bj >= 0 && bi > bj) return;
    facets_.emplace_back(bi, bj, volume, midpoint);
#else
    if (bi > bj) std::swap(bi, bj);
    auto it = facets_.find({bi, bj});
    if (it == facets_.end()) facets_.insert({{bi, bj}, volume});
#endif
  }

  void add_triangle(std::array<uint32_t, 3> triangle) {
    std::sort(triangle.begin(), triangle.end());
    delaunay_.insert(triangle);
  }

  const auto& facets() const { return facets_; }
  auto& facets() { return facets_; }

  const auto& delaunay() const { return delaunay_; }
  auto& delaunay() { return delaunay_; }

  void append(const VoronoiMesh& mesh) {
#if NEW_FACETS
    for (const auto& facet : mesh.facets()) facets_.push_back(facet);
#else
    for (const auto& [b, volume] : mesh.facets())
      add_facet(b.first, b.second, volume);
#endif
    n_incomplete_ += mesh.n_incomplete();
    n_boundary_facets_ += mesh.n_boundary_facets();
    boundary_area_ += mesh.boundary_area();
  }

  void append_triangles(const VoronoiMesh& mesh) {
    for (const auto& t : mesh.delaunay()) add_triangle(t);
  }

  size_t& n_incomplete() { return n_incomplete_; }
  size_t n_incomplete() const { return n_incomplete_; }

  size_t& n_boundary_facets() { return n_boundary_facets_; }
  size_t n_boundary_facets() const { return n_boundary_facets_; }

  double& boundary_area() { return boundary_area_; }
  double boundary_area() const { return boundary_area_; }

  const auto& properties() const { return properties_; }
  auto& properties() { return properties_; }

 protected:
  bool save_mesh_{false};
  bool save_facets_{false};
  bool save_delaunay_{false};
  std::vector<VoronoiCellProperties> properties_;
#if NEW_FACETS
  std::vector<VoronoiFacetData> facets_;
#else
  absl::flat_hash_map<std::pair<uint32_t, uint32_t>, double> facets_;
#endif
  absl::flat_hash_set<std::array<uint32_t, 3>> delaunay_;
  size_t n_incomplete_{0};
  size_t n_boundary_facets_{0};
  double boundary_area_{0};
};

class VoronoiDiagram : public VoronoiMesh {
 public:
  /// @brief Saves the dimension, number of sites and pointer to site
  /// coordinates.
  /// @param dim Dimension of the sites (should be 2, 3 or 4).
  /// @param sites Pointer to the sites.
  /// @param n_sites Number of sites.
  VoronoiDiagram(int dim, const coord_t* sites, uint64_t n_sites);
  VoronoiDiagram(const VoronoiDiagram&) = delete;

  /// @brief Calculates the Voronoi diagram of the saved sites restricted to
  /// some domain.
  /// @tparam Domain_t type of the domain (SquareDomain, SphereDomain,
  /// TriangulationDomain).
  /// @param domain domain to restrict the Voronoi diagram to.
  /// @param options (see above).
  template <typename Domain_t>
  void compute(const Domain_t& domain,
               VoronoiDiagramOptions options = VoronoiDiagramOptions());
  auto& status() { return status_; }
  const auto& status() const { return status_; }
  auto& weights() { return weights_; }
  const auto& weights() const { return weights_; }
  void smooth(Vertices& sites, bool on_sphere) const;

  /// @brief Calculate and return global properties like the area and energy.
  VoronoiDiagramProperties analyze() const;

  /// @brief Merge vertices with the same symbolic information (they represent
  /// the same Delaunay triangle).
  void merge();

  double max_radius() const {
    // auto props = *std::max_element(
    //     properties_.begin(), properties_.end(),
    //     [](const auto& pa, const auto& pb) { return pa.rmax > pb.rmax; });
    // return props.rmax;
    return max_radius_;
  }

  uint64_t n_sites() const { return n_sites_; }
  const coord_t* sites() const { return sites_; }
  int dim() const { return dim_; }

  void create_sqtree(int ns);
  const auto& statistics() const { return statistics_; }
  auto& statistics() { return statistics_; }
  const auto& average_statistics() const { return average_statistics_; }
  const auto& statistics_history() const { return statistics_history_; }
  void track_statistics_history(bool x) { track_statistics_history_ = x; }

 private:
  int dim_;
  const coord_t* sites_;
  uint64_t n_sites_;
  std::vector<VoronoiStatusCode> status_;
  std::vector<double> weights_;
  VoronoiNeighbors neighbors_;
  std::unique_ptr<SphereQuadtree> sqtree_;
  VoronoiStatistics statistics_;
  VoronoiStatistics average_statistics_;
  std::vector<VoronoiStatistics> statistics_history_;
  bool track_statistics_history_{true};
  double max_radius_{0};
};

/// @brief Lift the sites to 4d where the fourth coordinate = sqrt(wmax -
/// weight)
/// @param sites Vertices storing the sites (sites.dim() should be 4).
/// @param weights array of weights (sites.n() == weights.size())
void lift_sites(Vertices& sites, const std::vector<coord_t>& weights);

/// @brief Represents a Voronoi vertex in a single Voronoi cell (not the
/// entire Voronoi diagram).
using VoronoiVertex = uint8_t;

struct VoronoiCellMemoryPool {
  static constexpr int MAX_VERTICES = 256;
  static constexpr int MAX_PLANES = 256;
  static const unsigned NUM_THREADS = VORTEX_NUM_CORES;

  VoronoiCellMemoryPool() {}

  thread_memory_pool<VoronoiVertex, NUM_THREADS, MAX_VERTICES> polygon;
  thread_memory_pool<vec4, NUM_THREADS, MAX_PLANES> planes;
  thread_memory_pool<int64_t, NUM_THREADS, MAX_PLANES> bisector_to_site;
};

struct SphereDomain;

/// @brief Used to analytically calculate a Voronoi polygon on a sphere.
struct SphericalVoronoiPolygon {
  typedef VoronoiVertex Vertex_t;
  vec3 center;  // this will be the site (see SphereDomain below).

  /// @brief Calculates the polygon vertex coordinates by intersecting the
  /// line of intersection from two planes (pi, pj) with the sphere.
  /// @param pi first plane equation
  /// @param pj second plane equation
  /// @return coordinates of the intersection point in homogeneous
  /// coordinates.
  vec4 compute(const vec4& pi, const vec4& pj) const;

  /// @brief Returns the side of a point (intersection of pi, pj and sphere)
  /// with respect to the plane p.
  /// @param pi first plane used to compute the intersection point
  /// @param pj second plane used to compute the intersection point
  /// @param p plane equation
  uint8_t side(const vec4& pi, const vec4& pj, const vec4& p) const;

  /// @brief Calculate the contribution of this polygon to the cell
  /// properties.
  /// @param polygon array of Voronoi vertices.
  /// @param planes array of plane equations.
  /// @param props properties of a particular Voronoi cell to update.
  void get_properties(const device_vector<Vertex_t>& polygon,
                      const device_vector<vec4>& planes,
                      VoronoiCellProperties& props) const;
};

template <typename Domain_t>
class VoronoiPolygon;

struct SphereDomain {
  typedef SphericalVoronoiPolygon Cell_t;
  /// @brief Constructs a sphere domain with a unit radius centered at the
  /// origin.
  SphereDomain() {}

  /// @brief Copies a sphere domain.
  SphereDomain(const SphereDomain& domain)
      : radius(domain.radius),
        initialization_fraction(domain.initialization_fraction) {}

  /// @brief Sets the fraction of the radius used to initialize the square
  /// around the site
  void set_initialization_fraction(double x) { initialization_fraction = x; }

  /// @brief Initializes the Voronoi polygon (cell) to a square centered
  /// around a point (z)
  /// @param z point to center the initial square around.
  /// @param cell polygon to initialize in the calculation.
  void initialize(vec4 z, VoronoiPolygon<SphereDomain>& cell) const;
  const double radius{1.0};
  double initialization_fraction{0.7};  // use 0.7 for power diagrams

  static vec3d random_point() {
    coord_t theta = 2.0 * M_PI * irand(0, 1);
    coord_t phi = acos(2.0 * irand(0, 1) - 1.0);
    return {cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi)};
  }

  static void project(const coord_t* x, const coord_t* u, coord_t* v) {
    vec3d n(x);
    n = normalize(n);
    double u_dot_n = u[0] * n[0] + u[1] * n[1] + u[2] * n[2];
    for (int i = 0; i < 3; i++) v[i] = u[i] - u_dot_n * n[i];
  }

  template <typename V>
  static double length(const V& x, const V& y) {
    // static_assert(V::DIM == 3);
    // return std::acos(std::clamp(dot(x, y), -1.0, 1.0));
    const V c = cross(x, y);
    return std::atan2(std::sqrt(dot(c, c)), dot(x, y));
  }

  static vec3 project(const vec3& x) { return unit_vector(x); }

  double area() const { return 4 * M_PI * radius * radius; }
};

struct PlanarVoronoiPolygon {
  typedef VoronoiVertex Vertex_t;
  vec4 base{0, 0, 1, 0};

  /// @brief Calculates the polygon vertex coordinates as the intersection of
  /// the three planes: (base, pi, pi).
  /// @return coordinates of the intersection point in homogeneous
  /// coordinates.
  vec4 compute(const vec4& pi, const vec4& pj) const;

  /// @brief Returns the side of a point (intersection of pi, pj and base)
  /// with respect to the plane p.
  /// @param pi first plane used to compute the intersection point
  /// @param pj second plane used to compute the intersection point
  /// @param p plane equation
  uint8_t side(const vec4& pi, const vec4& pj, const vec4& p) const;

  /// @brief Initializes the polygon vertices and planes from the cell
  /// described by a list of points in CCW order.
  /// @param points pointer to the list of points
  /// @param n_points number of points in the cell (square = 4, triangle = 3)
  /// @param polygon output list of Voronoi polygon vertices
  /// @param planes output list of planes
  void initialize(const vec3* points, const size_t n_points,
                  device_vector<Vertex_t>& polygon,
                  device_vector<vec4>& planes);

  /// @brief Calculate the contribution of this polygon to the cell
  /// properties.
  /// @param polygon array of Voronoi vertices.
  /// @param planes array of plane equations.
  /// @param props properties of a particular Voronoi cell to update.
  void get_properties(const device_vector<Vertex_t>& polygon,
                      const device_vector<vec4>& planes,
                      VoronoiCellProperties& props) const;
};

struct SquareDomain {
  typedef PlanarVoronoiPolygon Cell_t;
  /// @brief Initializes a square domain centered at (0.5, 0.5, 0) with side
  /// lengths of 1.
  SquareDomain() {}

  /// @ brief Initializes a square domain from lower-left and upper-right
  /// corners.
  SquareDomain(vec3 pll, vec3 pur) {
    points_[0] = pll;
    points_[1] = {pur[0], pll[1], 0};
    points_[2] = pur;
    points_[3] = {pll[0], pur[1], 0};
  }

  /// @ brief Initializes a square domain from a point and two vectors defining
  /// the plane
  SquareDomain(vec3 p, vec3 u, vec3 v) {
    points_[0] = p;
    points_[1] = p + u;
    points_[2] = p + u + v;
    points_[3] = p + v;
  }

  double xlength() const { return points_[1][0] - points_[0][0]; }
  double ylength() const { return points_[2][1] - points_[0][1]; }
  double area() const {
    const vec3 n = cross(points_[1] - points_[0], points_[3] - points_[0]);
    return std::sqrt(dot(n, n));
  }

  /// @brief Copies a square domain
  SquareDomain(const SquareDomain& domain) {
    for (int i = 0; i < 4; i++) points_[i] = domain.points_[i];
  }

  vec3 random_point() const {
    const double x = points_[0][0] + rand() * xlength() / RAND_MAX;
    const double y = points_[0][1] + rand() * ylength() / RAND_MAX;
    return {x, y, 0};
  }

  static void project(const coord_t*, const coord_t* u, coord_t* v) {
    for (int i = 0; i < 3; i++) v[i] = u[i];
  }

  template <typename V>
  static double length(const V& x, const V& y) {
    return std::sqrt(dot(y - x, y - x));
  }

  static vec3 project(const vec3& x) { return x; }

  /// @brief Initializes the Voronoi polygon (cell) to the entire square
  /// domain.
  /// @param cell polygon to initialize in the calculation.
  void initialize(vec4 z, VoronoiPolygon<SquareDomain>& cell) const;
  vec3 points_[4] = {{0., 0., 0.}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
};

struct TriangulationDomain {
  typedef PlanarVoronoiPolygon Cell_t;

  /// @brief Initializes a triangulation domain by saving the vertices and
  /// triangles.
  TriangulationDomain(const coord_t* p, uint64_t np, const index_t* t,
                      uint64_t nt)
      : points(p), n_points(np), triangles(t), n_triangles(nt) {}

  /// @brief Initializes the Voronoi polygon (cell) to a particular triangle.
  /// @param cell polygon to initialize in the calculation.
  void initialize(int64_t elem,
                  VoronoiPolygon<TriangulationDomain>& cell) const;
  size_t n_elems() const { return n_triangles; }
  double area() const { return -1; }

  const coord_t* points{nullptr};
  uint64_t n_points{0};
  const index_t* triangles{nullptr};
  uint64_t n_triangles{0};
};

template <int dim>
std::shared_ptr<trees::KdTreeNd<coord_t, index_t>> get_nearest_neighbors(
    const coord_t* p, uint64_t np, const coord_t* q, uint64_t nq,
    std::vector<index_t>& knn, size_t n_neighbors,
    const VoronoiDiagramOptions& options, VoronoiStatistics& stats,
    std::shared_ptr<trees::KdTreeNd<coord_t, index_t>> ptree = nullptr);

}  // namespace vortex