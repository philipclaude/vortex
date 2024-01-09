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

#include <cassert>
#include <cstdint>

#include "defs.h"
#include "log.h"
#include "mesh.h"

namespace vortex {

static constexpr index_t kMaxSite = std::numeric_limits<index_t>::max();
enum { INSIDE = 0, OUTSIDE = 1 };
using vint_t = uint16_t;
enum class VoronoiStatusCode : uint8_t {
  kIncomplete,
  kSuccess,
  kRadiusNotReached,
  kNeedPredicates
};

/// @brief Helper class similar to a vector but accepts a buffer for storing
/// data. On the CPU, this would be the memory allocated by std::vector, but
/// on the GPU this could be a pointer to shared memory.
/// @tparam T type of the data to store.
template <typename T>
class pool {
 public:
  pool(T* data, int n) : data_(data), capacity_(n) {}

  T& operator[](size_t k) { return data_[k]; }
  const T& operator[](size_t k) const { return data_[k]; }
  void push_back(const T& x) {
    assert(index_ < capacity_);
    data_[index_++] = x;
  }

  size_t capacity() const { return capacity_; }
  size_t size() const { return index_; }
  void set_size(size_t n) {
    ASSERT(n <= capacity_) << fmt::format("n = {}, capacity_ = {}", n,
                                          capacity_);
    index_ = n;
  }
  bool empty() const { return index_ == 0; }
  void resize(size_t n) {
    assert(n < capacity_);
    index_ = n;
  }
  void clear() { index_ = 0; }
  T& back() { return data_[index_ - 1]; }
  T*& data() { return data_; }

#define SWAP(X, Y) \
  {                \
    auto t = Y;    \
    Y = X;         \
    X = t;         \
  }
  void swap(pool<T>& v) {
    SWAP(v.data_, data_);
    SWAP(v.capacity_, capacity_);
    SWAP(v.index_, index_);
  }
#undef SWAP

 private:
  T* data_{nullptr};
  size_t index_{0};
  size_t capacity_{0};
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
  coord_t d1 = (d * d + wi - wj) / (2.0 * d);
  vec4 m = ui + d1 * (uj - ui) / d;  // (uj - ui) / d is unit_vector(M2 - M1)
  vec3 n = unit_vector((ui - uj).xyz());  // normal points *into* cell i
  return {n.x, n.y, n.z, -dot(n, m.xyz())};
}

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
  bool interleave_neighbors{false};  // NOT IMPLEMENTED
  bool allow_reattempt{true};        // NOT IMPLEMENTED
  Mesh* mesh{nullptr};  // destination of the mesh when store_mesh is true
};

struct VoronoiCellProperties {
  VoronoiCellProperties() { reset(); }
  uint64_t site;         // which site/cell is this for?
  vec3 moment{0, 0, 0};  // centroid * mass
  double mass{0};        // mass
  void reset() {
    moment = {0, 0, 0};
    mass = 0;
  }
};

struct VoronoiDiagramProperties {
  double area{0};    // total area of the Voronoi diagram
  double energy{0};  // \sum_{site \in sites} [ \int_{\vec{x} \in cell(site)}
                     // \rho \lVert \vec{x} - \vec{site}\rVert - w_{site}
                     // \mathrm{d}\vec{x} + \nu_{site} w_{site} ]
};

class VoronoiDiagram : public Mesh {
 public:
  /// @brief Saves the dimension, number of sites and pointer to site
  /// coordinates.
  /// @param dim Dimension of the sites (should be 2, 3 or 4).
  /// @param sites Pointer to the sites.
  /// @param n_sites Number of sites.
  VoronoiDiagram(int dim, const coord_t* sites, uint64_t n_sites)
      : Mesh(3), dim_(dim), sites_(sites), n_sites_(n_sites) {}
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

  const auto& properties() const { return properties_; }
  auto& status() { return status_; }
  const auto& status() const { return status_; }
  auto& weights() { return weights_; }
  const auto& weights() const { return weights_; }
  void smooth(Vertices& sites) const;

  /// @brief Calculate and return global properties like the area and energy.
  VoronoiDiagramProperties analyze() const;

  /// @brief Merge vertices with the same symbolic information (they represent
  /// the same Delaunay triangle).
  void merge();

 private:
  int dim_;
  const coord_t* sites_;
  uint64_t n_sites_;
  std::vector<VoronoiCellProperties> properties_;
  std::vector<VoronoiStatusCode> status_;
  std::vector<double> weights_;
};

/// @brief Lift the sites to 4d where the fourth coordinate = sqrt(wmax -
/// weight)
/// @param sites Vertices storing the sites (sites.dim() should be 4).
/// @param weights array of weights (sites.n() == weights.size())
void lift_sites(Vertices& sites, const std::vector<coord_t>& weights);

/// @brief Represents a Voronoi vertex in a single Voronoi cell (not the entire
/// Voronoi diagram).
struct VoronoiVertex {
  uint8_t bl;  // left bisector
  uint8_t br;  // right bisector
};

struct SphereDomain;

/// @brief Used to analytically calculate a Voronoi polygon on a sphere.
struct SphericalVoronoiPolygon {
  typedef VoronoiVertex Vertex_t;
  vec3 center;  // this will be the site (see SphereDomain below).

  /// @brief Calculates the polygon vertex coordinates by intersecting the line
  /// of intersection from two planes (pi, pj) with the sphere.
  /// @param pi first plane equation
  /// @param pj second plane equation
  /// @return coordinates of the intersection point in homogeneous coordinates.
  vec4 compute(const vec4& pi, const vec4& pj) const;

  /// @brief Returns the side of a point (intersection of pi, pj and sphere)
  /// with respect to the plane p.
  /// @param pi first plane used to compute the intersection point
  /// @param pj second plane used to compute the intersection point
  /// @param p plane equation
  uint8_t side(const vec4& pi, const vec4& pj, const vec4& p) const {
    return plane_side(compute(pi, pj), p);
  }

  /// @brief Calculate the contribution of this polygon to the cell properties.
  /// @param polygon array of Voronoi vertices.
  /// @param planes array of plane equations.
  /// @param props properties of a particular Voronoi cell to update.
  void get_properties(const pool<Vertex_t>& polygon, const pool<vec4>& planes,
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
  SphereDomain(const SphereDomain& domain) : radius(domain.radius) {}

  /// @brief Initializes the Voronoi polygon (cell) to a square centered around
  /// a point (z)
  /// @param z point to center the initial square around.
  /// @param cell polygon to initialize in the calculation.
  void initialize(vec4 z, VoronoiPolygon<SphereDomain>& cell) const;
  const double radius{1.0};
};

struct PlanarVoronoiPolygon {
  typedef VoronoiVertex Vertex_t;
  vec4 base{0, 0, 1, 0};

  /// @brief Calculates the polygon vertex coordinates as the intersection of
  /// the three planes: (base, pi, pi).
  /// @return coordinates of the intersection point in homogeneous coordinates.
  vec4 compute(const vec4& pi, const vec4& pj) const;

  /// @brief Returns the side of a point (intersection of pi, pj and base)
  /// with respect to the plane p.
  /// @param pi first plane used to compute the intersection point
  /// @param pj second plane used to compute the intersection point
  /// @param p plane equation
  uint8_t side(const vec4& pi, const vec4& pj, const vec4& p) const;

  /// @brief Initializes the polygon vertices and planes from the cell described
  /// by a list of points in CCW order.
  /// @param points pointer to the list of points
  /// @param n_points number of points in the cell (square = 4, triangle = 3)
  /// @param polygon output list of Voronoi polygon vertices
  /// @param planes output list of planes
  void initialize(const vec3* points, const size_t n_points,
                  pool<Vertex_t>& polygon, pool<vec4>& planes);

  /// @brief Calculate the contribution of this polygon to the cell properties.
  /// @param polygon array of Voronoi vertices.
  /// @param planes array of plane equations.
  /// @param props properties of a particular Voronoi cell to update.
  void get_properties(const pool<Vertex_t>& polygon, const pool<vec4>& planes,
                      VoronoiCellProperties& props) const;
};

struct SquareDomain {
  typedef PlanarVoronoiPolygon Cell_t;
  /// @brief Initializes a square domain centered at (0.5, 0.5, 0) with side
  /// lengths of 1.
  SquareDomain() {}

  /// @brief Copies a square domain
  SquareDomain(const SquareDomain& domain) {}

  /// @brief Initializes the Voronoi polygon (cell) to the entire square domain.
  /// @param cell polygon to initialize in the calculation.
  void initialize(vec4 z, VoronoiPolygon<SquareDomain>& cell) const;
  const vec3 points_[4] = {{0., 0., 0.}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
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

  const coord_t* points{nullptr};
  uint64_t n_points{0};
  const index_t* triangles{nullptr};
  uint64_t n_triangles{0};
};

}  // namespace vortex