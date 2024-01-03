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

template <typename T, int dim>
struct vec;

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

inline uint8_t plane_side(vec4 p, const vec4& eqn) {
  p = p / p.w;
  auto s = p.x * eqn.x + p.y * eqn.y + p.z * eqn.z + eqn.w;
  return (s < 0) ? OUTSIDE : INSIDE;
}

struct VoronoiDiagramOptions {
  bool store_mesh{false};
  int max_kdtree_axis_dim{-1};
  int n_neighbors{50};
  bool verbose{true};
  bool parallel{true};
  bool interleave_neighbors{false};
  bool allow_reattempt{true};
  Mesh* mesh{nullptr};
};

struct VoronoiCellProperties {
  VoronoiCellProperties() { reset(); }
  uint64_t site;
  uint64_t elem;
  vec3 moment{0, 0, 0};
  double mass{0};
  void reset() {
    moment = {0, 0, 0};
    mass = 0;
  }
};

struct VoronoiDiagramProperties {
  double area{0};
  double energy{0};
};

class VoronoiDiagram : public Mesh {
 public:
  VoronoiDiagram(int dim, const coord_t* sites, uint64_t n_sites)
      : Mesh(3), dim_(dim), sites_(sites), n_sites_(n_sites) {}

  template <typename Domain_t>
  void compute(const Domain_t& domain,
               VoronoiDiagramOptions options = VoronoiDiagramOptions());

  const auto& properties() const { return properties_; }
  auto& status() { return status_; }
  const auto& status() const { return status_; }
  auto& weights() { return weights_; }
  const auto& weights() const { return weights_; }
  void smooth(Vertices& sites) const;
  VoronoiDiagramProperties analyze() const;
  void merge();

 private:
  int dim_;
  const coord_t* sites_;
  uint64_t n_sites_;
  std::vector<VoronoiCellProperties> properties_;
  std::vector<VoronoiStatusCode> status_;
  std::vector<double> weights_;
};

void lift_sites(Vertices& sites, const std::vector<coord_t>& weights);

struct VoronoiVertex {
  uint8_t bl;  // left bisector
  uint8_t br;  // right bisector
};

struct SphereDomain;
struct SphericalVoronoiPolygon {
  typedef VoronoiVertex Vertex_t;
  vec3 center;

  vec4 compute(const vec4& pi, const vec4& pj) const;
  uint8_t side(const vec4& pi, const vec4& pj, const vec4& p) const {
    return plane_side(compute(pi, pj), p);
  }
  vec4 plane_equation(const vec4& ui, const vec4& uj, const coord_t& wi,
                      const coord_t& wj);
  void get_properties(const pool<Vertex_t>& polygon, const pool<vec4>& planes,
                      VoronoiCellProperties& props) const;
};

template <typename Domain_t>
class VoronoiPolygon;

struct SphereDomain {
  typedef SphericalVoronoiPolygon Cell_t;
  SphereDomain(double r) : radius(r) {}
  SphereDomain(const SphereDomain& domain) : radius(domain.radius) {}
  void initialize(vec4 z, VoronoiPolygon<SphereDomain>& cell) const;
  double radius{1.0};
};

struct PlanarVoronoiPolygon {
  typedef VoronoiVertex Vertex_t;
  vec4 base{0, 0, 1, 0};

  vec4 compute(const vec4& pi, const vec4& pj) const;
  uint8_t side(const vec4& pi, const vec4& pj, const vec4& p) const;
  void initialize(const vec3* points, const size_t n_points,
                  pool<Vertex_t>& polygon, pool<vec4>& planes);
  vec4 plane_equation(const vec4& ui, const vec4& uj, const coord_t& wi,
                      const coord_t& wj);
  void get_properties(const pool<Vertex_t>& polygon, const pool<vec4>& planes,
                      VoronoiCellProperties& props) const;
};

struct SquareDomain {
  typedef PlanarVoronoiPolygon Cell_t;
  SquareDomain() {}
  SquareDomain(const SquareDomain& domain) {}
  void initialize(vec4 z, VoronoiPolygon<SquareDomain>& cell) const;
  const vec3 points_[4] = {{0., 0., 0.}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
};

struct TriangulationDomain {
  typedef PlanarVoronoiPolygon Cell_t;

  TriangulationDomain(const coord_t* p, uint64_t np, const index_t* t,
                      uint64_t nt)
      : points(p), n_points(np), triangles(t), n_triangles(nt) {}
  void initialize(int64_t elem,
                  VoronoiPolygon<TriangulationDomain>& cell) const;
  size_t n_elems() const { return n_triangles; }

  const coord_t* points{nullptr};
  uint64_t n_points{0};
  const index_t* triangles{nullptr};
  uint64_t n_triangles{0};
};

}  // namespace vortex