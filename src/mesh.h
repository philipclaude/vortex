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

#include <array>
#include <unordered_map>

#include "array2d.h"
#include "defs.h"
#include "elements.h"
#include "field.h"

namespace vortex {

class TopologyBase : public array2d<index_t> {
 public:
  using array2d<index_t>::n;
  using array2d<index_t>::length;

 protected:
  TopologyBase(int stride) : array2d<index_t>(stride) {}

 public:
  void reserve(int64_t m) {
    array2d<index_t>::reserve(m);
    group_.reserve(m);
  }

  template <typename R>
  void add(const R* x, int m = -1) {
    (m < 0) ? array2d<index_t>::template add<R>(x)
            : array2d<index_t>::template add<R>(x, m);
    group_.push_back(-1);
  }

  int group(index_t k) const {
    ASSERT(k < n());
    return group_[k];
  }

  void set_group(index_t k, int32_t value) {
    ASSERT(k < n());
    group_[k] = value;
  }

  const auto& groups() const { return group_; }
  auto& groups() { return group_; }

 protected:
  std::vector<int32_t> group_;
};

template <typename T>
class Topology : public TopologyBase {
 public:
  using TopologyBase::length;
  using TopologyBase::n;

  Topology() : TopologyBase(T::n_vertices) {}

  void append_edges(std::vector<Edge>& edges) const;
  void flip_orientation();
};

class Entity;
class Vertices : public array2d<coord_t> {
 public:
  static constexpr int max_dim = 4;
  Vertices(int dim) : array2d<coord_t>(dim), param_(dim - 1) {
    ASSERT(dim <= max_dim);
  }

  int dim() const { return array2d<coord_t>::stride(); }
  void set_dim(int dim) { array2d<coord_t>::set_stride(dim); }

  template <typename R>
  void add(const R* x, int32_t id = -1) {
    array2d<coord_t>::template add<R>(x);
    group_.push_back(id);
    entity_.push_back(nullptr);
    std::array<double, max_dim - 1> u;
    std::fill(u.begin(), u.end(), 0);
    param_.add(u.data());
  }

  int32_t group(size_t k) const {
    ASSERT(k < n());
    return group_[k];
  }

  void set_group(size_t k, int value) {
    ASSERT(k < n());
    group_[k] = value;
  }

  void set_entity(size_t k, Entity* entity) {
    ASSERT(k < n());
    entity_[k] = entity;
  }

  void set_param(size_t k, const coord_t* u, int nu) {
    ASSERT(k < n());
    ASSERT(k < param_.n());
    ASSERT(nu < dim());
    for (int d = 0; d < nu; d++) param_[k][d] = u[d];
  }

  const std::vector<int>& group() const { return group_; }
  std::vector<int32_t>& group() { return group_; }

  Entity* entity(size_t k) const {
    ASSERT(k < n());
    return entity_[k];
  }

  void print() const;
  const auto& groups() const { return group_; }
  auto& params() { return param_; }

  void allocate(size_t n) {
    entity_.resize(n);
    group_.resize(n);
  }

  void copy(Vertices& dst) const {
    array2d<coord_t>::copy(dst);
    dst.allocate(n());
    param_.copy(dst.params());
    for (size_t k = 0; k < n(); k++) {
      dst.set_entity(k, entity_[k]);
      dst.set_group(k, group_[k]);
      dst.set_param(k, param_[k], param_.stride());
    }
  }

 private:
  std::vector<int32_t> group_;
  std::vector<vortex::Entity*> entity_;
  array2d<double> param_;
};

class Mesh {
 public:
  Mesh(int dim) : vertices_(dim) {}

  Topology<Line>& lines() { return lines_; }
  const Topology<Line>& lines() const { return lines_; }

  Topology<Triangle>& triangles() { return triangles_; }
  const Topology<Triangle>& triangles() const { return triangles_; }

  Topology<Quad>& quads() { return quads_; }
  const Topology<Quad>& quads() const { return quads_; }

  Topology<Polygon>& polygons() { return polygons_; }
  const Topology<Polygon>& polygons() const { return polygons_; }

  Vertices& vertices() { return vertices_; }
  const Vertices& vertices() const { return vertices_; }

  void get_edges(std::vector<Edge>& edges) const;

  template <typename T>
  const Topology<T>& get() const;

  template <typename T>
  Topology<T>& get();

  const FieldLibrary& fields() const { return fields_; }
  FieldLibrary& fields() { return fields_; }

  int get_surface_connected_components(std::vector<int>& components) const;
  void merge(double tol = 1e-10);

 protected:
  Vertices vertices_;
  Topology<Line> lines_;
  Topology<Triangle> triangles_;
  Topology<Quad> quads_;
  Topology<Polygon> polygons_;

  FieldLibrary fields_;
};

}  // namespace vortex
