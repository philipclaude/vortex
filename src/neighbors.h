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

#include <algorithm>
#include <cmath>
#include <memory>
#include <unordered_map>
#include <vector>

#include "array2d.h"
#include "defs.h"
#include "mesh.h"

namespace vortex {

class VoronoiDiagram;

template <typename T>
class queue {
 public:
  void reserve(size_t n) { data_.reserve(n); }

  void clear() {
    data_.clear();
    front_ = 0;
  }
  void push(const T& x) { data_.push_back(x); }

  void pop() { front_++; }
  const T& front() const { return data_[front_]; }

  bool empty() const { return front_ == data_.size(); }
  auto& data() { return data_; }
  const auto& data() const { return data_; }
  size_t size() const { return data_.size(); }

 private:
  size_t front_{0};
  std::vector<T> data_;
};

struct NearestNeighborsWorkspace {
  NearestNeighborsWorkspace(int k, int nmax = 200) : n_neighbors(k) {
    neighbors.reserve(nmax);
  }

  void reset() {
    neighbors.clear();
    sites.clear();
    max_distance = 0;
  }

  void sort() {
    std::sort(sites.data().begin(), sites.data().end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });
  }

  uint32_t next() {
    uint32_t n = sites.front().first;
    sites.pop();
    return n;
  }

  void add(uint32_t n, uint8_t l, double d) {
    if (d > max_distance) max_distance = d;
    sites.push({n, d});
    neighbors.insert({n, l});
  }

  int n_neighbors;
  std::unordered_map<uint32_t, uint8_t> neighbors;
  queue<std::pair<uint32_t, double>> sites;
  size_t total_neighbors{0};
  double max_distance{0};
  int max_level{2};
};

class VoronoiNeighbors {
 public:
  VoronoiNeighbors(const VoronoiDiagram& voronoi, const coord_t* points,
                   int dim);

  void build();
  void knearest(uint32_t p, NearestNeighborsWorkspace& search) const;

 private:
  const int dim_;
  const VoronoiDiagram& voronoi_;
  const coord_t* points_;
  array2d<uint32_t> ring_;
};

struct SphereQuadtreeWorkspace {
  SphereQuadtreeWorkspace(int n) : n_neighbors(n) {}
  void reset() {
    neighbors.clear();
    max_distance = 0;
  }

  void sort() {
    if (neighbors.size() < n_neighbors) {
      std::sort(
          neighbors.begin(), neighbors.end(),
          [](const auto& a, const auto& b) { return a.second < b.second; });
      n_neighbors = neighbors.size();
    } else {
      std::partial_sort(
          neighbors.begin(), neighbors.begin() + n_neighbors, neighbors.end(),
          [](const auto& a, const auto& b) { return a.second < b.second; });
    }
  }

  void add(uint32_t n, double d) {
    // if (neighbors.size() > 50 && d > max_distance) return;
    // if (d > max_distance) max_distance = d;
    neighbors.push_back({n, d});
  }

  double max_distance;
  std::vector<std::pair<uint32_t, double>> neighbors;
  int n_neighbors;
};

class SphereQuadtree {
 public:
  SphereQuadtree(const coord_t* points, size_t n_points, int dim, int ns);

  void setup();
  void build();
  void knearest(uint32_t p, SphereNeighborsWorkspace& search) const;

  size_t n_triangles() const { return search_triangles_.size(); }
  const auto& mesh() const { return mesh_; }

 private:
  struct Subdivision : public Mesh {
    Subdivision(int ns);
    array2d<uint32_t> children;
    int n_levels;
  };

  // TODO make a lookup table for this
  // 8 * \sum_{i = 0}^n 4^i
  size_t t_first(int ns) {
    return Octahedron::n_faces * (std::pow(4, ns) - 1) / 3;
  }
  size_t t_last(int ns) {
    return Octahedron::n_faces * (std::pow(4, ns + 1) - 1) / 3;
  }

  const coord_t* points_;
  size_t n_points_;
  int dim_;
  std::vector<uint32_t> point2triangle_;
#if 1
  std::unordered_map<uint32_t, std::vector<uint32_t>> triangle2points_;
  std::unordered_map<uint32_t, std::vector<uint32_t>> search_triangles_;
#else
  absl::flat_hash_map<uint32_t, std::vector<uint32_t>> triangle2points_;
  absl::flat_hash_map<uint32_t, std::vector<uint32_t>> search_triangles_;
#endif
  Subdivision mesh_;
};

}  // namespace vortex