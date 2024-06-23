#pragma once

#include <absl/container/flat_hash_map.h>
#include <absl/container/flat_hash_set.h>

#include <queue>

#include "array2d.h"
#include "defs.h"

namespace vortex {

class VoronoiDiagram;
class Vertices;

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

  void add(uint32_t n, uint8_t l, double d) {
    if (d > max_distance) max_distance = d;
    sites.push({n, d});
    neighbors.insert({n, l});
  }

  int n_neighbors;
  // absl::flat_hash_map<uint32_t, uint8_t> neighbors;
  std::unordered_map<uint32_t, uint8_t> neighbors;
  queue<std::pair<uint32_t, double>> sites;
  size_t n_avg{0};
  double max_distance{0};
};

class VoronoiNeighbors {
 public:
  VoronoiNeighbors(const VoronoiDiagram& voronoi, const coord_t* points);

  void knearest(uint32_t p, NearestNeighborsWorkspace& search) const;

 private:
  void build();

  const VoronoiDiagram& voronoi_;
  const coord_t* points_;
  array2d<uint32_t> ring_;
};

}  // namespace vortex