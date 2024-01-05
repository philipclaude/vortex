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
#include "mesh.h"

#include <stlext.h>
#include <trees/kdtree.h>

#include <set>

#include "elements.h"
#include "math/vec.hpp"

namespace vortex {

template <typename T>
void get_element_edges(const Topology<T>& elems, int k,
                       std::vector<Edge>& edges) {
  for (int j = 0; j < T::n_edges; j++) {
    auto e0 = elems(k, T::edges[2 * j]);
    auto e1 = elems(k, T::edges[2 * j + 1]);

    if (e0 > e1) std::swap(e0, e1);
    edges.push_back({e0, e1});
  }
}

template <>
void get_element_edges<Polygon>(const Topology<Polygon>& elems, int k,
                                std::vector<Edge>& edges) {
  auto n_vertices = elems.length(k);
  for (int j = 0; j < n_vertices; j++) {
    auto e0 = elems(k, j);
    auto e1 = elems(k, (j + 1) % n_vertices);

    if (e0 > e1) std::swap(e0, e1);
    edges.push_back({e0, e1});
  }
}

template <typename T>
void Topology<T>::append_edges(std::vector<Edge>& edges) const {
  std::set<Edge> E;
  std::vector<Edge> element_edges;

  for (size_t k = 0; k < n(); k++) {
    element_edges.clear();
    get_element_edges(*this, k, element_edges);
    for (const auto& e : element_edges) {
      if (E.find(e) == E.end()) {
        E.insert(e);
        edges.push_back(e);
      }
    }
  }
}

template <>
void Topology<Triangle>::flip_orientation() {
  for (size_t k = 0; k < n(); k++) {
    index_t t1 = (*this)(k, 1);
    index_t t2 = (*this)(k, 2);
    (*this)(k, 2) = t1;
    (*this)(k, 1) = t2;
  }
}

void Mesh::get_edges(std::vector<Edge>& edges) const {
  edges.clear();
  triangles_.append_edges(edges);
  polygons_.append_edges(edges);
}

template <>
const Topology<Triangle>& Mesh::get<Triangle>() const {
  return triangles_;
}

template <>
Topology<Triangle>& Mesh::get<Triangle>() {
  return triangles_;
}

template <>
const Topology<Polygon>& Mesh::get<Polygon>() const {
  return polygons_;
}

template <>
Topology<Polygon>& Mesh::get<Polygon>() {
  return polygons_;
}

template <>
const Topology<Quad>& Mesh::get<Quad>() const {
  return quads_;
}

template <>
Topology<Quad>& Mesh::get<Quad>() {
  return quads_;
}

void Vertices::print() const {
  for (size_t k = 0; k < n(); k++) {
    std::cout << fmt::format("v[{}] = (", k);
    for (int d = 0; d < dim(); d++) {
      std::cout << (*this)[k][d];
      if (d + 1 < dim())
        std::cout << ", ";
      else
        std::cout << ") ";
    }
    std::cout << std::endl;
  }
}

void Mesh::merge(double tol) {
  ASSERT(vertices_.dim() == 3) << "kdtree dim set to 3 for now";
  trees::KdTree<3, coord_t, index_t> tree(vertices_[0], vertices_.n());

  std::vector<bool> visited(vertices_.n(), false);
  std::vector<index_t> vmap(vertices_.n());

  struct MergeWorkspace {
    std::vector<index_t> neighbors;
    std::vector<coord_t> distance;
    std::vector<index_t> mergelist;
  };
  size_t n_threads = std::thread::hardware_concurrency();
  std::vector<MergeWorkspace> workspace(n_threads);

  std::mutex lock;
  std::parafor_i(
      0, vertices_.n(),
      [&workspace, this, tol, &visited, &vmap, &tree, &lock](int tid, int k) {
        if (visited[k]) return;
        const auto* pk = vertices_[k];
        int n_neighbors = 2;
        auto& neighbors = workspace[tid].neighbors;
        auto& distance = workspace[tid].distance;
        auto& mergelist = workspace[tid].mergelist;
        while (true) {
          neighbors.resize(n_neighbors);
          distance.resize(n_neighbors);
          trees::NearestNeighborSearch<index_t, coord_t> search(
              n_neighbors, neighbors.data(), distance.data());
          tree.knearest(pk, search);
          mergelist.clear();
          mergelist.push_back(k);
          for (size_t j = 1; j < n_neighbors; j++) {
            const auto* pj = vertices_[neighbors[j]];
            double d_squared = 0;
            for (int d = 0; d < vertices_.dim(); d++)
              d_squared += (pj[d] - pk[d]) * (pj[d] - pk[d]);
            if (std::sqrt(d_squared) > tol) break;
            mergelist.push_back(neighbors[j]);
          }
          if (mergelist.size() < n_neighbors) break;
          n_neighbors += 5;
        }
        ASSERT(mergelist.size() < n_neighbors)
            << n_neighbors << ", " << mergelist.size();
        if (mergelist.size() == 1) {
          vmap[k] = k;
          return;
        }
        std::sort(mergelist.begin(), mergelist.end());
        lock.lock();
        for (size_t i = 0; i < mergelist.size(); i++) {
          vmap[mergelist[i]] = mergelist[0];
          visited[mergelist[i]] = true;
        }
        lock.unlock();
      });

  // map the element indices
  for (size_t k = 0; k < lines_.n(); k++)
    for (int j = 0; j < 2; j++) lines_[k][j] = vmap[lines_[k][j]];

  for (size_t k = 0; k < polygons_.n(); k++)
    for (int j = 0; j < polygons_.length(k); j++)
      polygons_[k][j] = vmap[polygons_[k][j]];

  for (size_t k = 0; k < triangles_.n(); k++)
    for (int j = 0; j < 3; j++) triangles_[k][j] = vmap[triangles_[k][j]];

  for (size_t k = 0; k < quads_.n(); k++)
    for (int j = 0; j < 4; j++) quads_[k][j] = vmap[quads_[k][j]];
}

}  // namespace vortex
