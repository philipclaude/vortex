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

#include <vector>

#include "halfedges.h"

namespace vortex {

class Mesh;
class Vertices;
template <typename T>
class Topology;

class OceanTriangulator {
 public:
  struct Options {};
  OceanTriangulator(Mesh& mesh, const Mesh& coast);
  void triangulate();

 private:
  void insert_points();
  void recover_edges();

  Mesh& mesh_;
  const Mesh& coast_;
  HalfMesh hmesh_;
};

class EarClipper {
 private:
  struct Node {
    uint8_t next;
    uint8_t prev;
    uint8_t indx;
  };

 public:
  void triangulate(const std::vector<vec3d>& points, const vec3d& normal);

  size_t n_triangles() const { return triangles_.size() / 3; }
  const index_t* triangle(size_t k) const { return triangles_.data() + 3 * k; }
  bool boundary(size_t k, int j) const { return boundary_[3 * k + j]; }

 private:
  std::vector<vec3d> points_;
  std::vector<index_t> triangles_;
  std::vector<bool> boundary_;
  std::vector<Node> nodes_;
};

enum class TangentSpaceType : uint8_t {
  kPlanar = 0,
  kSphere = 1,
  kGeneral = 2
};

class PolygonTriangulation {
 public:
  PolygonTriangulation(const Vertices& vertices,
                       const Topology<Polygon>& polygons);

  void triangulate(TangentSpaceType type, size_t m, size_t n);
  size_t n() const { return triangles_.size() / 3; }
  const index_t* triangle(size_t k) const { return triangles_.data() + 3 * k; }
  int group(size_t k) const { return group_[k]; }
  bool edge(size_t k, int j) const { return edge_[3 * k + j]; }

 private:
  const Vertices& vertices_;
  const Topology<Polygon>& polygons_;
  std::vector<index_t> triangles_;
  std::vector<bool> edge_;
  std::vector<int> group_;
};

}  // namespace vortex