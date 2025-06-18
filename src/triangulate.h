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
#pragma once

#include <vector>

#include "halfedges.h"

namespace vortex {

class Mesh;
class Vertices;
template <typename T>
class Topology;

/// @brief Triangulates the oceans by inserting points from the coast and
/// recovering the edges stored in the coast.lines(). [WARNING] This is a
/// work-in-progress and should not be used yet.
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

/// @brief Triangulates a general polygon using ear clipping.
class EarClipper {
 private:
  /// @brief Helper structure to represent a node in the linked list defining
  /// the current polygon.
  struct Node {
    uint8_t next;
    uint8_t prev;
    uint8_t indx;
  };

 public:
  /// @brief Triangulates the polygon.
  /// @param points array of points to triangulate.
  /// @param normal Normal vector to use to determine the orientation of
  /// triangles.
  /// @return Whether the polygon was determine to be convex (true) or concave
  /// (false) during the triangulation procedure.
  bool triangulate(const std::vector<vec3d>& points, const vec3d& normal);

  /// @brief Returns the number of triangles in the triangulation of the
  /// polygon.
  size_t n_triangles() const { return triangles_.size() / 3; }

  /// @brief Returns a pointer to the first index of triangle k in the
  /// triangulation.
  /// @param k triangle index.
  const index_t* triangle(size_t k) const { return triangles_.data() + 3 * k; }

  /// @brief Returns whether the edge opposite vertex j of triangle k is on the
  /// boundary of the polygon.
  /// @param k triangle index
  /// @param j vertex index (edge is opposite the vertex).
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

  /// @brief Triangulates polygons m through n (not including n).
  /// @param type How the tangent space should be computed (see enum above).
  void triangulate(TangentSpaceType type, size_t m, size_t n);

  /// @brief Returns the total number of triangles in the triangulation of all
  /// polygons.
  size_t n() const { return triangles_.size() / 3; }

  /// @brief Returns a pointer to the first index of triangle k in the
  /// triangulation.
  /// @param k triangle index.
  const index_t* triangle(size_t k) const { return triangles_.data() + 3 * k; }

  /// @brief Returns the index of the polygon for triangle k.
  /// @param k triangle index.
  int group(size_t k) const { return group_[k]; }

  /// @brief Returns whether the edge opposite vertex j of triangle k is on the
  /// boundary of the polygon.
  /// @param k triangle index
  /// @param j vertex index (edge is opposite the vertex).
  bool edge(size_t k, int j) const { return edge_[3 * k + j]; }

 private:
  const Vertices& vertices_;
  const Topology<Polygon>& polygons_;
  std::vector<index_t> triangles_;
  std::vector<bool> edge_;
  std::vector<int> group_;
};

}  // namespace vortex