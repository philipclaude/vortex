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

#include <stack>
#include <vector>

#include "log.h"
#include "mesh.h"

namespace vortex {

using half_t = int64_t;
static constexpr half_t halfnull_t = -1;

class HalfMesh;

class HalfEdge;
class HalfNode {
 public:
  HalfNode(HalfMesh& m, const double* x, int64_t index);
  half_t index() const { return index_; }
  half_t edge() const { return edge_; }
  HalfEdge& get_edge();
  const HalfEdge& get_edge() const;
  void set_index(int64_t id) { index_ = id; }
  void set_edge(const HalfEdge& e);
  const double* point() const;
  double* point();
  void deactivate() { index_ = halfnull_t; }
  bool active() const { return index_ != halfnull_t; }

  template <typename T>
  void get_onering(std::vector<T*>& ring) const;

 private:
  HalfMesh& mesh_;             // reference to the HalfMesh container
  int64_t index_{halfnull_t};  // index of this node in the HalfMesh::nodes_ as
                               // well as the coordinates in HalfMesh::vertices_
  half_t edge_{halfnull_t};    // index of an edge in HalfMesh::edges_
};

class HalfFace;
class HalfEdge {
 public:
  HalfEdge(HalfMesh& m, int64_t id) : mesh_(m), index_(id) {}
  void deactivate() { index_ = halfnull_t; }
  bool active() const { return index_ != halfnull_t; }
  bool boundary() const;
  HalfNode& get_node();
  const HalfNode& get_node() const;

  HalfEdge& get_twin();
  const HalfEdge& get_twin() const;
  HalfEdge& get_prev();
  const HalfEdge& get_prev() const;

  HalfNode& get_triangle_left_node();
  HalfNode& get_triangle_right_node();
  const HalfNode& get_triangle_left_node() const;
  const HalfNode& get_triangle_right_node() const;
  HalfFace& get_face();
  const HalfFace& get_face() const;
  HalfEdge& get_next();
  const HalfEdge& get_next() const;

  void set_index(int64_t id) { index_ = id; }
  void set_face(const HalfFace& f);
  void set_next(const HalfEdge& e) { next_ = e.index(); }
  void set_prev(const HalfEdge& e) { prev_ = e.index(); }
  void set_node(const HalfNode& n) { node_ = n.index(); }
  void set_twin(const HalfEdge& e) { twin_ = e.index(); }
  void make_pair(HalfEdge& e) {
    set_twin(e);
    e.set_twin(*this);
  }

  half_t index() const { return index_; }
  half_t twin() const { return twin_; }
  half_t node() const { return node_; }
  half_t face() const { return face_; }
  half_t prev() const { return prev_; }
  half_t next() const { return next_; }

 private:
  HalfMesh& mesh_;            // reference to the HalfMesh container
  half_t index_{halfnull_t};  // index of this edge in HalfMesh::edges_
  half_t twin_{halfnull_t};   // index of the twin in HalfMesh::edges_
  half_t node_{halfnull_t};   // index of the origin node in HalfMesh::nodes_
  half_t face_{halfnull_t};   // index of the left face in HalfMesh::faces_
  half_t prev_{halfnull_t};  // index of the previous edge (in CCW order) around
                             // the face in HalfMesh::edges_
  half_t next_{halfnull_t};  // index of the next edge (in CCW order) around the
                             // face in HalfMesh::edges_
};

class HalfFace {
 public:
  HalfFace(HalfMesh& mesh, int64_t index, uint8_t n)
      : mesh_(mesh), index_(index), n_(n) {}
  void deactivate() { index_ = halfnull_t; }
  bool active() const { return index_ != halfnull_t; }

  half_t index() const { return index_; }
  half_t edge() const { return edge_; }
  int32_t group() const { return group_; }
  uint8_t n() const { return n_; }
  void set_index(int64_t id) { index_ = id; }
  void set_group(int group) { group_ = group; }
  void set_edge(const HalfEdge& e) { edge_ = e.index(); }
  HalfEdge& get_edge();
  const HalfEdge& get_edge() const;

 private:
  HalfMesh& mesh_;            // reference to the HalfMesh container
  half_t edge_{halfnull_t};   // index of one of the edges on this face in
                              // HalfMesh::edges_
  half_t index_{halfnull_t};  // index of this face in HalfMesh::faces_
  int32_t group_{-1};         // group index
  uint8_t n_{0};              // number of vertices in this face
};

/// @brief Represents a mesh using half-edges.
/// Entities (HalfNode, HalfEdge, HalfFace) reference each other using indices
/// into the arrays stored in this HalfMesh container. The type of these indices
/// is specified by the definition of half_t above.
class HalfMesh {
 public:
  HalfMesh(const Mesh& mesh) : vertices_(mesh.vertices().dim()) { build(mesh); }

  HalfNode& create_node(const double* x) {
    ASSERT(nodes_.size() == size_t(vertices_.n()));
    if (available_node_.size() > 0) {
      size_t id = available_node_.top();
      ASSERT(!nodes_[id].active());
      available_node_.pop();
      nodes_[id].set_index(id);
      for (int d = 0; d < 3; d++) vertices_[id][d] = x[d];
      return nodes_[id];
    }
    size_t id = nodes_.size();
    nodes_.emplace_back(*this, x, id);
    return nodes_[id];
  }

  HalfEdge& create_edge() {
    if (available_edge_.size() > 0) {
      size_t id = available_edge_.top();
      ASSERT(!edges_[id].active());
      available_edge_.pop();
      edges_[id].set_index(id);
      return edges_[id];
    }
    size_t id = edges_.size();
    edges_.emplace_back(*this, id);
    return edges_[id];
  }

  HalfFace& create_face(int n) {
    if (available_face_.size() > 0) {
      size_t id = available_face_.top();
      ASSERT(!faces_[id].active());
      available_face_.pop();
      faces_[id].set_index(id);
      return faces_[id];
    }
    size_t id = faces_.size();
    faces_.emplace_back(*this, id, n);
    return faces_[id];
  }

  void extract(Mesh& mesh) const;
  bool check() const;

  void flip(half_t iedge) { flip(edges_[iedge]); }
  void flip(HalfEdge& edge);
  void split(half_t iedge, const coord_t* x);
  void insert(half_t iface, const double* x);
  void collapse(half_t iedge);

  size_t n_faces() const { return faces_.size(); }
  int get_connected_components(std::vector<int>& components) const;

  const auto& nodes() const { return nodes_; }
  const auto& edges() const { return edges_; }
  const auto& faces() const { return faces_; }

  auto& nodes() { return nodes_; }
  auto& edges() { return edges_; }
  auto& faces() { return faces_; }
  auto& vertices() { return vertices_; }

  void deactivate(HalfNode& node) {
    // available_node_.push(node.index());
    node.deactivate();
  }

  void deactivate(HalfEdge& edge) {
    // available_edge_.push(edge.index());
    edge.deactivate();
  }

  void deactivate(HalfFace& face) {
    // available_face_.push(face.index());
    face.deactivate();
  }

  template <typename fn>
  void deactivate_by(const fn& mask) {
    std::vector<HalfFace*> faces;
    for (auto& n : nodes_) {
      if (!n.active()) continue;
      if (mask(n)) {
        n.deactivate();
        n.get_onering(faces);
        for (auto* f : faces) f->deactivate();
      }
    }
  }

  template <typename fn>
  void activate_faces_by(const fn& predicate) {
    for (auto& f : faces_) f.deactivate();
    for (auto& n : nodes_) n.deactivate();
    for (size_t k = 0; k < faces_.size(); k++) {
      auto& f = faces_[k];
      if (predicate(k)) {
        f.set_index(k);
        // activate nodes of this face
        half_t e = f.edge();
        half_t first = e;
        do {
          half_t id = edges_[e].node();
          nodes_[id].set_index(id);
          e = edges_[e].next();
        } while (first != e);
      }
    }
  }

 private:
  void build(const Mesh& mesh);
  Vertices vertices_;
  std::vector<HalfNode> nodes_;
  std::vector<HalfEdge> edges_;
  std::vector<HalfFace> faces_;

  std::stack<size_t> available_node_;
  std::stack<size_t> available_edge_;
  std::stack<size_t> available_face_;
};

}  // namespace vortex