#pragma once

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
  void set_edge(const HalfEdge& e);
  const double* point() const;

  template <typename T> void get_onering(std::vector<T*>& ring) const;

 private:
  int64_t index_{halfnull_t};
  half_t edge_{halfnull_t};
  HalfMesh& mesh_;
};

class HalfFace;
class HalfEdge {
 public:
  HalfEdge(HalfMesh& m, int64_t id) : mesh_(m), index_(id) {}
  HalfNode& get_node();
  const HalfNode& get_node() const;

  HalfEdge& get_twin();
  const HalfEdge& get_twin() const;

  HalfNode& get_triangle_left_node();
  HalfNode& get_triangle_right_node();
  const HalfNode& get_triangle_left_node() const;
  const HalfNode& get_triangle_right_node() const;
  HalfFace& get_face();
  const HalfFace& get_face() const;
  HalfEdge& get_next();
  const HalfEdge& get_next() const;

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
  HalfMesh& mesh_;
  half_t index_{halfnull_t};
  half_t twin_{halfnull_t};
  half_t node_{halfnull_t};
  half_t face_{halfnull_t};
  half_t prev_{halfnull_t};
  half_t next_{halfnull_t};
};

class HalfFace {
 public:
  HalfFace(HalfMesh& mesh, int64_t index, uint8_t n)
      : mesh_(mesh), index_(index), n_(n) {}

  half_t index() const { return index_; }
  half_t edge() const { return edge_; }
  int32_t group() const { return group_; }
  uint8_t n() const { return n_; }
  void set_group(int group) { group_ = group; }
  void set_edge(const HalfEdge& e) { edge_ = e.index(); }
  HalfEdge& get_edge();
  const HalfEdge& get_edge() const;

 private:
  half_t edge_{halfnull_t};
  half_t index_{halfnull_t};
  int32_t group_{-1};
  uint8_t n_{0};
  HalfMesh& mesh_;
};

class HalfMesh {
 public:
  HalfMesh(const Mesh& mesh) : vertices_(mesh.vertices().dim()) { build(mesh); }

  HalfNode& create_node(const double* x) {
    ASSERT(nodes_.size() == size_t(vertices_.n()));
    size_t id = nodes_.size();
    nodes_.emplace_back(*this, x, id);
    return nodes_[id];
  }

  HalfEdge& create_edge() {
    size_t id = edges_.size();
    edges_.emplace_back(*this, id);
    return edges_[id];
  }

  HalfFace& create_face(int n) {
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

 private:
  void build(const Mesh& mesh);
  Vertices vertices_;
  std::vector<HalfNode> nodes_;
  std::vector<HalfEdge> edges_;
  std::vector<HalfFace> faces_;
};

}  // namespace vortex