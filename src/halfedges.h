#pragma once

#include <vector>

#include "log.h"
#include "mesh.h"

namespace vortex {

class HalfMesh {
 public:
  struct HalfEdge;
  struct HalfNode {
    int64_t index{-1};
    HalfEdge* edge{nullptr};
  };
  struct HalfFace;
  struct HalfEdge {
    int64_t index{-1};
    HalfFace* face{nullptr};
    HalfNode* node{nullptr};
    HalfEdge* prev{nullptr};
    HalfEdge* next{nullptr};
    HalfEdge* twin{nullptr};
  };
  struct HalfFace {
    HalfEdge* edge{nullptr};
    int64_t index{-1};
    int64_t group{0};
    uint8_t n{0};
  };

  HalfMesh(const Mesh& mesh) : vertices_(mesh.vertices().dim()) { build(mesh); }

  HalfNode& create_node(const double* x) {
    ASSERT(nodes_.size() == size_t(vertices_.n()));
    size_t id = nodes_.size();
    vertices_.add(x);
    nodes_.emplace_back();
    nodes_[id].index = id;
    return nodes_[id];
  }

  HalfEdge& create_edge() {
    size_t id = edges_.size();
    edges_.emplace_back();
    return edges_[id];
  }

  HalfFace& create_face(int n) {
    size_t id = faces_.size();
    faces_.emplace_back();
    faces_[id].index = id;
    faces_[id].n = n;
    return faces_[id];
  }

  void extract(Mesh& mesh) const;
  bool check() const;
  void flip(HalfEdge* edge);
  void split(HalfEdge* edge, const coord_t* x);
  void insert(HalfFace* face, const coord_t* x);
  void collapse(HalfEdge* edge);

  template <typename T>
  void get_onering(const HalfNode* node, std::vector<T*>& ring) const;

  size_t n_faces() const { return faces_.size(); }
  int get_connected_components(std::vector<int>& components) const;

  const auto& nodes() const { return nodes_; }
  const auto& edges() const { return edges_; }
  const auto& faces() const { return faces_; }

  auto& nodes() { return nodes_; }
  auto& edges() { return edges_; }
  auto& faces() { return faces_; }

 private:
  void build(const Mesh& mesh);
  Vertices vertices_;
  std::vector<HalfNode> nodes_;
  std::vector<HalfEdge> edges_;
  std::vector<HalfFace> faces_;
};

}  // namespace vortex
