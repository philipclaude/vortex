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
#include "halfedges.h"

#include "numerics.h"
#include "stlext.h"

namespace vortex {

void HalfMesh::build(const Mesh& mesh) {
  // allocate memory for the mesh
  const Topology<Triangle>& triangles = mesh.triangles();
  const Topology<Quad>& quads = mesh.quads();
  const Topology<Polygon>& polygons = mesh.polygons();
  edges_.reserve(6 * triangles.n() + 8 * quads.n() + 12 * polygons.n());
  nodes_.reserve(mesh.vertices().n());
  faces_.reserve(triangles.n() + quads.n() + polygons.n());

  // create nodes
  vertices_.reserve(mesh.vertices().n());
  for (index_t k = 0; k < mesh.vertices().n(); k++) {
    (void)create_node(mesh.vertices()[k]);
  }

  // create faces and edges
  std::unordered_map<std::pair<int, int>, half_t> twins;
  twins.reserve(edges_.capacity());

  auto create_half_edges = [&](const auto& topology) {
    for (size_t k = 0; k < topology.n(); k++) {
      const int n_vertices = topology.length(k);
      auto& face = create_face(n_vertices);
      face.set_group(topology.group(k));

      half_t prev{halfnull_t};
      half_t root{halfnull_t};
      for (int i = 0; i < n_vertices; i++) {
        // retrieve indices of the edge endpoint vertices
        int p = topology(k, i);
        int q = (i + 1 < n_vertices) ? topology(k, i + 1) : topology(k, 0);

        // create a new edge
        auto& edge = create_edge();
        edge.set_face(face);

        // connect the edge with the origin vertex
        HalfNode& node = nodes_[p];
        edge.set_node(node);
        node.set_edge(edge);

        // look for a twin
        auto it = twins.find({q, p});
        if (it != twins.end()) {
          // connect the edge with its twin
          half_t twin = it->second;
          edges_[twin].make_pair(edge);
        }
        twins.insert({{p, q}, edge.index()});

        // adjust the linked list around the face
        if (prev != halfnull_t) {
          edges_[prev].set_next(edge);
          edge.set_prev(edges_[prev]);
        } else {
          ASSERT(root == halfnull_t);
          root = edge.index();
        }
        prev = edge.index();
      }

      ASSERT(root != halfnull_t);
      face.set_edge(edges_[root]);          // first edge around the face
      edges_[root].set_prev(edges_[prev]);  // connect the first and last edges
      edges_[prev].set_next(edges_[root]);
    }
  };
  create_half_edges(triangles);
  create_half_edges(quads);
  create_half_edges(polygons);

  // go back through the edges and find those that do not have a twin
  int64_t n_boundary = 0;
  for (auto& iedge : edges_) {
    if (iedge.twin() != halfnull_t) continue;  // not a boundary
    n_boundary++;

    // create a new fictitious half edge
    // it is more useful to have this boundary half edge point to a null face
    HalfEdge& bedge = create_edge();

    bedge.set_node(edges_[iedge.next()].get_node());
    bedge.make_pair(iedge);
    ASSERT(bedge.face() == halfnull_t);
    ASSERT(bedge.next() == halfnull_t);
  }

  // go through the edges again, this time we will construct the next and prev
  // information for the boundary half edges
  for (auto& he : edges_) {
    if (he.face() != halfnull_t) continue;  // not a boundary edge
    ASSERT(he.next() == halfnull_t);

    // keep looking for another boundary edge
    half_t he_next = he.twin();
    while (edges_[he_next].face() != halfnull_t) {
      he_next = edges_[he_next].get_prev().twin();
    }

    he.set_next(edges_[he_next]);
    edges_[he_next].set_prev(he);
  }
}

bool HalfMesh::check() const {
  bool ok = true;

  for (const auto& node : nodes_) {
    if (!node.active()) continue;
    ASSERT(node.edge() != halfnull_t);

    int count = 1;
    half_t e = node.edge();
    while (true) {
      ASSERT(edges_[e].node() != halfnull_t);
      ASSERT(node.index() != halfnull_t);
      e = edges_[e].next();
      if (edges_[e].node() == node.index()) break;
      count++;
    }
    if (count != edges_[e].get_face().n()) ok = false;
  }

  for (const auto& face : faces_) {
    if (!face.active()) continue;
    if (face.edge() == halfnull_t) ok = false;
  }

  size_t n_issues = 0;
  size_t n_zero_length = 0;
  for (auto& edge : edges_) {
    if (!edge.active()) continue;
    if (!edge.boundary() && edge.get_triangle_left_node().index() ==
                                edge.get_triangle_right_node().index()) {
      n_issues++;
      LOG << fmt::format("non-manifold edge {} {}, L/R = {}", edge.node(),
                         edge.get_twin().node(),
                         edge.get_triangle_left_node().index());
      ok = false;
    }
    vec3d p(edge.get_node().point());
    vec3d q(edge.get_twin().get_node().point());
    double l = length(q - p);
    if (l == 0) {
      LOG << fmt::format("zero-length edge {} {}", edge.node(),
                         edge.get_twin().node());
      n_zero_length++;
      ok = false;
    }
  }
  if (n_issues != 0) LOG << fmt::format("# non-manifold edges = {}", n_issues);
  if (n_zero_length != 0)
    LOG << fmt::format("detected {} zero-length edges", n_zero_length);
  return ok;
}

void HalfMesh::extract(Mesh& mesh) const {
  mesh.vertices().clear();
  mesh.triangles().clear();
  mesh.quads().clear();
  mesh.polygons().clear();

  // create the vertices
  mesh.vertices().set_dim(vertices_.dim());
  std::unordered_map<half_t, index_t> node_map;
  node_map.reserve(nodes_.size());
  for (const auto& node : nodes_) {
    if (node.index() == halfnull_t) continue;
    mesh.vertices().add(vertices_[node.index()]);
    node_map.insert({node.index(), node_map.size()});
  }

  // create the faces
  std::vector<index_t> elem(16);
  for (const auto& face : faces_) {
    if (face.index() == halfnull_t) continue;
    half_t e = face.edge();
    elem.clear();
    elem.resize(face.n());
    for (int j = 0; j < face.n(); j++) {
      ASSERT(node_map.find(edges_[e].node()) != node_map.end())
          << edges_[e].node();
      elem[j] = node_map.at(edges_[e].node());
      e = edges_[e].next();
    }
    ASSERT(face.n() > 2);
    if (elem.size() == 3)
      mesh.triangles().add(elem.data());
    else if (elem.size() == 4)
      mesh.quads().add(elem.data());
    else
      mesh.polygons().add(elem.data(), elem.size());
  }

  // map lines
  for (size_t k = 0; k < mesh.lines().n(); k++) {
    auto* e = mesh.lines()[k];
    ASSERT(node_map.find(e[0]) != node_map.end()) << e[0];
    ASSERT(node_map.find(e[1]) != node_map.end()) << e[1];
    e[0] = node_map.at(e[0]);
    e[1] = node_map.at(e[1]);
  }
}

void HalfMesh::flip(HalfEdge& edge) {
  auto& twin = edge.get_twin();

  auto& p = edge.get_node();
  auto& q = twin.get_node();

  ASSERT(p.index() != q.index())
      << fmt::format("error flipping edge {} with vertices {}-{}", edge.index(),
                     p.index(), q.index());

  auto& vL = edge.get_triangle_left_node();
  auto& vR = edge.get_triangle_right_node();

  // get the faces adjacent to this edge
  ASSERT(edge.face() != halfnull_t);  // for now assume a closed mesh
  ASSERT(twin.face() != halfnull_t);
  auto& fL = edge.get_face();
  auto& fR = twin.get_face();

  // get the edges
  auto& eLR = twin.get_next();
  auto& eUR = eLR.get_next();
  auto& eUL = edge.get_next();
  auto& eLL = eUL.get_next();

  // top
  edge.set_next(eUR);
  edge.set_node(vL);
  edge.set_face(fL);
  edge.set_twin(twin);

  eUR.set_next(eUL);
  eUR.set_node(vR);
  eUR.set_face(fL);
  // eUR.twin (same)

  eUL.set_next(edge);
  eUL.set_node(q);
  eUL.set_face(fL);
  // eUL.twin (same)

  // bottom
  twin.set_next(eLL);
  twin.set_node(vR);
  twin.set_face(fR);
  twin.set_twin(edge);

  eLL.set_next(eLR);
  eLL.set_node(vL);
  eLL.set_face(fR);
  // eLL.twin (same)

  eLR.set_next(twin);
  eLR.set_node(p);
  eLR.set_face(fR);
  // eLR.twin (same)

  // vertices
  vL.set_edge(edge);
  vR.set_edge(eUR);
  p.set_edge(eLR);
  q.set_edge(eUL);

  // faces
  fL.set_edge(edge);
  fR.set_edge(twin);
}

void HalfMesh::split(half_t iedge, const coord_t* x) {
  // create a new half vertex at the requested point
  HalfNode& v = create_node(x);

  // create new half edges
  HalfEdge& h = create_edge();
  HalfEdge& ht = create_edge();
  HalfEdge& hUL = create_edge();
  HalfEdge& hLL = create_edge();
  HalfEdge& hUR = create_edge();
  HalfEdge& hLR = create_edge();

  // create new half faces
  HalfFace& fUL = create_face(3);
  HalfFace& fUR = create_face(3);

  HalfEdge& edge = edges_[iedge];
  HalfEdge& twin = edge.get_twin();

  HalfNode& p = edge.get_node();
  HalfNode& q = twin.get_node();

  HalfNode& vL = edge.get_triangle_left_node();
  HalfNode& vR = edge.get_triangle_right_node();

  // get the faces adjacent to this edge, assuming a closed mesh
  ASSERT(edge.face() != halfnull_t);
  ASSERT(edge.face() != halfnull_t);
  HalfFace& fL = edge.get_face();
  HalfFace& fR = twin.get_face();

  // get the lower edges
  HalfEdge& eLR = twin.get_next();
  HalfEdge& eUR = eLR.get_next();
  HalfEdge& eUL = edge.get_next();
  HalfEdge& eLL = eUL.get_next();

  // left side
  edge.set_next(hLL);
  edge.set_face(fL);
  edge.set_node(p);
  edge.set_twin(twin);

  hLL.set_next(eLL);
  hLL.set_face(fL);
  hLL.set_node(v);
  hLL.set_twin(hUL);

  eLL.set_next(edge);
  eLL.set_node(vL);
  eLL.set_face(fL);
  // eLL.twin (same)

  hUL.set_next(h);
  hUL.set_face(fUL);
  hUL.set_node(vL);
  hUL.set_twin(hLL);

  h.set_next(eUL);
  h.set_face(fUL);
  h.set_node(v);
  h.set_twin(ht);

  eUL.set_next(hUL);
  eUL.set_node(q);
  eUL.set_face(fUL);
  // eUL.twin (same)

  // right side
  twin.set_next(eLR);
  twin.set_face(fR);
  twin.set_node(v);
  twin.set_twin(edge);

  eLR.set_next(hLR);
  eLR.set_face(fR);
  eLR.set_node(p);
  // eLR.twin (same)

  hLR.set_next(twin);
  hLR.set_face(fR);
  hLR.set_node(vR);
  hLR.set_twin(hUR);

  hUR.set_next(eUR);
  hUR.set_face(fUR);
  hUR.set_node(v);
  hUR.set_twin(hLR);

  eUR.set_next(ht);
  eUR.set_face(fUR);
  eUR.set_node(vR);
  // eUR.twin (same)

  ht.set_next(hUR);
  ht.set_face(fUR);
  ht.set_node(q);
  ht.set_twin(h);

  // vertices
  p.set_edge(edge);
  vL.set_edge(eLL);
  vR.set_edge(eUR);
  q.set_edge(ht);
  v.set_edge(h);

  // faces
  fL.set_edge(edge);
  fR.set_edge(twin);
  fUL.set_edge(h);
  fUR.set_edge(ht);
  fUL.set_group(fL.group());
  fUR.set_group(fR.group());
}

void HalfMesh::insert(half_t iface, const double* x) {
  HalfEdge& e3 = create_edge();
  HalfEdge& e4 = create_edge();
  HalfEdge& e5 = create_edge();
  HalfEdge& e6 = create_edge();
  HalfEdge& e7 = create_edge();
  HalfEdge& e8 = create_edge();

  HalfFace& t0 = create_face(3);
  HalfFace& t1 = create_face(3);
  HalfNode& n = create_node(x);

  HalfFace& face = faces_[iface];
  auto& e0 = face.get_edge();
  auto& e1 = e0.get_next();
  auto& e2 = e1.get_next();

  auto& n0 = e0.get_node();
  auto& n1 = e1.get_node();
  auto& n2 = e2.get_node();

  n.set_edge(e7);

  // face
  face.set_edge(e0);
  n0.set_edge(e0);
  e0.set_node(n0);
  e0.set_next(e5);
  e0.set_face(face);
  e5.set_node(n1);
  e5.set_twin(e6);
  e5.set_next(e7);
  e5.set_face(face);
  e7.set_node(n);
  e7.set_twin(e8);
  e7.set_next(e0);
  e7.set_face(face);

  // t0
  t0.set_group(face.group());
  t0.set_edge(e1);
  n1.set_edge(e1);
  e1.set_node(n1);
  e1.set_next(e3);
  e1.set_face(t0);
  e3.set_node(n2);
  e3.set_twin(e4);
  e3.set_next(e6);
  e3.set_face(t0);
  e6.set_node(n);
  e6.set_twin(e5);
  e6.set_next(e1);
  e6.set_face(t0);

  // t1
  t1.set_group(face.group());
  t1.set_edge(e2);
  n2.set_edge(e2);
  e2.set_node(n2);
  e2.set_next(e8);
  e2.set_face(t1);
  e8.set_node(n0);
  e8.set_twin(e7);
  e8.set_face(t1);
  e8.set_next(e4);
  e4.set_node(n);
  e4.set_twin(e3);
  e4.set_face(t1);
  e4.set_next(e2);
}

void HalfMesh::collapse(half_t iedge) {
  HalfEdge& edge = edges_[iedge];
  HalfEdge& twin = edge.get_twin();

  HalfNode& p = edge.get_node();  // vertex that will be removed
  HalfNode& q = twin.get_node();  // receiving vertex

  HalfNode& vL = edge.get_triangle_left_node();
  HalfNode& vR = edge.get_triangle_right_node();

  // get the faces adjacent to this edge
  HalfFace& fL = edge.get_face();
  HalfFace& fR = twin.get_face();

  // we are assuming a closed mesh
  ASSERT(fL.index() != halfnull_t);
  ASSERT(fR.index() != halfnull_t);

  // loop through the one-ring of the vertex
  // and re-assign the vertex of the affected half edges
  HalfEdge* e = &edge;
  int n_edges = 0;
  do {
    // assign the vertex of this half edge to the receiving vertex
    e->set_node(q);

    // go to the next edge
    e = &e->get_twin().get_next();

    ASSERT(size_t(n_edges++) < edges_.size());
  } while (e->index() != edge.index());

  // get the lower edges
  HalfEdge& edgeLR = twin.get_next();
  HalfEdge& edgeLL = edge.get_next().get_next();

  // get the upper edges
  HalfEdge& edgeUR = twin.get_next().get_next();
  HalfEdge& edgeUL = edge.get_next();

  // glue the lower and upper edges together on the right side
  edgeUR.get_twin().set_twin(edgeLR.get_twin());
  edgeLR.get_twin().set_twin(edgeUR.get_twin());

  // glue the lower and upper edges together on the left side
  edgeUL.get_twin().set_twin(edgeLL.get_twin());
  edgeLL.get_twin().set_twin(edgeUL.get_twin());

  // re-assign the vertex edges in case they pointed to deleted ones
  vL.set_edge(edgeUL.get_twin());
  vR.set_edge(edgeLR.get_twin());
  q.set_edge(edgeUR.get_twin());

  ASSERT(edgeUR.get_twin().node() == q.index());
  ASSERT(edgeLL.get_twin().node() == q.index());

  // deactivate deleted edges
  deactivate(edge);
  deactivate(twin);
  deactivate(edgeLR);
  deactivate(edgeLL);
  deactivate(edgeUR);
  deactivate(edgeUL);

  // deactivate deleted faces
  deactivate(fL);
  deactivate(fR);

  // deactivate the deleted vertex
  deactivate(p);
}

int HalfMesh::get_connected_components(std::vector<int>& components) const {
  int n_components = 0;
  components.reserve(faces_.size());
  std::unordered_map<half_t, int> face_component;
  face_component.reserve(faces_.size());
  while (true) {
    const HalfFace* face = nullptr;
    for (const auto& f : faces_) {
      if (face_component.find(f.index()) != face_component.end()) continue;
      face = &f;
      break;
    }
    if (face == nullptr) break;

    LOG << "processing component " << n_components;
    std::vector<const HalfFace*> stack;
    stack.reserve(64);
    stack.push_back(face);
    while (!stack.empty()) {
      auto* f = stack.back();
      stack.pop_back();
      face_component.insert({f->index(), n_components});
      half_t e = f->edge();
      do {
        half_t etf = edges_[e].get_twin().face();
        if (face_component.find(etf) == face_component.end())
          stack.push_back(&faces_[etf]);
        e = edges_[e].next();
      } while (e != f->edge());
    }
    n_components++;
  }

  // save the component indices
  for (auto& [f, c] : face_component) components[f] = c;

  return n_components;
}

HalfNode::HalfNode(HalfMesh& mesh, const double* x, int64_t index)
    : mesh_(mesh), index_(index) {
  mesh_.vertices().add(x);
}
void HalfNode::set_edge(const HalfEdge& e) { edge_ = e.index(); }
const double* HalfNode::point() const { return mesh_.vertices()[index_]; }
double* HalfNode::point() { return mesh_.vertices()[index_]; }
HalfEdge& HalfNode::get_edge() { return mesh_.edges()[edge_]; }
const HalfEdge& HalfNode::get_edge() const { return mesh_.edges()[edge_]; }

HalfEdge& HalfFace::get_edge() { return mesh_.edges()[edge_]; }
const HalfEdge& HalfFace::get_edge() const { return mesh_.edges()[edge_]; }
HalfEdge& HalfEdge::get_next() { return mesh_.edges()[next_]; }
const HalfEdge& HalfEdge::get_next() const { return mesh_.edges()[next_]; }
HalfEdge& HalfEdge::get_prev() { return mesh_.edges()[prev_]; }
const HalfEdge& HalfEdge::get_prev() const { return mesh_.edges()[prev_]; }
HalfNode& HalfEdge::get_node() { return mesh_.nodes()[node_]; }
const HalfNode& HalfEdge::get_node() const { return mesh_.nodes()[node_]; }

HalfEdge& HalfEdge::get_twin() { return mesh_.edges()[twin_]; }
const HalfEdge& HalfEdge::get_twin() const { return mesh_.edges()[twin_]; }

bool HalfEdge::boundary() const {
  return (face_ == halfnull_t) || (get_twin().face() == halfnull_t);
}

HalfNode& HalfEdge::get_triangle_left_node() {
  return mesh_.nodes()[get_next().get_next().node()];
}
HalfNode& HalfEdge::get_triangle_right_node() {
  return mesh_.nodes()[get_twin().get_next().get_next().node()];
}

const HalfNode& HalfEdge::get_triangle_left_node() const {
  return mesh_.nodes()[get_next().get_next().node()];
}
const HalfNode& HalfEdge::get_triangle_right_node() const {
  return mesh_.nodes()[get_twin().get_next().get_next().node()];
}

HalfFace& HalfEdge::get_face() { return mesh_.faces()[face_]; }
const HalfFace& HalfEdge::get_face() const { return mesh_.faces()[face_]; }
void HalfEdge::set_face(const HalfFace& f) { face_ = f.index(); }

template <typename T>
void HalfNode::get_onering(std::vector<T*>& ring) const {
  ring.clear();
  half_t first = edge_;
  half_t edge = first;
  do {
    if constexpr (std::is_same<T, HalfEdge>::value)
      ring.push_back(&mesh_.edges()[edge]);
    else if constexpr (std::is_same<T, HalfFace>::value)
      ring.push_back(&mesh_.faces()[mesh_.edges()[edge].face()]);
    else if constexpr (std::is_same<T, HalfNode>::value) {
      const auto& t = mesh_.edges()[edge].get_twin();
      ring.push_back(&mesh_.nodes()[t.node()]);
    }
    edge = mesh_.edges()[edge].get_twin().next();
  } while (edge != first);
}

template void HalfNode::get_onering(std::vector<HalfNode*>&) const;
template void HalfNode::get_onering(std::vector<HalfEdge*>&) const;
template void HalfNode::get_onering(std::vector<HalfFace*>&) const;

}  // namespace vortex