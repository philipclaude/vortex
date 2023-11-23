#include "halfedges.h"

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
    auto& node = create_node(mesh.vertices()[k]);
    ASSERT(node.index == k);
    vertices_.set_entity(k, mesh.vertices().entity(k));
    node.edge = nullptr;
  }

  // create faces and edges
  std::unordered_map<std::pair<int, int>, HalfEdge*> twins;
  twins.reserve(edges_.capacity());

  auto create_half_edges = [&](const auto& topology) {
    for (int k = 0; k < topology.n(); k++) {
      const int n_vertices = topology.length(k);
      auto& face = create_face(n_vertices);
      face.group = topology.group(k);

      HalfEdge* prev{nullptr};
      HalfEdge* root{nullptr};
      for (int i = 0; i < n_vertices; i++) {
        // retrieve indices of the edge endpoint vertices
        int p = topology(k, i);
        int q = (i + 1 < n_vertices) ? topology(k, i + 1) : topology(k, 0);

        // create a new edge
        HalfEdge& edge = create_edge();
        edge.face = &face;

        // connect the edge with the origin vertex
        HalfNode& node = nodes_[p];
        edge.node = &node;
        node.edge = &edge;  // written many times

        // look for a twin
        auto it = twins.find({q, p});
        if (it != twins.end()) {
          // connect the edge with its twin
          HalfEdge* twin = it->second;
          twin->twin = &edge;
          edge.twin = twin;
        }
        twins.insert({{p, q}, &edge});

        // adjust the linked list around the face
        if (prev != nullptr) {
          prev->next = &edge;
          edge.prev = prev;
        } else {
          ASSERT(root == nullptr);
          root = &edge;
        }
        prev = &edge;
      }

      ASSERT(root != nullptr);
      face.edge = root;   // first edge in the linked list around the face
      root->prev = prev;  // connect the first and last edges
      prev->next = root;
    }
  };
  create_half_edges(triangles);
  create_half_edges(quads);
  create_half_edges(polygons);

  // go back through the edges and find those that do not have a twin
  int64_t n_boundary = 0;
  for (auto& iedge : edges_) {
    if (iedge.twin != nullptr) continue;  // not a boundary
    n_boundary++;

    // create a new fictitious half edge
    // it is more useful to have this boundary half edge point to a null face
    HalfEdge& bedge = create_edge();

    bedge.node = iedge.next->node;
    bedge.twin = &iedge;
    iedge.twin = &bedge;
    bedge.face = nullptr;

    // we cannot set the next yet as they may not be constructed
    bedge.next = nullptr;
  }

  // go through the edges again, this time we will construct the next and prev
  // information for the boundary half edges
  for (auto& he : edges_) {
    if (he.face != nullptr) continue;  // not a boundary edge
    ASSERT(he.next == nullptr);

    // keep looking for another boundary edge
    HalfEdge* he_next = he.twin;
    while (he_next->face != nullptr) {
      he_next = he_next->prev->twin;
    }

    he.next = he_next;
    he_next->prev = &he;
  }
}

bool HalfMesh::check() const {
  bool ok = true;

  for (const auto& node : nodes_) {
    if (node.edge == nullptr) ok = false;

    int count = 1;
    HalfEdge* e = node.edge;
    while (true) {
      e = e->next;
      if (e->node == &node) break;
      count++;
    }
    if (count != e->face->n) ok = false;
  }

  for (const auto& face : faces_) {
    if (face.edge == nullptr) ok = false;
  }
  return ok;
}

void HalfMesh::extract(Mesh& mesh) const {
  mesh.vertices().clear();
  mesh.triangles();
  mesh.quads().clear();
  mesh.polygons().clear();

  // create the vertices
  mesh.vertices().set_dim(vertices_.dim());
  std::unordered_map<const HalfNode*, index_t> node_map;
  node_map.reserve(nodes_.size());
  for (const auto& node : nodes_) {
    mesh.vertices().add(vertices_[node.index]);
    node_map.insert({&node, node_map.size()});
  }

  // create the faces
  std::vector<index_t> elem(16);
  for (const auto& face : faces_) {
    HalfEdge* e = face.edge;
    elem.clear();
    elem.resize(face.n);
    for (int j = 0; j < face.n; j++) {
      elem[j] = node_map.at(e->node);
      e = e->next;
    }
    ASSERT(face.n > 2);
    if (elem.size() == 3)
      mesh.triangles().add(elem.data());
    else if (elem.size() == 4)
      mesh.quads().add(elem.data());
    else
      mesh.polygons().add(elem.data(), elem.size());
  }
}

void HalfMesh::flip(HalfEdge* edge) {
  HalfNode* p = edge->node;
  HalfNode* q = edge->twin->node;

  ASSERT(p != q);

  HalfNode* vL = edge->next->next->node;
  HalfNode* vR = edge->twin->next->next->node;

  HalfEdge* twin = edge->twin;

  // get the faces adjacent to this edge
  HalfFace* fL = edge->face;
  HalfFace* fR = twin->face;

  // we are assuming a closed mesh
  ASSERT(fL != nullptr);
  ASSERT(fR != nullptr);

  // get the lower edges
  HalfEdge* eLR = twin->next;
  HalfEdge* eLL = edge->next->next;

  // get the upper edges
  HalfEdge* eUR = twin->next->next;
  HalfEdge* eUL = edge->next;

  // top
  edge->next = eUR;
  edge->node = vL;
  edge->face = fL;
  edge->twin = twin;

  eUR->next = eUL;
  eUR->node = vR;
  eUR->face = fL;
  // eUR->twin (same)

  eUL->next = edge;
  eUL->node = q;
  eUL->face = fL;
  // eUL->twin (same)

  // bottom
  twin->next = eLL;
  twin->node = vR;
  twin->face = fR;
  twin->twin = edge;

  eLL->next = eLR;
  eLL->node = vL;
  eLL->face = fR;
  // eLL->twin (same)

  eLR->next = twin;
  eLR->node = p;
  eLR->face = fR;
  // eLR->twin (same)

  // vertices
  vL->edge = edge;
  vR->edge = eUR;
  p->edge = eLR;
  q->edge = eUL;

  // faces
  fL->edge = edge;
  fR->edge = twin;
}

void HalfMesh::split(HalfEdge* edge, const coord_t* x) {
  HalfNode* p = edge->node;
  HalfNode* q = edge->twin->node;

  HalfNode* vL = edge->next->next->node;
  HalfNode* vR = edge->twin->next->next->node;

  HalfEdge* twin = edge->twin;

  // get the faces adjacent to this edge
  HalfFace* fL = edge->face;
  HalfFace* fR = twin->face;

  // we are assuming a closed mesh
  ASSERT(fL != nullptr);
  ASSERT(fR != nullptr);

  // get the lower edges
  HalfEdge* eLR = twin->next;
  HalfEdge* eLL = edge->next->next;

  // get the upper edges
  HalfEdge* eUR = twin->next->next;
  HalfEdge* eUL = edge->next;

  // create a new half vertex at the requested point
  HalfNode* v = &create_node(x);

  // create new half edges
  HalfEdge* h = &create_edge();
  HalfEdge* ht = &create_edge();
  HalfEdge* hUL = &create_edge();
  HalfEdge* hLL = &create_edge();
  HalfEdge* hUR = &create_edge();
  HalfEdge* hLR = &create_edge();

  // create new half faces
  HalfFace* fUL = &create_face(3);
  HalfFace* fUR = &create_face(3);

  // left side
  edge->next = hLL;
  edge->face = fL;
  edge->node = p;
  edge->twin = twin;

  hLL->next = eLL;
  hLL->face = fL;
  hLL->node = v;
  hLL->twin = hUL;

  eLL->next = edge;
  eLL->node = vL;
  eLL->face = fL;
  // eLL->twin (same)

  hUL->next = h;
  hUL->face = fUL;
  hUL->node = vL;
  hUL->twin = hLL;

  h->next = eUL;
  h->face = fUL;
  h->node = v;
  h->twin = ht;

  eUL->next = hUL;
  eUL->node = q;
  eUL->face = fUL;
  // eUL->twin (same)

  // right side
  twin->next = eLR;
  twin->face = fR;
  twin->node = v;
  twin->twin = edge;

  eLR->next = hLR;
  eLR->face = fR;
  eLR->node = p;
  // eLR->twin (same)

  hLR->next = twin;
  hLR->face = fR;
  hLR->node = vR;
  hLR->twin = hUR;

  hUR->next = eUR;
  hUR->face = fUR;
  hUR->node = v;
  hUR->twin = hLR;

  eUR->next = ht;
  eUR->face = fUR;
  eUR->node = vR;
  // eUR->twin (same)

  ht->next = hUR;
  ht->face = fUR;
  ht->node = q;
  ht->twin = h;

  // vertices
  p->edge = edge;
  vL->edge = eLL;
  vR->edge = eUR;
  q->edge = ht;
  v->edge = h;

  // faces
  fL->edge = edge;
  fR->edge = twin;
  fUL->edge = h;
  fUR->edge = ht;
}

void HalfMesh::collapse(HalfEdge* edge) {
  HalfNode* p = edge->node;        // vertex that will be removed
  HalfNode* q = edge->twin->node;  // receiving vertex

  // do not collapse boundary edges for now
  //   if (p->index != 0 || q->index != 0) {
  //     heap.erase(edge);
  //     return false;
  //   }

  HalfNode* vL = edge->next->next->node;
  HalfNode* vR = edge->twin->next->next->node;

  HalfEdge* twin = edge->twin;

  // get the faces adjacent to this edge
  HalfFace* fL = edge->face;
  HalfFace* fR = twin->face;

  // we are assuming a closed mesh
  ASSERT(fL != nullptr);
  ASSERT(fR != nullptr);

  // determine if the collapse is geometrically valid
  //   std::vector<HalfFace*> faces;
  //   get_onering(p, faces);
  //   bool flips = false;
  //   // check the first ring
  //   for (int i = 0; i < faces.size(); i++) {
  //     if (faces[i] == fL || faces[i] == fR) continue;
  //     if (face_flips(faces[i], p, q->point)) {
  //       flips = true;
  //       break;
  //     }
  //   }
  //   if (flips) {
  //     heap.erase(edge);
  //     return false;
  //   }

  // loop through the one-ring of the vertex
  // and re-assign the vertex of the affected half edges
  HalfEdge* e = edge;
  int nb_edges = 0;
  do {
    // assign the vertex of this half edge to the receiving vertex
    e->node = q;

    // go to the next edge
    e = e->twin->next;

    nb_edges++;
    if (size_t(nb_edges) > edges_.size()) NOT_POSSIBLE;

  } while (e != edge);

  // get the lower edges
  HalfEdge* edgeLR = twin->next;
  HalfEdge* edgeLL = edge->next->next;

  // get the upper edges
  HalfEdge* edgeUR = twin->next->next;
  HalfEdge* edgeUL = edge->next;

  // glue the lower and upper edges together on the right side
  edgeUR->twin->twin = edgeLR->twin;
  edgeLR->twin->twin = edgeUR->twin;

  // glue the lower and upper edges together on the left side
  edgeUL->twin->twin = edgeLL->twin;
  edgeLL->twin->twin = edgeUL->twin;

  // re-assign the vertex edges in case they previously pointed to deleted ones
  vL->edge = edgeUL->twin;
  vR->edge = edgeLR->twin;
  q->edge = edgeUR->twin;

  ASSERT(edgeUR->twin->node == q);
  ASSERT(edgeLL->twin->node == q);

  // remove the deleted edges
  HalfEdge* deleted_edges[6] = {edge, twin, edgeLR, edgeLL, edgeUR, edgeUL};
  for (int k = 0; k < 6; k++) {
    // mesh.remove(deleted_edges[k]);
  }

  // remove the deleted faces
  // mesh.remove(fL);
  // mesh.remove(fR);

  // remove the unused vertex
  // mesh.remove(p);

  // remove the edges from the heap
  //   heap.erase(edgeLR);
  //   heap.erase(edgeLL);
  //   heap.erase(edgeUR);
  //   heap.erase(edgeUL);
  //   heap.erase(twin);
  //   heap.erase(edge);
}

int HalfMesh::get_connected_components(std::vector<int>& components) const {
  int n_components = 0;
  components.reserve(faces_.size());
  std::unordered_map<const HalfFace*, int> face_component;
  face_component.reserve(faces_.size());
  while (true) {
    const HalfFace* face = nullptr;
    for (auto& f : faces_) {
      if (face_component.find(&f) != face_component.end()) continue;
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
      face_component.insert({f, n_components});
      auto* e = f->edge;
      do {
        if (e->twin->face &&
            face_component.find(e->twin->face) == face_component.end()) {
          stack.push_back(e->twin->face);
        }
        e = e->next;
      } while (e != f->edge);
    }
    n_components++;
  }

  // save the component indices
  for (auto& [f, c] : face_component) components[f->index] = c;

  return n_components;
}

template <typename T>
void HalfMesh::get_onering(const HalfNode* node, std::vector<T*>& ring) const {
  ring.clear();
  HalfEdge* first = node->edge;
  HalfEdge* edge = first;
  do {
    if constexpr (std::is_same<T, HalfEdge>::value)
      ring.push_back(edge);
    else if constexpr (std::is_same<T, HalfFace>::value)
      ring.push_back(edge->face);
    else if constexpr (std::is_same<T, HalfNode>::value)
      ring.push_back(edge->twin->node);
    edge = edge->twin->next;
  } while (edge != first);
}

template void HalfMesh::get_onering(const HalfNode*,
                                    std::vector<HalfNode*>&) const;
template void HalfMesh::get_onering(const HalfNode*,
                                    std::vector<HalfEdge*>&) const;
template void HalfMesh::get_onering(const HalfNode*,
                                    std::vector<HalfFace*>&) const;

}  // namespace vortex