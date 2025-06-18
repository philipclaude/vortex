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
#include "triangulate.h"

#include <queue>
#include <unordered_set>

#include "io.h"
#include "library.h"
#include "math/linalg.h"
#include "mesh.h"
#include "numerics.h"
#include "predicates.h"
#include "trees/kdtree.h"

namespace vortex {

using Face_t = HalfFace;

OceanTriangulator::OceanTriangulator(Mesh& mesh, const Mesh& coast)
    : mesh_(mesh), coast_(coast), hmesh_(mesh_) {}

void OceanTriangulator::triangulate() {
  insert_points();
  recover_edges();
  hmesh_.extract(mesh_);
}

vec3d get_triangle_barycentric(const Face_t& face, const double* x) {
  // extract 3d coordinates of triangle and query point
  const auto& e0 = face.get_edge();
  const auto& e1 = e0.get_next();
  const auto& e2 = e1.get_next();

  vec3d pa(e0.get_node().point(), 3);
  vec3d pb(e1.get_node().point(), 3);
  vec3d pc(e2.get_node().point(), 3);
  vec3d pd(x, 3);

  vec3d alpha = {-1, -1, -1};
  if (length(pd - pa) > 1) return alpha;
  if (length(pd - pb) > 1) return alpha;
  if (length(pd - pc) > 1) return alpha;

  // cast a ray from the sphere center to the point
  mat3 a;
  vec3d b;
  for (int i = 0; i < 3; i++) {
    a(i, 0) = pa[i] - pc[i];
    a(i, 1) = pb[i] - pc[i];
    a(i, 2) = -pd[i];
    b[i] = -pc[i];
  }
  vec3d solution = vortex::inverse(a) * b;

  alpha[0] = solution[0];
  alpha[1] = solution[1];
  alpha[2] = 1 - alpha[0] - alpha[1];
  return alpha;
}

double face_area(const HalfFace& face) {
  const auto& e0 = face.get_edge();
  const auto& e1 = e0.get_next();
  const auto& e2 = e1.get_next();
  const auto& n0 = e0.get_node();
  const auto& n1 = e1.get_node();
  const auto& n2 = e2.get_node();
  return face_area(n0.point(), n1.point(), n2.point());
}

Face_t* search_for(HalfMesh& mesh, half_t root, const double* x) {
  std::queue<half_t> queue;
  std::unordered_set<index_t> visited;

  queue.push(root);
  visited.insert(root);
  while (!queue.empty()) {
    half_t iface = queue.front();
    Face_t& face = mesh.faces()[iface];
    queue.pop();
    double tol = 1e-6;
    vec3d alpha = get_triangle_barycentric(face, x);
    bool inface = (alpha[0] >= -tol && alpha[0] <= 1 + tol) &&
                  (alpha[1] >= -tol && alpha[1] <= 1 + tol) &&
                  (alpha[2] >= -tol && alpha[2] <= 1 + tol);
    if (inface) return &face;

    // check which neighbors should be searched
    const auto& e0 = face.get_edge();
    const auto& e1 = e0.get_next();
    const auto& e2 = e1.get_next();
    const auto& face0 = e0.get_twin().get_face();
    const auto& face1 = e1.get_twin().get_face();
    const auto& face2 = e2.get_twin().get_face();

    half_t faces[3] = {face0.index(), face1.index(), face2.index()};
    for (int i = 0; i < 3; i++) {
      if (faces[i] == halfnull_t) continue;
      if (visited.find(faces[i]) == visited.end()) {
        visited.insert(faces[i]);
        queue.push(faces[i]);
      }
    }
  }

  return nullptr;
}

double quality(const vec3d& pa, const vec3d& pb, const vec3d& pc) {
  double area = 0.5 * length(cross(pb - pa, pc - pa));
  double lab = length(pb - pa);
  double lbc = length(pc - pb);
  double lac = length(pc - pa);
  return area / (lab * lab + lbc * lbc + lac * lac);
}

double quality(const Face_t& face) {
  const auto& e0 = face.get_edge();
  const auto& e1 = e0.get_next();
  const auto& e2 = e1.get_next();
  vec3d pa(e0.get_node().point(), 3);
  vec3d pb(e1.get_node().point(), 3);
  vec3d pc(e2.get_node().point(), 3);
  return quality(pa, pb, pc);
}

bool valid_flip(const HalfEdge& edge, bool check_quality = true) {
  // extract 3d coordinates of edge endpoints
  const HalfEdge& twin = edge.get_twin();

  vec3d pa(edge.get_node().point(), 3);
  vec3d pb(twin.get_node().point(), 3);

  // extract 3d coordinates of left & right vertices
  vec3d pl(edge.get_triangle_left_node().point(), 3);
  vec3d pr(edge.get_triangle_right_node().point(), 3);

  auto normal = [](const vec3d& a, const vec3d& b, const vec3d& c) {
    return normalize(cross(b - a, c - a));
  };

  if (check_quality) {
    double q0 = std::min(quality(pa, pb, pl), quality(pa, pr, pb));
    double q1 = std::min(quality(pl, pr, pb), quality(pa, pr, pl));
    if (q1 < q0) return false;  // quality should improve
  }

  double a0 = face_area(&pl[0], &pr[0], &pb[0]);
  double a1 = face_area(&pa[0], &pr[0], &pl[0]);
  if (a0 <= 0 || a1 <= 0) return false;

  // double d0 = dot(normal(pa, pb, pl), normal(pa, pr, pb));
  double d1 = dot(normal(pl, pr, pb), normal(pa, pr, pl));
  if (d1 < 1e-2) return false;
  // if (d1 < d0) return false;  // dot product should increase

  return true;
}

size_t optimize(HalfMesh& mesh, HalfNode& n) {
  // list all nodes
  std::vector<HalfNode*> nodes0;
  n.get_onering(nodes0);
  // nodes.push_back(&n);
  std::set<HalfNode*> nodes1;
  for (auto node : nodes0) {
    std::vector<HalfNode*> nodes2;
    node->get_onering(nodes2);
    for (auto* node2 : nodes2) nodes1.insert(node2);
  }
  std::vector<HalfNode*> nodes(nodes1.begin(), nodes1.end());

  auto cmp = [&mesh](half_t a, half_t b) {
    return quality(mesh.faces()[a]) < quality(mesh.faces()[b]);
  };
  std::set<half_t, decltype(cmp)> faces(cmp);

  size_t n_flips = 0;
  bool flipped = true;
  do {
    // build list of faces in order of quality (lowest to highest)
    faces.clear();
    for (const auto* node : nodes) {
      half_t e = node->edge();
      half_t first = e;
      do {
        faces.insert(mesh.edges()[e].face());
        e = mesh.edges()[e].next();
      } while (e != first);
    }

    // try flipping faces
    flipped = false;
    for (auto f : faces) {
      // retrieve edge and neighbor information
      auto& face = mesh.faces()[f];
      const auto& e0 = face.get_edge();
      const auto& e1 = e0.get_next();
      const auto& e2 = e1.get_next();
      const auto& f0 = e0.get_twin().get_face();
      const auto& f1 = e1.get_twin().get_face();
      const auto& f2 = e2.get_twin().get_face();

      // determine which edge would be best to flip
      double q0 = quality(f0);
      double q1 = quality(f1);
      double q2 = quality(f2);
      flipped = true;
      if (q0 < q1 && q0 < q2 && valid_flip(e0)) {
        mesh.flip(e0.index());
      } else if (q1 < q0 && q1 < q2 && valid_flip(e1)) {
        mesh.flip(e1.index());
      } else if (q2 < q0 && q2 < q1 && valid_flip(e2)) {
        mesh.flip(e2.index());
      } else
        flipped = false;
      if (flipped) {
        n_flips++;
        break;  // need to restart
      }
    }
  } while (flipped);
  return n_flips;
}

size_t optimize_mesh(HalfMesh& mesh) {
  // create initial list of faces
  auto cmp = [&mesh](half_t a, half_t b) {
    return quality(mesh.faces()[a]) < quality(mesh.faces()[b]);
  };
  std::set<half_t, decltype(cmp)> faces(cmp);
  for (auto& f : mesh.faces()) faces.insert(f.index());

  size_t n_flips = 0;
  bool flipped = true;
  do {
    auto it = faces.begin();
    half_t f = *it;
    faces.erase(it);

    // retrieve edge and neighbor information
    auto& face = mesh.faces()[f];
    const auto& e0 = face.get_edge();
    const auto& e1 = e0.get_next();
    const auto& e2 = e1.get_next();
    const auto& f0 = e0.get_twin().get_face();
    const auto& f1 = e1.get_twin().get_face();
    const auto& f2 = e2.get_twin().get_face();

    // determine which edge would be best to flip
    double q0 = quality(f0);
    double q1 = quality(f1);
    double q2 = quality(f2);
    flipped = true;
    half_t t = halfnull_t;
    if (q0 < q1 && q0 < q2 && valid_flip(e0)) {
      mesh.flip(e0.index());
      t = f0.index();
    } else if (q1 < q0 && q1 < q2 && valid_flip(e1)) {
      mesh.flip(e1.index());
      t = f1.index();
    } else if (q2 < q0 && q2 < q1 && valid_flip(e2)) {
      mesh.flip(e2.index());
      t = f2.index();
    } else
      flipped = false;
    if (flipped) {
      faces.erase(t);
      faces.insert(f);
      faces.insert(t);
      n_flips++;
    }
  } while (faces.size() > 0);
  return n_flips;
}

void OceanTriangulator::insert_points() {
  // build a kdtree of the original surface points
  trees::KdTree<3, coord_t, index_t, true> tree(mesh_.vertices()[0],
                                                mesh_.vertices().n());

  // TODO: allocate faces and edges as well
  hmesh_.nodes().reserve(hmesh_.nodes().size() + coast_.vertices().n());

  std::vector<index_t> coast_map(coast_.vertices().n());

  // insert each point in the coast
  size_t n_flips = 0;
  size_t n_inserted = 0;
  size_t n_edge = 0;
  LOG << "inserting " << coast_.vertices().n();
  for (size_t i = 0; i < coast_.vertices().n(); i++) {
    // find the closest node (in the original mesh)
    const auto* qi = coast_.vertices()[i];
    index_t n = tree.nearest(qi);

    // find the face containing this point (in a triangulation of the sphere)
    auto& node = hmesh_.nodes()[n];
    half_t root = node.get_edge().face();
    auto* face = search_for(hmesh_, root, qi);
    ASSERT(face);

    const auto& e0 = face->get_edge();
    const auto& e1 = e0.get_next();
    const auto& e2 = e1.get_next();

    // determine if this is an insertion on an edge (skip for now)
    vec3d alpha = get_triangle_barycentric(*face, qi);
    constexpr double tol = 1e-6;
    size_t n_zero = std::count_if(&alpha[0], &alpha[2], [](const double& x) {
      return std::fabs(x) < tol;
    });
    if (n_zero == 2) {
      int64_t id = -1;
      // alpha should be close to 1: 0.5 is conservative
      if (alpha[0] > 0.5) id = e0.node();
      if (alpha[1] > 0.5) id = e1.node();
      if (alpha[2] > 0.5) id = e2.node();
      ASSERT(id >= 0);
      coast_map[i] = id;
      continue;
    }

    size_t id = hmesh_.nodes().size();
    coast_map[i] = id;

    bool edge = true;
    if (std::fabs(alpha[0]) <= tol) {
      hmesh_.split(e1.index(), qi);
    } else if (std::fabs(alpha[1]) <= tol) {
      hmesh_.split(e2.index(), qi);
    } else if (std::fabs(alpha[2]) <= tol) {
      hmesh_.split(e0.index(), qi);
    } else {
      // insert the point and flip edges to optimize quality
      edge = false;
      hmesh_.insert(face->index(), qi);
    }
    if (edge) n_edge++;
    n_flips += optimize(hmesh_, hmesh_.nodes()[id]);
    n_inserted++;
  }
  LOG << fmt::format("inserted {} (/ {}) ({} on edges), # flips = {}",
                     n_inserted, coast_.vertices().n(), n_edge, n_flips);

  n_flips = optimize_mesh(hmesh_);
  LOG << fmt::format("mesh optimization used {} flips", n_flips);

  // add the lines using the map for the inserted point indices
  for (size_t i = 0; i < coast_.lines().n(); i++) {
    index_t e[2];
    for (int j = 0; j < 2; j++) e[j] = coast_map[coast_.lines()(i, j)];
    mesh_.lines().add(e);
    mesh_.lines().set_group(i, coast_.lines().group(i));
  }
}

bool edge_intersects(HalfNode& n0, HalfNode& n1, HalfEdge* e) {
  auto& nl = e->get_node();
  auto& nr = e->get_twin().get_node();

  ASSERT(n0.index() != nl.index());

  vec3d p0(n0.point());
  vec3d p1(n1.point());
  vec3d pl(nl.point());
  vec3d pr(nr.point());

  vec3d u0, u1, ul, ur;
  get_params(p0, p1, pl, pr, u0, u1, ul, ur);

  double a0 = orient2d(&u0[0], &u1[0], &ul[0]);
  double a1 = orient2d(&u0[0], &ur[0], &u1[0]);
  double a2 = orient2d(&ul[0], &ur[0], &u1[0]);
  double a3 = orient2d(&u0[0], &ur[0], &ul[0]);

  if (a0 * a1 < 0) return false;
  if (a2 * a3 < 0) return false;

  return true;
}

bool recover_edge(HalfMesh& mesh, index_t e0, index_t e1) {
  auto& n0 = mesh.nodes()[e0];
  // auto& n1 = mesh.nodes()[e1];
  std::vector<HalfEdge*> edges;

  // extract the edges around e0
  n0.get_onering(edges);
  for (auto* e : edges) {
    // check if any edge is already the edge e0 - e1
    if (e->get_twin().node() == e1) {
      return true;
    }
  }

  // TODO(philip) recover the edge
  // LOG << fmt::format("recovering edge ({}, {})", e0, e1);

  return false;
}

void OceanTriangulator::recover_edges() {
  exactinit();
  size_t n_recovered = 0;
  for (size_t i = 0; i < mesh_.lines().n(); i++) {
    auto e0 = mesh_.lines()(i, 0);
    auto e1 = mesh_.lines()(i, 1);
    // auto& node = hmesh_.nodes()[e0];
    if (e0 == e1) continue;
    if (recover_edge(hmesh_, e0, e1)) {
      n_recovered++;
    }
  }
  LOG << fmt::format("recovered {} edges (/ {})", n_recovered,
                     mesh_.lines().n());
}

bool EarClipper::triangulate(const std::vector<vec3d>& points,
                             const vec3d& normal) {
  size_t n_points = points.size();
  triangles_.clear();
  boundary_.clear();
  nodes_.clear();
  triangles_.reserve(3 * n_points);
  boundary_.reserve(triangles_.size());
  nodes_.resize(n_points);

  // first check if the polygon is convex while building the nodes
  bool convex = true;
  size_t prev = n_points - 1;
  for (size_t k = 0; k < n_points; k++) {
    nodes_[k].prev = prev;
    nodes_[k].next = (k + 1) == n_points ? 0 : k + 1;
    nodes_[k].indx = k;

    // calculate the normal
    const vec3d& a = points[prev];
    const vec3d& b = points[k];
    const vec3d& c = points[nodes_[k].next];
    vec3d nk = cross(b - a, c - a);
    if (dot(nk, normal) < 0) convex = false;

    prev = k;
  }

  if (convex) {
    for (size_t k = 2; k < n_points; k++) {
      triangles_.push_back(0);
      triangles_.push_back(k - 1);
      triangles_.push_back(k);

      bool e1 = (k + 1 == n_points);
      bool e2 = (k == 2);
      boundary_.push_back(true);
      boundary_.push_back(e1);
      boundary_.push_back(e2);
    }
    return true;
  }

  auto is_boundary = [n_points](size_t i, size_t j) {
    if (i + 1 == j) return true;
    if (i + 1 == n_points && j == 0) return true;
    return false;
  };

  auto inside = [](const vec3d& a, const vec3d& b, const vec3d& c,
                   const vec3d& d) {
    vec3d n_abc = cross(b - a, c - a);
    vec3d n_abd = cross(b - a, d - a);
    vec3d n_bcd = cross(c - b, d - b);
    vec3d n_cad = cross(a - c, d - c);
    if (dot(n_abc, n_abd) > 0 && dot(n_abc, n_bcd) > 0 && dot(n_abc, n_cad) > 0)
      return true;
    return false;
  };

  size_t n_edges = n_points;
  size_t root = 0;
  while (n_edges >= 3) {
    // find an ear
    size_t n = root;
    do {
      auto prev = nodes_[n].prev;
      auto indx = nodes_[n].indx;
      auto next = nodes_[n].next;
      const vec3d& a = points[prev];
      const vec3d& b = points[indx];
      const vec3d& c = points[next];

      // check if this triangle is inverted
      if (dot(cross(b - a, c - a), normal) < 0) {
        n = nodes_[n].next;
        continue;
      }

      // check if any other vertex is contained in the triangle
      bool ear = true;
      size_t m = nodes_[next].next;  // start with the next of next
      do {
        // check if we looped back to the prev of the current vertex
        if (m == prev) break;

        const vec3d& d = points[nodes_[m].indx];
        if (inside(a, b, c, d)) {
          ear = false;
          break;
        }
        m = nodes_[m].next;
      } while (m != root);
      if (ear) break;

      n = nodes_[n].next;
    } while (n != root);

    // found an ear
    size_t prev = nodes_[n].prev;
    size_t next = nodes_[n].next;
    triangles_.push_back(prev);
    triangles_.push_back(n);
    triangles_.push_back(next);

    boundary_.push_back(is_boundary(n, next));
    boundary_.push_back(false);
    boundary_.push_back(is_boundary(prev, n));

    // update the linked list
    nodes_[prev].next = next;
    nodes_[next].prev = prev;
    root = next;

    // we add one edge but remove two, so there is one less edge
    n_edges -= 1;
  }
  return false;
}

PolygonTriangulation::PolygonTriangulation(const Vertices& vertices,
                                           const Topology<Polygon>& polygons)
    : vertices_(vertices), polygons_(polygons) {}

void PolygonTriangulation::triangulate(TangentSpaceType type, size_t m,
                                       size_t n) {
  std::vector<vec3d> points;
  int dim = vertices_.dim();

  triangles_.reserve((n - m) * 10 * 3);
  group_.reserve(n - m);
  edge_.reserve(triangles_.capacity());

  EarClipper clipper;
  size_t n_convex = 0;
  for (size_t k = m; k < n; k++) {
    // save the points and calculate average (center)
    points.resize(polygons_.length(k));
    vec3d normal;
    vec3d center;
    const auto* pk = polygons_[k];
    for (int j = 0; j < polygons_.length(k); j++) {
      for (int d = 0; d < dim; d++) points[j][d] = vertices_[pk[j]][d];
      center = center + points[j];
    }
    center = (1.0 / points.size()) * center;

    // calculate plane normal
    if (type == TangentSpaceType::kPlanar ||
        type == TangentSpaceType::kGeneral) {
      const auto& p0 = points[0];
      for (size_t i = 2; i < points.size(); i++)
        normal = normal + cross(points[i - 1] - p0, points[i] - p0);
      normal = normalize(normal);
    } else if (type == TangentSpaceType::kSphere) {
      normal = normalize(center);
    }

    // project points to the tangent plane
    for (size_t i = 0; i < points.size(); i++)
      points[i] = points[i] - dot(points[i] - center, normal) * normal;

    // triangulate the polygon
    bool convex = clipper.triangulate(points, normal);
    if (convex) n_convex++;

    // save the result
    for (size_t i = 0; i < clipper.n_triangles(); i++) {
      for (int j = 0; j < 3; j++) {
        triangles_.push_back(pk[clipper.triangle(i)[j]]);
        edge_.push_back(clipper.boundary(i, j));
      }
      group_.push_back(k);
    }
  }
  LOG << fmt::format("detected {} convex cells out of {}", n_convex, n - m);
}

}  // namespace vortex