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
#include "mesher.h"

#include <stlext.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>

#include "halfedges.h"
#include "library.h"
#include "numerics.h"
#include "texture.h"

namespace vortex {

namespace {

double size(const Texture& texture, const MeshingParameters& opts, double* p) {
  double tp;
  vec3d vp(p), up;
  sphere_params(vp, up);
  texture.sample(up[0], up[1], &tp);
  tp /= 255.;
  return opts.h_min + tp * (opts.h_max - opts.h_min);
}

double size(const Texture& texture, const MeshingParameters& opts, double* p,
            double* q) {
  // size at p and q
  double hp = size(texture, opts, p);
  double hq = size(texture, opts, q);

  // calculate the interpolated size
  double r, hm;
  if (hp > hq) {
    r = hp / hq;
    hm = hp;
  } else {
    r = hq / hp;
    hm = hq;
  }
  if (std::fabs(r - 1.0) < 1e-12) return hp;
  return hm * (r - 1.0) / (r * std::log(r));
}

bool triangle_points_out(const vec3d& a, const vec3d& b, const vec3d& c) {
  vec3d n = normalize(cross(b - a, c - a));
  double tol = 1e-3;  // todo pass allowed dot product based on current mesh
  if (dot(a, n) < tol) return false;
  if (dot(b, n) < tol) return false;
  if (dot(c, n) < tol) return false;
  return true;
}

bool valid_collapse(const HalfMesh& mesh, const HalfEdge& edge) {
  // check if the collapse of p onto q is valid
  auto& p = edge.get_node();
  auto& q = edge.get_twin().get_node();
  vec3d x(q.point());
  half_t fl = edge.face();
  half_t fr = edge.get_twin().face();

  half_t e = edge.index();
  half_t first = e;
  do {
    const HalfFace& f = mesh.edges()[e].get_face();
    if (f.index() == fl || f.index() == fr) continue;

    // extract the coordinates of the face
    auto& e0 = f.get_edge();
    auto& e1 = e0.get_next();
    auto& e2 = e1.get_next();

    vec3d p0(e0.get_node().point());
    vec3d p1(e1.get_node().point());
    vec3d p2(e2.get_node().point());

    if (e0.node() == p.index()) {
      if (!triangle_points_out(x, p1, p2)) return false;
    }
    if (e1.node() == p.index()) {
      if (!triangle_points_out(p0, x, p2)) return false;
    }
    if (e2.node() == p.index()) {
      if (!triangle_points_out(p0, p1, x)) return false;
    }
    e = mesh.edges()[e].get_twin().next();
  } while (e != first);

  return true;
}

void collapse(HalfMesh& mesh, const Texture& texture,
              const MeshingParameters& opts) {
  // metric edge length
  auto metric_length = [&mesh, &texture, &opts](const half_t iedge) {
    auto& edge = mesh.edges()[iedge];
    ASSERT(edge.node() != edge.get_twin().node());
    vec3d p(edge.get_node().point());
    vec3d q(edge.get_twin().get_node().point());
    double l = length(q - p);
    ASSERT(l > 0);
    double h = size(texture, opts, &p[0], &q[0]);
    ASSERT(h > 0);
    return l / h;
  };

  // insert short edges into the queue
  const double lmin = 0.5 * std::sqrt(2.0);
  auto cmp = [&metric_length](half_t i, half_t j) {
    return metric_length(i) < metric_length(j);
  };
  std::vector<half_t> queue;
  std::vector<bool> visited;
  std::vector<HalfNode*> np, nq;

  bool collapsed = true;
  size_t n_collapsed = 0;
  size_t n_rejected = 0;
  size_t n_nonmanifold = 0;
  int pass = 0;
  do {
    pass++;
    collapsed = false;
    visited.resize(mesh.edges().size(), false);
    queue.clear();
    for (auto& e : mesh.edges()) {
      if (!e.active()) continue;
      if (!e.get_twin().active()) continue;
      double l = metric_length(e.index());
      if (l < lmin) queue.push_back(e.index());
    }
    std::sort(queue.begin(), queue.end(), cmp);

    for (auto& iedge : queue) {
      // retrieve edge and twin information
      HalfEdge& edge = mesh.edges()[iedge];
      HalfEdge& twin = edge.get_twin();
      if (!edge.active()) continue;
      if (!twin.active()) continue;
      ASSERT(iedge == edge.index());
      ASSERT(edge.node() != twin.node());

      // check if we have already visited this edge, or if the length changed
      if (visited[edge.index()]) continue;
      if (visited[twin.index()]) continue;
      visited[edge.index()] = true;
      visited[twin.index()] = true;
      double l = metric_length(iedge);
      if (l > lmin) continue;

      // the number of common vertices between p and q must be 2
      // this is an O(N^2) loop but N is small so it's fine
      edge.get_node().get_onering(np);
      twin.get_node().get_onering(nq);
      int n_common = 0;
      for (auto* u : np) {
        for (auto* v : nq) {
          if (u->index() == v->index()) n_common++;
        }
      }
      if (n_common != 2) {
        n_nonmanifold++;
        continue;
      }

      // try to collapse this edge
      if (!valid_collapse(mesh, edge)) {
        n_rejected++;
        continue;
      }
      mesh.collapse(iedge);
      n_collapsed++;
      collapsed = true;
    }
  } while (collapsed);
  LOG << fmt::format(
      "--> collapsed {} in {} passes, rejected {}, # nonmanifold = {} ",
      n_collapsed, pass, n_rejected, n_nonmanifold);
}

bool split(HalfMesh& mesh, const Texture& texture,
           const MeshingParameters& opts, double lmax) {
  // metric edge length
  auto metric_length = [&mesh, &texture, &opts](const half_t iedge) {
    auto& edge = mesh.edges()[iedge];
    auto& twin = edge.get_twin();
    vec3d p(edge.get_node().point());
    vec3d q(edge.get_twin().get_node().point());
    ASSERT(edge.node() != twin.node());
    double l = length(q - p);
    ASSERT(l > 0) << fmt::format("l = {}", l);
    double h = size(texture, opts, &p[0], &q[0]);
    ASSERT(h > 0);
    return l / h;
  };

  // insert long edges into the queue
  double lmin = 1e-8;
  auto cmp = [&metric_length](half_t i, half_t j) {
    return metric_length(i) > metric_length(j);
  };
  std::vector<half_t> queue;

  bool inserted = true;
  size_t n_inserted = 0;
  int pass = 0;
  do {
    pass++;
    inserted = false;
    queue.clear();
    for (auto& e : mesh.edges()) {
      if (!e.active()) continue;
      if (!e.get_twin().active()) continue;
      double l = metric_length(e.index());
      if (l > lmax) queue.push_back(e.index());
    }
    std::sort(queue.begin(), queue.end(), cmp);

    for (auto& iedge : queue) {
      HalfEdge& edge = mesh.edges()[iedge];
      HalfEdge& twin = edge.get_twin();
      if (!edge.active()) continue;
      if (!twin.active()) continue;
      ASSERT(iedge == edge.index()) << iedge << ", " << edge.index();
      double l = metric_length(iedge);
      if (l < lmax) continue;

      ASSERT(edge.node() != twin.node());

      // coordinates of point to insert
      vec3d p(edge.get_node().point());
      vec3d q(twin.get_node().point());
      ASSERT(length(p - q) > lmin);
      vec3d x = normalize(0.5 * (p + q));

      // TODO check that points are not too close
      vec3d xl(edge.get_triangle_left_node().point());
      vec3d xr(edge.get_triangle_right_node().point());

      double hl = length(xl - x);  // / size(texture, opts, &xl[0], &x[0]);
      double hr = length(xr - x);  // / size(texture, opts, &xr[0], &x[0]);
      if (hl < lmin || hr < lmin) continue;

      // insert the point on the edge
      mesh.split(iedge, &x[0]);
      //(mesh.check());

      n_inserted++;
      inserted = true;
    }
  } while (inserted);
  LOG << fmt::format("--> inserted {} points in {} passes", n_inserted, pass);
  return true;
}

double quality(const vec3d& pa, const vec3d& pb, const vec3d& pc) {
  double area = spherical_triangle_area(pa, pb, pc);
  double lab = length(pb - pa);
  double lbc = length(pc - pb);
  double lac = length(pc - pa);
  return 4.0 * std::sqrt(3.0) * area / (lab * lab + lbc * lbc + lac * lac);
}

double quality(const HalfFace& face) {
  const auto& e0 = face.get_edge();
  const auto& e1 = e0.get_next();
  const auto& e2 = e1.get_next();
  vec3d pa(e0.get_node().point());
  vec3d pb(e1.get_node().point());
  vec3d pc(e2.get_node().point());
  return quality(pa, pb, pc);
}

bool valid_flip(const HalfEdge& edge) {
  const HalfEdge& twin = edge.get_twin();

  // extract 3d coordinates of edge endpoints
  vec3d pa(edge.get_node().point());
  vec3d pb(twin.get_node().point());

  // extract 3d coordinates of left & right vertices
  ASSERT(edge.get_triangle_left_node().index() !=
         edge.get_triangle_right_node().index());
  vec3d pl(edge.get_triangle_left_node().point());
  vec3d pr(edge.get_triangle_right_node().point());

  auto normal = [](const vec3d& a, const vec3d& b, const vec3d& c) {
    return normalize(cross(b - a, c - a));
  };

  double q0 = std::min(quality(pa, pb, pl), quality(pa, pr, pb));
  double q1 = std::min(quality(pl, pr, pb), quality(pa, pr, pl));
  if (q1 < q0) return false;  // quality should improve

  // double d0 = dot(normal(pa, pb, pl), normal(pa, pr, pb));
  double d1 = dot(normal(pl, pr, pb), normal(pa, pr, pl));
  if (d1 < 0.1) return false;

  return true;
}

void flips(HalfMesh& mesh, const Texture& texture) {
  auto q_edge = [&mesh](half_t e) -> double {
    const auto& edge = mesh.edges()[e];
    return std::min(quality(edge.get_face()),
                    quality(edge.get_twin().get_face()));
  };

  // create initial list of faces
  auto cmp = [&q_edge](half_t a, half_t b) { return q_edge(a) < q_edge(b); };
  std::vector<half_t> queue;
  std::vector<bool> visited;
  std::vector<HalfFace*> faces;

  double qmin = 0.8;

  bool flipped = true;
  size_t n_flipped = 0;
  size_t n_rejected = 0;
  size_t n_nonmanifold = 0;
  int pass = 0;
  do {
    pass++;
    flipped = false;
    visited.resize(mesh.edges().size(), false);
    queue.clear();
    double qmax = -1;
    for (auto& e : mesh.edges()) {
      if (!e.active()) continue;
      double q = q_edge(e.index());
      if (q > qmax) qmax = q;
      if (q < qmin) queue.push_back(e.index());
    }
    std::sort(queue.begin(), queue.end(), cmp);

    for (auto& iedge : queue) {
      HalfEdge& edge = mesh.edges()[iedge];
      HalfEdge& twin = edge.get_twin();
      if (!edge.active()) continue;
      if (!twin.active()) continue;
      ASSERT(iedge == edge.index()) << iedge << ", " << edge.index();
      if (visited[edge.index()]) continue;
      if (visited[twin.index()]) continue;
      visited[edge.index()] = true;
      visited[twin.index()] = true;
      double q = q_edge(iedge);
      if (q >= qmin) continue;

      // check the topology of the flip
      edge.get_node().get_onering(faces);
      if (faces.size() <= 3) {
        n_nonmanifold++;
        continue;
      }
      twin.get_node().get_onering(faces);
      if (faces.size() <= 3) {
        n_nonmanifold++;
        continue;
      }

      // try to flip this edge
      if (!valid_flip(edge)) {
        n_rejected++;
        continue;
      }
      mesh.flip(iedge);

      n_flipped++;
      flipped = true;
    }
  } while (flipped);
  LOG << fmt::format(
      "--> flipped {} in {} passes, rejected {}, # nonmanifold = {}", n_flipped,
      pass, n_rejected, n_nonmanifold);
}

void smooth(HalfMesh& mesh, int max_iter) {
  int n_threads = std::thread::hardware_concurrency();
  std::vector<index_t> n_smoothed(n_threads, 0);
  std::vector<double> distance(n_threads, 0);
  std::vector<double> coordinates(mesh.nodes().size() * 3);
  for (int iter = 0; iter < max_iter; iter++) {
    std::parafor_i(
        0, mesh.nodes().size(),
        [&mesh, &n_smoothed, &distance, &coordinates](int tid, int inode) {
          HalfNode& n = mesh.nodes()[inode];
          if (!n.active()) return;
          n_smoothed[tid]++;
          double area = 0;
          half_t e = n.edge();
          half_t first = e;
          vec3d centroid;
          const double alpha = 0.0;
          vec3d p(n.point());
          vec3d moment;
          do {
            HalfFace& f = mesh.edges()[e].get_face();
            const auto& e0 = f.get_edge();
            const auto& e1 = e0.get_next();
            const auto& e2 = e1.get_next();
            vec3d pa(e0.get_node().point());
            vec3d pb(e1.get_node().point());
            vec3d pc(e2.get_node().point());
            double at = spherical_triangle_area(pa, pb, pc);
            vec3d c = normalize((1.0 / 3.0) * (pa + pb + pc));
            moment = moment + at * c;
            area += at;
            e = mesh.edges()[e].get_twin().next();
          } while (e != first);
          centroid = normalize(moment / area);
          distance[tid] += length(centroid - p);
          for (int d = 0; d < 3; d++)
            coordinates[n.index() * 3 + d] =
                alpha * p[d] + (1 - alpha) * centroid[d];
        });
    std::parafor_i(0, mesh.nodes().size(),
                   [&mesh, coordinates](int tid, int inode) {
                     if (!mesh.nodes()[inode].active()) return;
                     for (int d = 0; d < 3; d++)
                       mesh.vertices()[inode][d] = coordinates[3 * inode + d];
                   });
  }
  double total_d = std::accumulate(distance.begin(), distance.end(), 0.0);
  double total_n = std::accumulate(n_smoothed.begin(), n_smoothed.end(), 0);
  LOG << fmt::format(
      "--> smoothed {} points with {} threads, displacement = {:.1f}",
      total_n / max_iter, n_threads, total_d);
}

}  // namespace

EarthMesher::EarthMesher(const Texture& texture) : texture_(texture) {}

void EarthMesher::generate(MeshingParameters params) {
  // estimate the number of levels from the mean of the min and max sizes
  double h_avg = params.h_max;  // 0.5 * (params.h_min + params.h_max);
  double at = std::sqrt(3.0) * h_avg * h_avg / 4.0;
  int n_triangles = std::floor(4.0 * M_PI / at);
  int n_level =
      int(std::log(n_triangles / Icosahedron::n_faces) / std::log(4.0));
  LOG << fmt::format("initializing sphere mesh with {} subdivisions", n_level);
  SubdividedIcosahedron sphere(n_level);

  mesh_ = std::make_unique<HalfMesh>(sphere);
  auto& mesh = *mesh_;

  for (int iter = 1; iter <= params.max_iter; iter++) {
    size_t n_nodes = std::count_if(mesh.nodes().begin(), mesh.nodes().end(),
                                   [](const auto& n) { return n.active(); });
    size_t n_faces = std::count_if(mesh.faces().begin(), mesh.faces().end(),
                                   [](const auto& f) { return f.active(); });
    LOG << fmt::format("adaptation {}: # nodes = {}, # faces = {}", iter,
                       n_nodes, n_faces);

    // collapse short edges
    collapse(mesh, texture_, params);
    if (!mesh.check()) break;
    flips(mesh, texture_);

    // split long edges without creating short edges
    split(mesh, texture_, params, 2.0);
    split(mesh, texture_, params, std::sqrt(2.0));
    ASSERT(mesh.check());

    // optimize by flips
    flips(mesh, texture_);
    ASSERT(mesh.check());

    // optimize by smoothing
    smooth(mesh, 10);
    // if (iter == params.max_iter) smooth(mesh, 100);
  }
}

void EarthMesher::extract(Mesh& mesh) const { mesh_->extract(mesh); }

EarthMesher::~EarthMesher() {}

}  // namespace vortex