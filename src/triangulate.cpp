#include "triangulate.h"

#include <mapletrees/kdtree.h>

#include <unordered_set>

#include "halfedges.h"
#include "library.h"
#include "mesh.h"
#include "numerics.h"
#include "predicates.h"

namespace vortex {

OceanTriangulator::OceanTriangulator(Mesh& mesh, const Mesh& coast)
    : mesh_(mesh), coast_(coast) {}

void OceanTriangulator::triangulate() {
  insert_points();
  recover_edges();
}

using vec2d = llama::vecs<2, double>;

void sphere_params(const vec3d& xyz, vec3d& uv) {
  constexpr double tol = 1e-12;
  uv[0] = 0.5 * (atan2(xyz[1], xyz[0]) + M_PI) / M_PI;
  if (std::fabs(xyz[0]) < tol && std::fabs(xyz[1]) < tol) uv[0] = 0.0;
  if (std::fabs(uv[0] - 1.0) < tol) uv[0] = 0.0;
  uv[1] = (M_PI - acos(xyz[2])) / M_PI;
  uv[2] = 0.0;
  ASSERT(uv[0] >= 0 && uv[0] <= 1);
  ASSERT(uv[1] >= 0 && uv[1] <= 1);
}

void triangle_params(const vec3d& pa, const vec3d& pb, const vec3d& pc,
                     vec3d& pd, vec3d& ua, vec3d& ub, vec3d& uc, vec3d& ud,
                     bool verbose = false) {
  sphere_params(pa, ua);
  sphere_params(pb, ub);
  sphere_params(pc, uc);
  sphere_params(pd, ud);

  // check the bounding box in the u-direction
  double umin = std::min(ua[0], std::min(ub[0], uc[0]));
  double umax = std::max(ua[0], std::max(ub[0], uc[0]));
  if (verbose) {
    LOG << fmt::format("checking query (u, v) = ({}, {}), urange = [{}, {}]",
                       ud[0], ud[1], umin, umax);
  }
  if (umax - umin > 0.5) {
    // LOG << fmt::format("spans period: umin = {}, umax = {}", umin, umax);
    //  the triangle points cross the periodic boundary, adjust them
    if (ua[0] > 0.5) ua[0] -= 1;
    if (ub[0] > 0.5) ub[0] -= 1;
    if (uc[0] > 0.5) uc[0] -= 1;
    // check if the query point needs to be adjusted
    if (ud[0] > umax) {
      ud[0] -= 1;
    }
    umin = std::min(ua[0], std::min(ub[0], uc[0]));
    umax = std::max(ua[0], std::max(ub[0], uc[0]));
    if (verbose) LOG << fmt::format("range = [{}, {}]", umin, umax);
  }
}

bool in_triangle(const Vertices& vertices, const HalfMesh::HalfFace* face,
                 const double* x, bool verbose = false) {
  ASSERT(face);
  // extract 3d coordinates of triangle and query point
  vec3d pa(vertices[face->edge->node->index], 3);
  vec3d pb(vertices[face->edge->next->node->index], 3);
  vec3d pc(vertices[face->edge->next->next->node->index], 3);
  vec3d pd(x, 3);

  // check distance to triangle vertices is less than radius
  // if (length(pa - pd) > 1) return false;
  // if (length(pb - pd) > 1) return false;
  // if (length(pc - pd) > 1) return false;

  // determine uv coordinates on the sphere, accounting for periodicity
  vec3d ua, ub, uc, ud;
  triangle_params(pa, pb, pc, pd, ua, ub, uc, ud, verbose);

  if (verbose) {
    LOG << fmt::format("ua = ({}, {})", ua[0], ua[1]);
    LOG << fmt::format("ub = ({}, {})", ub[0], ub[1]);
    LOG << fmt::format("uc = ({}, {})", uc[0], uc[1]);
    LOG << fmt::format("ud = ({}, {})", ud[0], ud[1]);
  }

  // calculate barycentric coordinates
  vec3d nt = orient2d(&ua[0], &ub[0], &uc[0]);
  vec3d na = orient2d(&ub[0], &uc[0], &ud[0]);
  vec3d nb = orient2d(&ud[0], &uc[0], &ua[0]);
  vec3d nc = orient2d(&ua[0], &ub[0], &ud[0]);
  double at = nt[2];
  double ka = na[2] / at;
  double kb = nb[2] / at;
  double kc = nc[2] / at;

  return (ka >= 0 && ka <= 1 && kb >= 0 && kb <= 1 && kc >= 0 && kc <= 1);
}

void OceanTriangulator::insert_points() {
  // build a kdtree of the surface points
  maple::KdTree<3, coord_t, index_t, true> tree(mesh_.vertices()[0],
                                                mesh_.vertices().n());

  // build a halfedge representation of the surface mesh
  HalfMesh hmesh(mesh_);

  std::unordered_set<HalfMesh::HalfFace*> visited;
  visited.reserve(hmesh.faces().size());

  exactinit();

  // for each point in the coast
  std::vector<HalfMesh::HalfFace*> faces;
  for (int i = 0; i < coast_.vertices().n(); i++) {
    index_t n = tree.nearest(coast_.vertices()[i]);
    LOG << fmt::format("closest = {} to vertex {}", n, i);

    // pick one face to start the search
    auto& node = hmesh.nodes()[n];
    // faces.clear();
    // hmesh.get_onering(&node, faces);
    HalfMesh::HalfFace* triangle = nullptr;
    for (auto& face : hmesh.faces()) {
      // LOG << fmt::format("checking face {} out of", face.index);
      if (in_triangle(mesh_.vertices(), &face, coast_.vertices()[i])) {
        triangle = &face;
        break;
      }
    }
    if (!triangle) {
      for (auto& face : hmesh.faces()) {
        LOG << fmt::format("checking triangle {}", face.index);
        in_triangle(mesh_.vertices(), &face, coast_.vertices()[i], true);
      }
    }
    ASSERT(triangle);
    LOG << fmt::format("found face {}", triangle->index);
    // hmesh.insert(triangle, coast_.vertices()[i]);
  }
}

void OceanTriangulator::recover_edges() {}

}  // namespace vortex