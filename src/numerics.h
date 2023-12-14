#include <algorithm>
#include <llama/mat.hpp>
#include <llama/vec.hpp>

#include "log.h"
#include "predicates.h"

namespace vortex {

using vec3d = llama::vec3d;
using mat3 = llama::mats<3, 3, double>;
using vec2d = llama::vecs<2, double>;

static void sphere_params(const vec3d& xyz, vec3d& uv) {
  constexpr double tol = 1e-12;
  uv[0] = 0.5 * (atan2(xyz[1], xyz[0]) + M_PI) / M_PI;
  if (std::fabs(xyz[0]) < tol && std::fabs(xyz[1]) < tol) uv[0] = 0.0;
  if (std::fabs(uv[0] - 1.0) < tol) uv[0] = 0.0;
  uv[1] = (M_PI - acos(xyz[2])) / M_PI;
  uv[2] = 0.0;
  ASSERT(uv[0] >= 0 && uv[0] <= 1);
  ASSERT(uv[1] >= 0 && uv[1] <= 1);
}

static bool get_params(const vec3d& pa, const vec3d& pb, const vec3d& pc,
                       vec3d& pd, vec3d& ua, vec3d& ub, vec3d& uc, vec3d& ud) {
  sphere_params(pa, ua);
  sphere_params(pb, ub);
  sphere_params(pc, uc);
  sphere_params(pd, ud);

  // check the bounding box in the u-direction
  double umin = std::min(ua[0], std::min(ub[0], std::min(uc[0], ud[0])));
  double umax = std::max(ua[0], std::max(ub[0], std::max(uc[0], ud[0])));
  if (umax - umin > 0.5) {
    //  the triangle points cross the periodic boundary, adjust them
    if (ua[0] > 0.5) ua[0] -= 1;
    if (ub[0] > 0.5) ub[0] -= 1;
    if (uc[0] > 0.5) uc[0] -= 1;
    if (ud[0] > 0.5) ud[0] -= 1;
    return false;
  }
  return true;
}

static bool get_params(const vec3d& pa, const vec3d& pb, const vec3d& pc,
                       vec3d& ua, vec3d& ub, vec3d& uc) {
  sphere_params(pa, ua);
  sphere_params(pb, ub);
  sphere_params(pc, uc);

  // check the bounding box in the u-direction
  double umin = std::min(ua[0], std::min(ub[0], uc[0]));
  double umax = std::max(ua[0], std::max(ub[0], uc[0]));
  if (umax - umin > 0.5) {
    //  the triangle points cross the periodic boundary, adjust them
    if (ua[0] > 0.5) ua[0] -= 1;
    if (ub[0] > 0.5) ub[0] -= 1;
    if (uc[0] > 0.5) uc[0] -= 1;
    return false;
  }
  return true;
}

static double face_area(const double* xa, const double* xb, const double* xc) {
  vec3d pa(xa);
  vec3d pb(xb);
  vec3d pc(xc);
  vec3d ua, ub, uc;
  get_params(pa, pb, pc, ua, ub, uc);
  return 0.5 * orient2d(&ua[0], &ub[0], &uc[0]);
}

}  // namespace vortex