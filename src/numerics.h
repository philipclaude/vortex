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
#pragma once
#include <algorithm>

#include "log.h"
#include "math/mat.hpp"
#include "math/vec.hpp"
#include "predicates.h"

namespace vortex {

using mat3 = vortex::mats<3, 3, double>;
using vec2d = vortex::vecs<2, double>;

/// @brief Determines the parameter space coordinates (u, v) from a point on the
/// sphere centered at the origin with a radius of 1.
/// @param xyz Coordinates of the point on the sphere.
/// @param uv Parameter space coordinates (theta, phi) divided by (2\pi, \pi).
inline void sphere_params(const vec3d& xyz, vec3d& uv) {
  uv[0] = 0.5 * (atan2(xyz[1], xyz[0]) + M_PI) / M_PI;
  uv[1] = 1.0 - acos(xyz[2]) / M_PI;
  uv[2] = 0.0;
  ASSERT(uv[0] >= 0 && uv[0] <= 1);
  ASSERT(uv[1] >= 0 && uv[1] <= 1);
}

inline bool get_params(const vec3d& pa, const vec3d& pb, const vec3d& pc,
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

inline bool get_params(const vec3d& pa, const vec3d& pb, const vec3d& pc,
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

inline double face_area(const double* xa, const double* xb, const double* xc) {
  vec3d pa(xa);
  vec3d pb(xb);
  vec3d pc(xc);
  vec3d ua, ub, uc;
  get_params(pa, pb, pc, ua, ub, uc);
  return 0.5 * orient2d(&ua[0], &ub[0], &uc[0]);
}

inline double spherical_triangle_area(const vec3d& a, const vec3d& b,
                                      const vec3d& c) {
  double num = std::fabs(dot(a, cross(b, c)));
  double den = 1.0 + dot(a, b) + dot(b, c) + dot(a, c);
  return 2.0 * std::atan2(num, den);
}

}  // namespace vortex