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
#include "elements.h"

#include "math/mat.hpp"
#include "math/vec.hpp"

namespace vortex {

int Quad::faces[8] = {0, 1, 1, 2, 2, 3, 3, 0};
int Triangle::faces[6] = {1, 2, 2, 0, 0, 1};
int Triangle::edges[6] = {0, 1, 1, 2, 2, 0};
const vec3d Triangle::center{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};

vec3d Triangle::get_physical_coordinates(const coord_t* pa, const coord_t* pb,
                                         const coord_t* pc, const coord_t* x) {
  vec3d point;
  for (int d = 0; d < 3; d++)
    point[d] = pa[d] * (1 - x[0] - x[1]) + pb[d] * x[0] + pc[d] * x[1];
  return point;
}

coord_t Triangle::area(const coord_t* xa, const coord_t* xb,
                       const coord_t* xc) {
  vec3d a(xa), b(xb), c(xc);
  return 0.5 * length(cross(b - a, c - a));
}

coord_t Triangle::jacobian(const coord_t* pa, const coord_t* pb,
                           const coord_t* pc, const coord_t* x) {
  vec3d a(pa), b(pb), c(pc);
  return length(cross(b - a, c - a));
}

void Triangle::get_refcoord_gradient(coord_t s, coord_t t, const coord_t* pa,
                                     const coord_t* pb, const coord_t* pc,
                                     vec3d& grads, vec3d& gradt) {
  // basis function derivatives
  vec3d dphi_ds{-1, 1, 0};
  vec3d dphi_dt{-1, 0, 1};

  // build jacobian (transpose) of transformation
  const coord_t* p[3] = {pa, pb, pc};
  vec3d dx_ds{0, 0, 0}, dx_dt{0, 0, 0};
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      dx_ds[j] += dphi_ds[i] * p[i][j];
      dx_dt[j] += dphi_dt[i] * p[i][j];
    }
  }

  vec3d n = cross(dx_ds, dx_dt);
  mats<3, 3, double> j;
  for (int i = 0; i < 3; i++) {
    j(0, i) = dx_ds[i];
    j(1, i) = dx_dt[i];
    j(2, i) = n[i];
  }

  auto jinv = inverse(j);
  for (int i = 0; i < 3; i++) {
    grads[i] = jinv(i, 0);
    gradt[i] = jinv(i, 1);
  }
}

void Triangle::get_basis_gradient(coord_t s, coord_t t, const coord_t* pa,
                                  const coord_t* pb, const coord_t* pc,
                                  vec3d& grad_fa, vec3d& grad_fb,
                                  vec3d& grad_fc) {
  vec3d grads, gradt;
  Triangle::get_refcoord_gradient(s, t, pa, pb, pc, grads, gradt);

  grad_fa = -1.0 * grads - 1.0 * gradt;
  grad_fb = grads;
  grad_fc = gradt;
}

vec3d SphericalTriangle::get_physical_coordinates(const coord_t* pa,
                                                  const coord_t* pb,
                                                  const coord_t* pc,
                                                  const coord_t* x) {
  vec3d point;
  for (int d = 0; d < 3; d++)
    point[d] = pa[d] * (1 - x[0] - x[1]) + pb[d] * x[0] + pc[d] * x[1];
  return normalize(point);
}

coord_t SphericalTriangle::area(const coord_t* xa, const coord_t* xb,
                                const coord_t* xc) {
  vec3d a(xa), b(xb), c(xc);
  coord_t num = std::fabs(dot(a, cross(b, c)));
  coord_t den = 1.0 + dot(a, b) + dot(b, c) + dot(a, c);
  return 2.0 * std::atan2(num, den);
}

coord_t SphericalTriangle::jacobian(const coord_t* pa, const coord_t* pb,
                                    const coord_t* pc, const coord_t* x) {
  vec3d x0(pa), x1(pb), x2(pc);
  coord_t s = x[0], t = x[1];
  vec3d p = (1 - s - t) * x0 + s * x1 + t * x2;
  double normp = length(p);
  vec3d dp_ds = x1 - x0;
  vec3d dp_dt = x2 - x0;

  // x = p / ||p||, where p is a function of (s, t), need dx/ds and dx/dt
  // https://math.stackexchange.com/questions/3700803/derivative-of-a-unit-vector
  mats<3, 3, double> ppt;
  for (int ii = 0; ii < 3; ii++)
    for (int jj = 0; jj < 3; jj++) ppt(ii, jj) = p[ii] * p[jj];
  vec3d dx_ds = (1.0 / normp) * dp_ds - (ppt * dp_ds) * std::pow(normp, -3);
  vec3d dx_dt = (1.0 / normp) * dp_dt - (ppt * dp_dt) * std::pow(normp, -3);

  return length(cross(dx_ds, dx_dt));
}

}  // namespace vortex