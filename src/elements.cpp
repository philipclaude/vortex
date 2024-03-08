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

#include "math/vec.hpp"

namespace vortex {

int Quad::faces[8] = {0, 1, 1, 2, 2, 3, 3, 0};
int Triangle::faces[6] = {1, 2, 2, 0, 0, 1};
int Triangle::edges[6] = {0, 1, 1, 2, 2, 0};

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

}  // namespace vortex