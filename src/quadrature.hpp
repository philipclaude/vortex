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
#include "quadrature.h"

namespace vortex {

template <typename T>
template <typename Integrand>
coord_t TriangleQuadrature<T>::integrate(const Integrand& fn, const coord_t* pa,
                                         const coord_t* pb,
                                         const coord_t* pc) const {
  coord_t integral = 0.0;
  double da = 2 * T::area(pa, pb, pc);
  for (size_t k = 0; k < points_.size(); k++) {
    // evaluate the point in physical space
    vec3d point = T::get_physical_coordinates(pa, pb, pc, &points_[k][0]);

    // evaluate the determinant of the jacobian at the reference point
    if (!use_constant_jacobian_) da = T::jacobian(pa, pb, pc, &points_[k][0]);

    // add the contribution to the integral
    integral += weights_[k] * fn(point) * da;
  }

  return integral;
}

}  // namespace vortex