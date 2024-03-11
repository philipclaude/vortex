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

#include <vector>

#include "defs.h"
#include "elements.h"
#include "math/vec.h"

namespace vortex {

template <typename T>
class TriangleQuadrature {
 public:
  static const int dim = T::dimension;
  using point_t = vecs<dim, coord_t>;

  TriangleQuadrature(int order) { define(order); }

  template <typename Integrand>
  coord_t integrate(const Integrand& fn, const coord_t* pa, const coord_t* pb,
                    const coord_t* pc) const;

  const auto& points() const { return points_; }
  const auto& weights() const { return weights_; }
  void use_constant_jacobian(bool x) { use_constant_jacobian_ = x; }

 private:
  bool use_constant_jacobian_{true};
  void define(int order);
  std::vector<point_t> points_;
  std::vector<coord_t> weights_;
};

}  // namespace vortex