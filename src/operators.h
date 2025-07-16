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
#pragma once

#include <vector>

#include "defs.h"

namespace vortex {

class VoronoiDiagram;

enum OperatorMethod : uint8_t { kKincl, kSpringel };

template <typename Domain_t>
class VoronoiOperators {
 public:
  VoronoiOperators(const VoronoiDiagram& voronoi);

  void set_boundary_value(double x) { boundary_value_ = x; }
  void set_project(bool x) { project_ = x; }
  void set_method(OperatorMethod m) { method_ = m; }
  void calculate_gradient(const coord_t* f, coord_t* grad_f);
  void calculate_divergence(const coord_t* u, coord_t* div_u);
  void calculate_curl(const coord_t* u, coord_t* w);

 private:
  const VoronoiDiagram& voronoi_;
  double boundary_value_{1e20};
  bool project_{true};
  OperatorMethod method_{OperatorMethod::kSpringel};
};

}  // namespace vortex