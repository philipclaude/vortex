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

namespace vortex {
class Mesh;
class Vertices;

template <typename coord_t>
void get_bounding_box(const coord_t* points, int64_t n_points, int8_t dim,
                      coord_t* xmin, coord_t* xmax);

template <typename coord_t, typename index_t>
void sort_points_on_zcurve(const coord_t* points, uint64_t n_points, int8_t dim,
                           std::vector<index_t>& order);

void sample_surface(const Mesh& mesh, Vertices& sites, index_t n);
void sample_volume(const Mesh& mesh, Vertices& sites, index_t n);

}  // namespace vortex