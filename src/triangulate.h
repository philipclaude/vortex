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

#include "halfedges.h"

namespace vortex {

class Mesh;

class OceanTriangulator {
 public:
  struct Options {};
  OceanTriangulator(Mesh& mesh, const Mesh& coast);
  void triangulate();

 private:
  void insert_points();
  void recover_edges();

  Mesh& mesh_;
  const Mesh& coast_;
  HalfMesh hmesh_;
};

}  // namespace vortex