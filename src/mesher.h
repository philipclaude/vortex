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

#include <memory>

namespace vortex {

class Texture;
class Mesh;
class HalfMesh;

struct MeshingParameters {
  double h_min{0.025};
  double h_max{0.05};
  int max_iter{5};
};

class EarthMesher {
 public:
  EarthMesher(const Texture& texture);
  ~EarthMesher();
  void generate(MeshingParameters params);
  auto& mesh() { return *mesh_; }
  void extract(Mesh& mesh) const;

 private:
  std::unique_ptr<HalfMesh> mesh_;
  const Texture& texture_;
};

}  // namespace vortex