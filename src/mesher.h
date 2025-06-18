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

#include <memory>

namespace vortex {

class Texture;
class Mesh;
class HalfMesh;

struct MeshingParameters {
  double h_min{0.025};  // minimum size in the mesh
  double h_max{0.05};   // maximum size in the mesh
  int max_iter{5};      // number of adaptation iterations
};

/// @brief Build a mesh of a sphere from an image (texture) to determine the
/// sizes that drive the mesh adaptation.
class EarthMesher {
 public:
  /// @brief Initializes the mesher, saving the texture image.
  /// @param texture image that determines the sizes used for the continents and
  /// oceans.
  EarthMesher(const Texture& texture);
  ~EarthMesher();

  /// @brief Generates a mesh of the Earth.
  /// @param params See the MeshingParameters above.
  void generate(MeshingParameters params);

  /// @brief Accesses a reference to the working half-edge mesh data structure.
  auto& mesh() { return *mesh_; }

  /// @brief Extracts the mesh from the working half-edge mesh data structure to
  /// an array-based mesh structure.
  /// @param mesh Destination of the mesh.
  void extract(Mesh& mesh) const;

 private:
  std::unique_ptr<HalfMesh> mesh_;
  const Texture& texture_;
};

}  // namespace vortex