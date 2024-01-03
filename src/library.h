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

#include "mesh.h"

namespace vortex {

/**
 * \brief represents a structured grid for any element type
 */
template <typename T>
class Grid : public Mesh {
 public:
  /**
   * \brief initializes and build a structured grid
   *
   * \param[in] sizes a vector with the number of divisions in each direction
   *            for a 1d mesh (Line), sizes.size() = 1
   *            for a 2d mesh (Triangle, Quad) sizes.size() = 2
   *            for a 3d mesh (Tet), sizes.size() = 3
   * \param[in] dim - the dimension of the vertices.
   *                  Sometimes you may want to create a mesh in 3d even if
   *                  the mesh is really in 2d. When the default of -1 is used
   *                  then the ambient dimension becomes the topological
   * dimension of the element.
   */
  Grid(const std::vector<int>& sizes, int dim = -1);

  /**
   * \brief builds the structured mesh
   */
  void build();

 private:
  const std::vector<int>& sizes_;  // number of sizes in each direction
};

class Sphere : public Mesh {
 public:
  Sphere(int n = 0) : Mesh(3) { build(n); }

  void build(int n);

  static const int n_icosahedron_triangles = 20;

 private:
  void subdivide();
};

}  // namespace vortex
