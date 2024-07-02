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
 * \brief Represents a structured grid for any element type
 */
template <typename T>
class Grid : public Mesh {
 public:
  /**
   * \brief Initializes and builds a structured grid
   *
   * \param[in] sizes a vector with the number of divisions in each direction
   *            for a 1d mesh (Line), sizes.size() = 1
   *            for a 2d mesh (Triangle, Quad) sizes.size() = 2
   *            for a 3d mesh (Tet), sizes.size() = 3
   */
  Grid(const std::vector<int>& sizes);

  /**
   * \brief Builds the structured mesh
   */
  void build();

 private:
  const std::vector<int>& sizes_;  // number of sizes in each direction
};

template <typename S>
class SubdividedSphere : public Mesh {
 public:
  SubdividedSphere(int n = 0) : Mesh(3) { build(n); }

 private:
  void build(int n);
  void subdivide();
};

/// @brief Represents a subdivided icosahedron mesh.
class SubdividedIcosahedron : public SubdividedSphere<Icosahedron> {
 public:
  /// @brief Constructs a subdivided icosahedron mesh.
  /// @param n number of subdivisions.
  SubdividedIcosahedron(int n = 0) : SubdividedSphere<Icosahedron>(n) {}
};

/// @brief Represents a circle mapped to a square.
/// http://squircular.blogspot.com/2015/09/mapping-circle-to-square.html
class Squircle : public Mesh {
 public:
  /// @brief Constructs the circle-in-square mesh.
  /// @param r radius of the inner circle.
  /// @param nr number of divisions in the radial direction.
  /// @param nt number of divisions in theta direction from 0 - PI.
  /// @param half whether to create a domain in [0, PI] (true) or [0, 2 * PI].
  Squircle(double r, int nr, int nt, bool half = false) : Mesh(3) {
    build(r, nr, nt, half);
  }

 private:
  void build(double r, int nr, int nt, bool half);
};

}  // namespace vortex
