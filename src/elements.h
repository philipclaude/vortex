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

#include <array>

#include "defs.h"

namespace vortex {

template <int M, typename T>
class vecs;
typedef vecs<3, coord_t> vec3d;

typedef std::array<index_t, 2> Edge;

struct Line {
  static const int dimension = 1;
  static const int n_vertices = 2;
  static const int n_edges = 1;
  static const int n_faces = 2;
  static constexpr int edges[2] = {0, 1};
};

struct Triangle {
  static const int dimension = 2;
  static const int n_vertices = 3;
  static const int n_edges = 3;
  static const int n_faces = 3;
  static int edges[6];
  static int faces[6];
  typedef Line face_type;
  static vec3d get_physical_coordinates(const coord_t* pa, const coord_t* pb,
                                        const coord_t* pc, const coord_t* x);
  static coord_t area(const coord_t* xa, const coord_t* xb, const coord_t* xc);
  static coord_t jacobian(const coord_t* pa, const coord_t* pb,
                          const coord_t* pc, const coord_t* x);
  static void get_refcoord_gradient(coord_t s, coord_t t, const coord_t* pa,
                                    const coord_t* pb, const coord_t* pc,
                                    vec3d& grads, vec3d& gradt);

  static void get_basis(coord_t s, coord_t t, double* basis);
  static void get_basis_gradient(coord_t s, coord_t t, const coord_t* pa,
                                 const coord_t* pb, const coord_t* pc,
                                 vec3d& grad_fa, vec3d& grad_fb,
                                 vec3d& grad_fc);
  static const vec3d center;
};

struct SphericalTriangle {
  static const int dimension = 2;  // topological dimension
  static const int n_vertices = 3;
  static const int n_edges = 3;
  static const int n_faces = 3;
  static int edges[6];
  static int faces[6];
  typedef Line face_type;
  static vec3d get_physical_coordinates(const coord_t* pa, const coord_t* pb,
                                        const coord_t* pc, const coord_t* x);
  static coord_t area(const coord_t* xa, const coord_t* xb, const coord_t* xc);
  static coord_t jacobian(const coord_t* pa, const coord_t* pb,
                          const coord_t* pc, const coord_t* x);
};

struct Quad {
  static const int dimension = 2;
  static const int n_vertices = 4;
  static const int n_edges = 4;
  static const int n_faces = 4;
  static constexpr int edges[8] = {0, 1, 1, 2, 2, 3, 3, 0};
  static int faces[8];
  typedef Line face_type;
};

struct Polygon {
  static const int dimension = 2;
  static const int n_vertices = -1;
  static const int n_edges = -1;
  static const int n_faces = -1;
  static constexpr int* edges = nullptr;
  typedef Line face_type;
};

struct Icosahedron {
  static const int n_faces = 20;
  static const int n_vertices = 12;
  static double coordinates[n_vertices][3];
  static int faces[n_faces][3];
};

struct Octahedron {
  static const int n_faces = 8;
  static const int n_vertices = 6;
  static double coordinates[n_vertices][3];
  static int faces[n_faces][3];
};

}  // namespace vortex
