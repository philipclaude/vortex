#pragma once

#include <array>

#include "defs.h"

namespace vortex {

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

}  // namespace vortex
