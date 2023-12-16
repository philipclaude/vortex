#pragma once

#include "mesh.h"

namespace vortex {

/**
 * \brief represents a structured grid for any element type
 */
template <typename T> class Grid : public Mesh {
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
