#include "library.h"

#include <cmath>
#include <map>
#include <vector>

#include "llama/vec.hpp"

namespace vortex {

template <typename type>
Grid<type>::Grid(const std::vector<int>& sizes, int dim)
    : Mesh((dim < 0) ? sizes.size() : dim), sizes_(sizes) {
  build();
}

template <>
void Grid<Triangle>::build() {
  int nx = sizes_[0];
  int ny = sizes_[1];

  double dx = 1. / double(nx);
  double dy = 1. / double(ny);

  vertices_.reserve((nx + 1) * (ny + 1));
  triangles_.reserve(2 * nx * ny);

  std::vector<double> x(vertices_.dim(), 0.0);
  for (int j = 0; j < ny + 1; j++) {
    for (int i = 0; i < nx + 1; i++) {
      x[0] = i * dx;
      x[1] = j * dy;

      vertices_.add(x.data());
    }
  }

  index_t t[3];
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      int i0 = j * (nx + 1) + i;
      int i1 = i0 + 1;
      int i2 = i1 + nx + 1;
      int i3 = i2 - 1;

      t[0] = i0;
      t[1] = i1;
      t[2] = i2;
      triangles_.add(t);

      t[0] = i0;
      t[1] = i2;
      t[2] = i3;
      triangles_.add(t);
    }
  }
}

template <>
void Grid<Quad>::build() {
  // these are actually quads
  int nx = sizes_[0];
  int ny = sizes_[1];

  double dx = 1. / double(nx);
  double dy = 1. / double(ny);

  vertices_.reserve((nx + 1) * (ny + 1));
  quads_.reserve(nx * ny);

  std::vector<double> x(vertices_.dim(), 0.0);
  for (int j = 0; j < ny + 1; j++) {
    for (int i = 0; i < nx + 1; i++) {
      x[0] = i * dx;
      x[1] = j * dy;

      vertices_.add(x.data());
    }
  }

  index_t p[4];
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      int i0 = j * (nx + 1) + i;
      int i1 = i0 + 1;
      int i2 = i1 + nx + 1;
      int i3 = i2 - 1;

      p[0] = i0;
      p[1] = i1;
      p[2] = i2;
      p[3] = i3;
      quads_.add(p);
    }
  }
}

template <>
void Grid<Polygon>::build() {
  // these are actually quads
  int nx = sizes_[0];
  int ny = sizes_[1];

  double dx = 1. / double(nx);
  double dy = 1. / double(ny);

  vertices_.reserve((nx + 1) * (ny + 1));
  polygons_.reserve(nx * ny);

  std::vector<double> x(vertices_.dim(), 0.0);
  for (int j = 0; j < ny + 1; j++) {
    for (int i = 0; i < nx + 1; i++) {
      x[0] = i * dx;
      x[1] = j * dy;

      vertices_.add(x.data());
    }
  }

  index_t p[4];
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      int i0 = j * (nx + 1) + i;
      int i1 = i0 + 1;
      int i2 = i1 + nx + 1;
      int i3 = i2 - 1;

      p[0] = i0;
      p[1] = i1;
      p[2] = i2;
      p[3] = i3;
      polygons_.add(p, 4);
    }
  }
}

void Sphere::build(int n) {
  // icosahedron vertices
  coord_t t = (1.0 + std::sqrt(5.0)) / 2.0;

  coord_t coordinates[12][3] = {{t, 1, 0}, {-t, 1, 0}, {t, -1, 0}, {-t, -1, 0},
                                {1, 0, t}, {1, 0, -t}, {-1, 0, t}, {-1, 0, -t},
                                {0, t, 1}, {0, -t, 1}, {0, t, -1}, {0, -t, -1}};

  for (int i = 0; i < 12; i++) {
    for (int d = 0; d < 3; d++) coordinates[i][d] /= std::sqrt(1 + t * t);
    vertices_.add(coordinates[i]);
  }

  index_t triangles[20][3] = {
      {0, 8, 4},   // 0
      {0, 5, 10},  // 1
      {2, 4, 9},   // 2
      {2, 11, 5},  // 3
      {1, 6, 8},   // 4
      {1, 10, 7},  // 5
      {3, 9, 6},   // 6
      {3, 7, 11},  // 7
      {0, 10, 8},  // 8
      {1, 8, 10},  // 9
      {2, 9, 11},  // 10
      {3, 11, 9},  // 11
      {4, 2, 0},   // 12
      {5, 0, 2},   // 13
      {6, 1, 3},   // 14
      {7, 3, 1},   // 15
      {8, 6, 4},   // 16
      {9, 4, 6},   // 17
      {10, 5, 7},  // 18
      {11, 7, 5}   // 19
  };

  for (int i = 0; i < 20; i++) {
    triangles_.add(triangles[i]);
    triangles_.set_group(i, 0);
  }

  for (int i = 0; i < n; i++) subdivide();
}

void Sphere::subdivide() {
  std::map<Edge, index_t> edges;

  Topology<Triangle> triangles;
  for (int k = 0; k < triangles_.n(); k++) {
    int edge_indices[3];
    for (int j = 0; j < 3; j++) {
      index_t e0 = triangles_(k, j);
      index_t e1 = triangles_(k, (j + 1) % 3);

      // does this edge exist?
      index_t idx;
      auto itr = edges.find({e1, e0});
      if (itr == edges.end()) {
        // create a new point on this edge
        llama::vec3d p0(vertices_[e0]);
        llama::vec3d p1(vertices_[e1]);
        llama::vec3d q = 0.5 * (p0 + p1);

        // normalize to place point on unit sphere
        q = normalize(q);

        idx = vertices_.n();
        vertices_.add(q.data());
        edges.insert({{e0, e1}, idx});
      } else {
        idx = itr->second;
      }
      edge_indices[j] = idx;
    }

    // create the four new triangles from the subdivision
    index_t t0 = triangles_(k, 0);
    index_t t1 = triangles_(k, 1);
    index_t t2 = triangles_(k, 2);

    index_t e0 = edge_indices[0];
    index_t e1 = edge_indices[1];
    index_t e2 = edge_indices[2];

    index_t triangle0[3] = {t0, e0, e2};
    index_t triangle1[3] = {e0, t1, e1};
    index_t triangle2[3] = {e2, e1, t2};
    index_t triangle3[3] = {e0, e1, e2};

    triangles.add(triangle0);
    triangles.add(triangle1);
    triangles.add(triangle2);
    triangles.add(triangle3);
  }

  triangles_.clear(false);
  for (index_t k = 0; k < triangles.n(); k++) {
    triangles_.add(triangles[k]);
    triangles_.set_group(k, 0);
  }
}

template class Grid<Triangle>;
template class Grid<Quad>;
template class Grid<Polygon>;

}  // namespace vortex
