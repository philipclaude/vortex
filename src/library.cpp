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
#include "library.h"

#include <cmath>
#include <map>
#include <unordered_map>
#include <vector>

#include "math/vec.hpp"
#include "stlext.h"

namespace vortex {

template <typename type>
Grid<type>::Grid(const std::vector<int>& sizes) : Mesh(3), sizes_(sizes) {
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

  double x[3] = {0, 0, 0};
  for (int j = 0; j < ny + 1; j++) {
    for (int i = 0; i < nx + 1; i++) {
      x[0] = i * dx;
      x[1] = j * dy;

      vertices_.add(x);
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

  double x[3] = {0, 0, 0};
  for (int j = 0; j < ny + 1; j++) {
    for (int i = 0; i < nx + 1; i++) {
      x[0] = i * dx;
      x[1] = j * dy;

      vertices_.add(x);
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

  double x[3] = {0, 0, 0};
  for (int j = 0; j < ny + 1; j++) {
    for (int i = 0; i < nx + 1; i++) {
      x[0] = i * dx;
      x[1] = j * dy;

      vertices_.add(x);
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

template <typename T>
void SubdividedSphere<T>::build(int n) {
  for (int i = 0; i < T::n_vertices; i++) {
    vertices_.add(T::coordinates[i]);
  }

  for (int i = 0; i < T::n_faces; i++) {
    triangles_.add(T::faces[i]);
    triangles_.set_group(i, 0);
  }

  for (int i = 0; i < n; i++) subdivide();
}

template <typename T>
void SubdividedSphere<T>::subdivide() {
  std::map<Edge, index_t> edges;

  Topology<Triangle> triangles;
  for (size_t k = 0; k < triangles_.n(); k++) {
    int edge_indices[3];
    for (int j = 0; j < 3; j++) {
      index_t e0 = triangles_(k, j);
      index_t e1 = triangles_(k, (j + 1) % 3);

      // does this edge exist?
      index_t idx;
      auto itr = edges.find({e1, e0});
      if (itr == edges.end()) {
        // create a new point on this edge
        vortex::vec3d p0(vertices_[e0]);
        vortex::vec3d p1(vertices_[e1]);
        vortex::vec3d q = 0.5 * (p0 + p1);

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

void Squircle::build(double R, int nr, int nt, bool half) {
  ASSERT(R < 1);
  double dr = (1.0 - R) / double(nr);
  double dt = M_PI / double(nt);

  if (half) {
    vertices_.reserve((nr + 1) * (nt + 1));
    triangles_.reserve(2 * nr * nt);
  } else {
    vertices_.reserve(2 * (nr + 1) * (nt + 1));
    triangles_.reserve(4 * nr * nt);
  }

  const double tr2 = 2.0 * std::sqrt(2);

  double x[3] = {0, 0, 0};
  for (int j = 0; j < nt + 1; j++) {
    for (int i = 0; i < nr + 1; i++) {
      double r = R + std::pow(i * dr, 1);
      double u = r * cos(j * dt);
      double v = r * sin(j * dt);

      x[0] = 0.5 * std::sqrt(std::fabs(2 + tr2 * u + u * u - v * v)) -
             0.5 * std::sqrt(std::fabs(2 - tr2 * u + u * u - v * v));
      x[1] = 0.5 * std::sqrt(std::fabs(2 + tr2 * v - u * u + v * v)) -
             0.5 * std::sqrt(std::fabs(2 - tr2 * v - u * u + v * v));
      vertices_.add(x);
    }
  }

  if (!half) {
    for (int j = 1; j < nt; j++) {
      for (int i = 0; i < nr + 1; i++) {
        double r = R + std::pow(i * dr, 1);
        double u = r * cos(M_PI + j * dt);
        double v = r * sin(M_PI + j * dt);

        x[0] = 0.5 * std::sqrt(std::fabs(2 + tr2 * u + u * u - v * v)) -
               0.5 * std::sqrt(std::fabs(2 - tr2 * u + u * u - v * v));
        x[1] = 0.5 * std::sqrt(std::fabs(2 + tr2 * v - u * u + v * v)) -
               0.5 * std::sqrt(std::fabs(2 - tr2 * v - u * u + v * v));
        vertices_.add(x);
      }
    }
  }

  index_t t[3];
  for (int j = 0; j < nt; j++) {
    for (int i = 0; i < nr; i++) {
      int i0 = j * (nr + 1) + i;
      int i1 = i0 + 1;
      int i2 = i1 + nr + 1;
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

  if (!half) {
    for (int j = nt; j < 2 * nt; j++) {
      for (int i = 0; i < nr; i++) {
        int i0 = j * (nr + 1) + i;
        int i1 = i0 + 1;
        int i2 = i1 + nr + 1;
        int i3 = i2 - 1;

        if (j + 1 == 2 * nt) {
          i2 = i + 1;
          i3 = i;
        }

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

  // extract boundary edges
  std::unordered_map<std::pair<int, int>, int> edges;
  for (size_t k = 0; k < triangles_.n(); k++) {
    for (int j = 0; j < 3; j++) {
      int p = triangles_[k][j];
      int q = triangles_[k][j == 2 ? 0 : j + 1];
      auto it = edges.find({q, p});
      if (it == edges.end()) {
        edges.insert({{p, q}, 2});
      } else
        edges.erase(it);
    }
  }

  // tag boundaries
  lines_.reserve(edges.size());
  for (const auto& [e, _] : edges) {
    vec3d e0(vertices_[e.first]);
    vec3d e1(vertices_[e.second]);
    if (length(0.5 * (e0 + e1)) < 1.5 * R) continue;  // on interior circle
    vec3d n = cross(e1 - e0, {0, 0, 1});
    if (std::fabs(n[0]) < 1e-6) continue;  // upper/lower boundaries
    if (n[0] > 0)
      edges[e] = 1;
    else
      edges[e] = 3;
  }

  for (auto& [e, bnd] : edges) {
    int edge[2] = {e.first, e.second};
    size_t nl = lines_.n();
    lines_.add(edge);
    lines_.set_group(nl, bnd);
  }
}

template class Grid<Triangle>;
template class Grid<Quad>;
template class Grid<Polygon>;
template class SubdividedSphere<Icosahedron>;
template class SubdividedSphere<Octahedron>;

}  // namespace vortex
