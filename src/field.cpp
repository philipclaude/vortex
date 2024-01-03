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
#include "field.h"

#include "mesh.h"

namespace vortex {

template <typename T>
static void allocate(const Topology<T>& topology, ElementField<T>& field) {
  ASSERT(field.n() == 0);
  ASSERT(n_basis<T>(field.order()) * field.ranks() == field.stride());
  field.set_layout(topology.layout());
  if (topology.stride() < 0) field.set_stride(-1);

  field.reserve(topology.n());
  std::vector<coord_t> data;
  for (index_t k = 0; k < topology.n(); k++) {
    data.resize(topology.length(k), 0.0);
    field.add(data.data(), data.size());
  }
}

void Field::initialize(const Mesh& mesh) {
  if (mesh.lines().n() > 0) allocate(mesh.lines(), lines_);
  if (mesh.triangles().n() > 0) allocate(mesh.triangles(), triangles_);
  if (mesh.quads().n() > 0) allocate(mesh.quads(), quads_);
  if (mesh.polygons().n() > 0) allocate(mesh.polygons(), polygons_);
}

template <typename T>
void evaluate_basis(int order, const coord_t* x, coord_t* phi);
template <typename T>
void polytope_correction(index_t n, std::vector<double>& phi) {
  return;
}

template <>
void polytope_correction<Polygon>(index_t n, std::vector<coord_t>& phi) {
  phi.resize(n, 1.0 / n);
}

template <>
void evaluate_basis<Line>(int order, const coord_t* x, coord_t* phi) {
  if (order == 0)
    phi[0] = 1.0;
  else if (order == 1) {
    phi[0] = 1.0 - x[0];
    phi[1] = x[0];
  } else
    NOT_IMPLEMENTED;
}

template <>
void evaluate_basis<Triangle>(int order, const coord_t* x, coord_t* phi) {
  if (order == 0)
    phi[0] = 1.0;
  else if (order == 1) {
    phi[0] = 1.0 - x[0] - x[1];
    phi[1] = x[0];
    phi[2] = x[1];
  } else
    NOT_IMPLEMENTED;
}

template <>
void evaluate_basis<Quad>(int order, const coord_t* x, coord_t* phi) {
  if (order == 0)
    phi[0] = 1.0;
  else if (order == 1) {
    phi[0] = (1.0 - x[0]) * (1.0 - x[1]);
    phi[1] = x[0] * (1.0 - x[1]);
    phi[2] = x[0] * x[1];
    phi[3] = (1.0 - x[0]) * x[1];
  } else
    NOT_IMPLEMENTED;
}

template <>
void evaluate_basis<Polygon>(int order, const coord_t* x, coord_t* phi) {
  if (order == 0)
    phi[0] = 1.0;
  else if (order == 1)
    return;
  else
    NOT_IMPLEMENTED;
}

template <typename T>
static void evaluate_topology(const Vertices& vertices,
                              const Topology<T>& topology,
                              ElementField<T>& field,
                              std::function<void(const double* x, double*)> f) {
  ReferenceElement<T> element;
  element.set_order(field.order());
  // using vec = typename ReferenceElement<T>::vec;

  ASSERT(n_basis<T>(field.order()) == int(element.nodes.size()));

  std::vector<coord_t> q(field.ranks());
  // for now only linear elements are supported
  std::vector<coord_t> phi(n_basis<T>(1));
  ASSERT(topology.n() == field.n());
  for (index_t k = 0; k < field.n(); k++) {
    for (size_t j = 0; j < element.nodes.size(); j++) {
      // evaluate the element basis functions at this reference coordinate
      // (only linear elements are currently supported, hence the 1)
      evaluate_basis<T>(1, element.nodes[j].data(), phi.data());
      polytope_correction<T>(topology.length(k), phi);

      // interpolate the mesh coordinates
      vortex::vec3d p;
      for (size_t i = 0; i < phi.size(); i++) {
        // vertex coordinate corresponding to this basis function
        vortex::vec3d v(vertices[topology[k][i]], vertices.dim());

        for (int d = 0; d < vertices.dim(); d++) p[d] += v[d] * phi[i];
      }

      // evaluate the function at the point
      f(p.data(), q.data());
      for (int r = 0; r < field.ranks(); r++) {
        index_t offset = r * (n_basis<T>(field.order()));
        field[k][offset + j] = q[r];
      }
    }
  }
}

void Field::evaluate(const Mesh& mesh,
                     std::function<void(const double*, double*)> f) {
  initialize(mesh);
  evaluate_topology(mesh.vertices(), mesh.lines(), lines_, f);
  evaluate_topology(mesh.vertices(), mesh.triangles(), triangles_, f);
  evaluate_topology(mesh.vertices(), mesh.quads(), quads_, f);
  evaluate_topology(mesh.vertices(), mesh.polygons(), polygons_, f);
}

template <typename T>
static void set_group(const Topology<T>& topology, ElementField<T>& field) {
  ASSERT(field.n() == topology.n()) << "must initialize field";
  for (size_t k = 0; k < field.n(); k++) {
    field[k][0] = topology.group(k);
  }
}

void FieldLibrary::set_defaults(const Mesh& mesh) {
  add("cell", 1, 0);
  auto f1 = [](const coord_t* x, coord_t* y) {
    y[0] = double(rand()) / double(RAND_MAX);
  };
  fields_["cell"].evaluate(mesh, f1);

  add("group", 1, 0);
  Field& field = fields_["group"];
  field.initialize(mesh);
  set_group(mesh.lines(), field.lines());
  set_group(mesh.triangles(), field.triangles());
  set_group(mesh.quads(), field.quads());
  set_group(mesh.polygons(), field.polygons());
}

}  // namespace vortex
