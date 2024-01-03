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
#include "triangulate.h"

#include <unordered_set>

#include "io.h"
#include "library.h"
#include "mesh.h"
#include "stlext.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(triangulate_tests)

UT_TEST_CASE(test1) {
  // input mesh written by mesher_test
  // TODO(philip) ensure mesher_test is always run first
  Mesh input_mesh(3);
  read_mesh("water.meshb", input_mesh);

  // extract the lines
  std::unordered_set<std::array<index_t, 2>> edges;
  for (size_t k = 0; k < input_mesh.triangles().n(); k++) {
    const auto* t = input_mesh.triangles()[k];
    for (int i = 0; i < 3; i++) {
      int j = (i + 1) == 3 ? 0 : i + 1;
      std::array<index_t, 2> e = {t[i], t[j]};
      if (e[0] > e[1]) std::swap(e[0], e[1]);
      if (edges.find(e) == edges.end()) {
        edges.insert(e);
        input_mesh.lines().add(e.data());
      }
    }
  }

  // only keep vertices on the lines
  Mesh coast(3);
  std::unordered_map<index_t, index_t> vertex_map;
  vertex_map.reserve(input_mesh.vertices().n());
  for (size_t k = 0; k < input_mesh.lines().n(); k++) {
    auto* e = input_mesh.lines()[k];
    for (int j = 0; j < 2; j++) {
      if (vertex_map.find(e[j]) == vertex_map.end()) {
        vertex_map.insert({e[j], vertex_map.size()});
        coast.vertices().add(input_mesh.vertices()[e[j]]);
      }
      e[j] = vertex_map.at(e[j]);
    }
  }

  Sphere oceans(4);
  OceanTriangulator triangulator(oceans, coast);
  triangulator.triangulate();

  meshb::write(oceans, "oceans.meshb");
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(triangulate_tests)