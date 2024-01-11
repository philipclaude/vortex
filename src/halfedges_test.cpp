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
#include "halfedges.h"

#include "io.h"
#include "library.h"
#include "numerics.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(halfedges_test_suite)

UT_TEST_CASE(test_closed) {
  SubdividedIcosahedron mesh(3);
  HalfMesh hmesh(mesh);
  UT_ASSERT(hmesh.check());
}
UT_TEST_CASE_END(test_closed)

UT_TEST_CASE(test_open) {
  Grid<Triangle> mesh({10, 10}, 3);
  HalfMesh hmesh(mesh);
  UT_ASSERT(hmesh.check());

  Mesh m(3);
  hmesh.extract(m);
  meshb::write(m, "grid.meshb");
}
UT_TEST_CASE_END(test_open)

UT_TEST_CASE(insert_test) {
  SubdividedIcosahedron sphere(6);
  HalfMesh hmesh(sphere);

  size_t n_faces = hmesh.faces().size();
  LOG << n_faces;
  for (size_t i = 0; i < n_faces; i++) {
    double x[3] = {0, 0, 0};
    for (int j = 0; j < 3; j++)
      for (int d = 0; d < 3; d++)
        x[d] += sphere.vertices()[sphere.triangles()(i, j)][d] / 3.0;
    hmesh.insert(i, x);
  }
  UT_ASSERT(hmesh.check());

  Mesh mesh(3);
  // hmesh.extract(mesh);
  //  meshb::write(mesh, "insert.meshb");
}
UT_TEST_CASE_END(insert_test)

UT_TEST_CASE(split_test) {
  SubdividedIcosahedron sphere(5);
  HalfMesh hmesh(sphere);

  size_t n_edges = hmesh.edges().size();
  for (size_t i = 0; i < n_edges; i++) {
    auto& e = hmesh.edges()[i];
    vec3d p(e.get_node().point());
    vec3d q(e.get_twin().get_node().point());
    vec3d m = 0.5 * (p + q);
    hmesh.split(i, &m[0]);
  }
  UT_ASSERT(hmesh.check());

  Mesh mesh(3);
  hmesh.extract(mesh);
  meshb::write(mesh, "split.meshb");
}
UT_TEST_CASE_END(split_test)

UT_TEST_SUITE_END(halfedges_test_suite)
