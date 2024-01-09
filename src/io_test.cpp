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
#include "io.h"

#include "library.h"
#include "mesh.h"
#include "numerics.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(io_tests)

UT_TEST_CASE(obj_test) {
  SubdividedIcosahedron sphere(4);
  obj::write(sphere, "sphere.obj");

  Mesh mesh(3);
  read_mesh("sphere.obj", mesh);
  UT_ASSERT_EQUALS(mesh.vertices().dim(), 3);
  UT_ASSERT_EQUALS(mesh.vertices().n(), sphere.vertices().n());
  UT_ASSERT_EQUALS(mesh.triangles().n(), sphere.triangles().n());
}
UT_TEST_CASE_END(obj_test)

UT_TEST_CASE(meshb_test) {
  SubdividedIcosahedron sphere(4);
  meshb::write(sphere, "sphere.meshb");

  Mesh mesh(3);
  read_mesh("sphere.meshb", mesh);
  UT_ASSERT_EQUALS(mesh.vertices().dim(), 3);
  UT_ASSERT_EQUALS(mesh.vertices().n(), sphere.vertices().n());
  UT_ASSERT_EQUALS(mesh.triangles().n(), sphere.triangles().n());
}
UT_TEST_CASE_END(meshb_test)

UT_TEST_SUITE_END(io_tests)