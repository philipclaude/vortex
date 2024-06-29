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
#include "graphics.h"

#include "io.h"
#include "library.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(graphics_test_suite)

UT_TEST_CASE_SKIP(test1) {
  // SubdividedIcosahedron mesh(0);
  //    mesh.vertices().print();
  //    mesh.triangles().print();
  //  Grid<Triangle> mesh({10, 10}, 3);
  //  Mesh mesh(3);
  //  meshb::read("water.meshb", mesh);
  Grid<Quad> mesh({10, 10});
  mesh.fields().set_defaults(mesh);
#if VORTEX_FULL_UNIT_TEST == 0
  ParticleAnimationParameters particle_params;
  Viewer viewer(mesh, particle_params, 7681);
#endif
}
UT_TEST_CASE_END(test1)

UT_TEST_CASE(animation_test) {
#if VORTEX_FULL_UNIT_TEST == 0
  Mesh mesh(3);  // empty mesh
  ParticleAnimationParameters particle_params;
  particle_params.points_prefix = "../../build/release/sphere100k/particles";
  Viewer viewer(mesh, particle_params, 7681);
#endif
}
UT_TEST_CASE_END(animation_test)

UT_TEST_SUITE_END(graphics_test_suite)