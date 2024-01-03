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
#include "numerics.h"

#include "halfedges.h"
#include "io.h"
#include "library.h"
#include "predicates.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(numerics_test_suite)

double face_area(const HalfFace& face) {
  const auto& e0 = face.get_edge();
  const auto& e1 = e0.get_next();
  const auto& e2 = e1.get_next();
  const auto& n0 = e0.get_node();
  const auto& n1 = e1.get_node();
  const auto& n2 = e2.get_node();
  return vortex::face_area(n0.point(), n1.point(), n2.point());
}

UT_TEST_CASE(test1) {
  exactinit();

  Sphere sphere(4);
  meshb::write(sphere, "sphere.meshb");

  double area = 0.0;
  size_t n_periodic = 0;
  for (index_t k = 0; k < sphere.triangles().n(); k++) {
    const auto* t = sphere.triangles()[k];
    vec3d p0(sphere.vertices()[t[0]]);
    vec3d p1(sphere.vertices()[t[1]]);
    vec3d p2(sphere.vertices()[t[2]]);

    vec3d u0, u1, u2;
    if (!get_params(p0, p1, p2, u0, u1, u2)) n_periodic++;

    double ak = 0.5 * orient2d(&u0[0], &u1[0], &u2[0]);

    if (ak < 0) {
      LOG << fmt::format("area = {}", ak);
      LOG << fmt::format("p0 = ({}, {}, {})", p0[0], p0[1], p0[2]);
      LOG << fmt::format("p1 = ({}, {}, {})", p1[0], p1[1], p1[2]);
      LOG << fmt::format("p2 = ({}, {}, {})", p2[0], p2[1], p2[2]);
    }

    UT_ASSERT(ak > 0);
    area += ak;
  }
  LOG << fmt::format("area = {}, n_periodic = {}", area, n_periodic);
}
UT_TEST_CASE_END(test1)

UT_TEST_CASE_SKIP(test2) {
  exactinit();

  Sphere sphere(0);
  HalfMesh mesh(sphere);

  for (auto& f : mesh.faces()) {
    double area = face_area(f);
    LOG << fmt::format("area {} = {}", f.index(), area);
    // UT_ASSERT(area > 0);
  }

  double area = 0.0;
  size_t n_periodic = 0;
  for (auto& e : mesh.edges()) {
    auto& p = e.get_node();
    auto& q = e.get_twin().get_node();
    auto& m = e.get_triangle_left_node();
    auto& n = e.get_triangle_right_node();

    vec3d xp(p.point());
    vec3d xq(q.point());
    vec3d xm(m.point());
    vec3d xn(n.point());

    vec3d up, uq, um, un;
    if (!get_params(xp, xq, xm, xn, up, uq, um, un)) n_periodic++;

    double a0 = 0.5 * orient2d(&up[0], &uq[0], &um[0]);
    UT_ASSERT(a0 > 0);

    double a1 = 0.5 * orient2d(&up[0], &un[0], &uq[0]);
    UT_ASSERT(a1 > 0);
    area += a0 + a1;
  }
  // every triangle is counted 6 times (3 for each edge x 2 for twins)
  LOG << fmt::format("area = {}, n_periodic = {}", area, n_periodic);
}
UT_TEST_CASE_END(test2)

UT_TEST_SUITE_END(numerics_test_suite)