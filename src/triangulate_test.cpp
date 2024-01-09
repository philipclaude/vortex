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

#include <trees/kdtree.h>

#include <unordered_set>

#include "halfedges.h"
#include "io.h"
#include "library.h"
#include "math/vec.hpp"
#include "mesh.h"
#include "stlext.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(triangulate_tests)

UT_TEST_CASE(test1) {
  // input mesh written by dev/examples.sh
  Mesh input_mesh(3);
  read_mesh("oceans_coarse.meshb", input_mesh);

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

  meshb::write(oceans, "inserted_coast.meshb");
}
UT_TEST_CASE_END(test1)

UT_TEST_CASE(ear_clipping_test) {
  std::vector<vec3d> points;
  points.push_back({0, 0, 0});
  points.push_back({1, 0, 0});
  points.push_back({0.75, 0.5, 0.0});
  points.push_back({0.75, 0.6, 0.0});
  points.push_back({1, 1, 0});
  points.push_back({0, 1, 0});
  points.push_back({0.74, 0.5, 0.0});
  EarClipper clipper;
  clipper.triangulate(points, {0, 0, 1});

  Mesh mesh(3);
  for (size_t k = 0; k < points.size(); k++) mesh.vertices().add(&points[k][0]);
  for (size_t k = 0; k < clipper.n_triangles(); k++) {
    mesh.triangles().add(clipper.triangle(k));
  }

  std::vector<index_t> polygon(points.size());
  for (size_t k = 0; k < polygon.size(); k++) polygon[k] = k;
  mesh.polygons().add(polygon.data(), polygon.size());

  meshb::write(mesh, "clip.meshb");
}
UT_TEST_CASE_END(ear_clipping_test)

UT_TEST_CASE(polygon_triangulation_test) {
  Mesh oceans(3);
  // meshb::read("oceans_voronoi_coarse.meshb", oceans);
  meshb::read("voronoi_oceans_coarse.meshb", oceans);

  // merge vertices
  Timer timer;
  timer.start();
  oceans.merge();
  timer.stop();
  LOG << fmt::format("done merging vertices in {} seconds", timer.seconds());

  // combine polygons with same group, then separate into connected components
  Mesh mesh(3);
  oceans.vertices().copy(mesh.vertices());
  LOG << "# of polygons = " << oceans.polygons().n();
  timer.start();
  oceans.separate_polygons_into_connected_components(mesh.polygons());
  timer.stop();
  LOG << fmt::format("found polygon connected components in {} seconds",
                     timer.seconds());

  // triangulate the polygons
  timer.start();
  PolygonTriangulation triangulation(mesh.vertices(), mesh.polygons());
  triangulation.triangulate(TangentSpaceType::kSphere, 0, mesh.polygons().n());
  timer.stop();
  LOG << fmt::format("triangulated polygons in {} seconds", timer.seconds());

  // save triangles
  for (size_t k = 0; k < triangulation.n(); k++) {
    mesh.triangles().add(triangulation.triangle(k));
    mesh.triangles().set_group(k, triangulation.group(k));
  }

  meshb::write(mesh, "earclip.meshb");
}
UT_TEST_CASE_END(polygon_triangulation_test)

UT_TEST_SUITE_END(triangulate_tests)