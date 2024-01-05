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

  meshb::write(mesh, "clip.meshb");
}
UT_TEST_CASE_END(ear_clipping_test)

UT_TEST_CASE(polygon_triangulation_test) {
  Mesh oceans(3);
  // meshb::read("water_voronoi.meshb", oceans);
  meshb::read("test.meshb", oceans);

  // merge vertices
  oceans.merge();
  LOG << "done merging vertices";

  Mesh mesh(3);
  oceans.vertices().copy(mesh.vertices());
  LOG << "# of polygons = " << oceans.polygons().n();

  // count number of groups
  int n_groups = 0;
  for (size_t k = 0; k < oceans.polygons().n(); k++)
    n_groups = std::max(oceans.polygons().group(k), n_groups);
  n_groups++;  // account for zero-indexing
  LOG << "detected " << n_groups << " groups";

  // combine polygons with the same group (Voronoi cell)
  std::vector<std::vector<index_t>> cell2polygon(n_groups);
  for (auto& cell : cell2polygon) cell.reserve(10);
  for (size_t k = 0; k < oceans.polygons().n(); k++) {
    int group = oceans.polygons().group(k);
    // LOG << fmt::format("poly {} in group {}", k, group);
    cell2polygon[group].push_back(k);
  }

  // merge polygons with the same group
  std::vector<index_t> polygon;
  std::unordered_map<std::array<index_t, 2>, std::array<int64_t, 2>> edges;
  std::vector<int64_t> components;
  std::queue<index_t> queue;
  for (size_t k = 0; k < cell2polygon.size(); k++) {
    auto& polygons = cell2polygon[k];
    if (polygons.empty()) continue;

    // build a list of edges to left & right polygons
    edges.clear();
    for (int64_t j = 0; j < polygons.size(); j++) {
      // loop through all vertices (edges) of this polygon
      auto p = polygons[j];
      auto n = oceans.polygons().length(p);
      for (size_t i = 0; i < n; i++) {
        index_t a = oceans.polygons()[p][i];
        index_t b = oceans.polygons()[p][(i + 1) == n ? 0 : i + 1];
        auto eit = edges.find({b, a});
        if (eit != edges.end()) {
          // edge has already been visited in opposite direction
          eit->second[1] = j;
        } else
          edges.insert({{a, b}, {j, -1}});  // right polygon is null (-1)
      }
    }

// identify connected components of this cell
#if 1
    components.resize(polygons.size());
    for (auto& c : components) c = -1;
    int component = 0;
    while (true) {
      // find a polygon
      // LOG << "building component " << component << " with " <<
      // polygons.size();
      int64_t p0 = -1;
      for (size_t j = 0; j < polygons.size(); j++) {
        if (components[j] < 0) {
          p0 = j;
          break;
        }
      }
      if (p0 < 0) {
        // LOG << "done";
        break;  // done
      }

      // BFS starting with p0
      queue.push(p0);
      while (!queue.empty()) {
        index_t j = queue.front();
        queue.pop();
        components[j] = component;

        // loop through all the edges
        auto p = polygons[j];
        auto n = oceans.polygons().length(p);
        for (size_t i = 0; i < n; i++) {
          index_t a = oceans.polygons()[p][i];
          index_t b = oceans.polygons()[p][(i + 1) == n ? 0 : i + 1];
          auto eit = edges.find({a, b});
          int64_t neighbor;
          if (eit == edges.end()) {
            auto tie = edges.find({b, a});
            ASSERT(tie != edges.end());
            ASSERT(tie->second[1] == j);
            neighbor = tie->second[0];
          } else {
            ASSERT(eit->second[0] == j);
            neighbor = eit->second[1];
          }

          if (components[neighbor] < 0) {
            components[neighbor] = component;
            queue.push(neighbor);
          }
        }
      }
      component++;
    }
    ASSERT(component > 0);
    if (component > 1) LOG << fmt::format("detected {} components", component);
#endif

    // rebuild the polygon
    polygon.clear();

    // find an edge to start with
    index_t p, q;
    for (const auto& [edge, faces] : edges) {
      if (faces[1] >= 0) continue;
      p = edge[0];
      q = edge[1];
      break;
    }

    polygon.push_back(p);
    index_t root = p;
    size_t n_visited = 0;
    do {
      n_visited++;
      bool found = false;
      index_t next;
      for (const auto& [edge, faces] : edges) {
        if (faces[1] >= 0) continue;
        if (edge[0] == q) {
          found = true;
          next = edge[1];
          break;
        }
      }
      ASSERT(found);
      polygon.push_back(q);
      p = q;
      q = next;
      if (n_visited == edges.size()) break;
    } while (q != root);
    if (n_visited == edges.size()) {
      LOG << "skip";
    }
    // std::reverse(polygon.begin(), polygon.end());
    mesh.polygons().add(polygon.data(), polygon.size());
  }

  PolygonTriangulationThread triangulation(mesh.vertices(), mesh.polygons());
  triangulation.triangulate(TangentSpaceType::kSphere, 0, mesh.polygons().n());

  for (size_t k = 0; k < triangulation.n(); k++) {
    mesh.triangles().add(triangulation.triangle(k));
    mesh.triangles().set_group(k, triangulation.group(k));
  }

  // mesh.polygons().clear();

  oceans.triangles().clear();
  meshb::write(oceans, "test2.meshb");

  meshb::write(mesh, "earclip.meshb");
}
UT_TEST_CASE_END(polygon_triangulation_test)

UT_TEST_SUITE_END(triangulate_tests)