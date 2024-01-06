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
  int max_components = 1, n_multiple_components = 0;
  using edge_t = std::array<index_t, 2>;
  using face_pair_t = std::array<int64_t, 2>;
  std::vector<index_t> polygon;
  std::unordered_map<edge_t, face_pair_t> edges;
  std::vector<int> components;
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

    // initialize component of each polygon to null (-1)
    components.resize(polygons.size());
    for (auto& c : components) c = -1;
    int component = 0;

    // identify connected components of this cell
    while (true) {
      // find an unset polygon (with a null component) to start with
      int64_t root = -1;
      for (size_t j = 0; j < polygons.size(); j++) {
        if (components[j] < 0) {
          root = j;
          break;
        }
      }
      if (root < 0) break;  // done, no more unset polygons

      // BFS starting with the first polygon index (root)
      queue.push(root);
      while (!queue.empty()) {
        // next polygon index in the queue
        index_t j = queue.front();
        queue.pop();
        components[j] = component;

        // loop through all the edges of this polygon
        auto p = polygons[j];
        auto n = oceans.polygons().length(p);
        for (size_t i = 0; i < n; i++) {
          index_t a = oceans.polygons()[p][i];
          index_t b = oceans.polygons()[p][(i + 1) == n ? 0 : i + 1];
          auto eit = edges.find({a, b});
          int64_t neighbor;
          if (eit == edges.end()) {
            // look for the reversed edge
            auto tie = edges.find({b, a});
            ASSERT(tie != edges.end());
            ASSERT(tie->second[1] == j);
            neighbor = tie->second[0];
          } else {
            ASSERT(eit->second[0] == j);
            neighbor = eit->second[1];
          }

          if (neighbor < 0) continue;  // boundary
          if (components[neighbor] < 0) {
            // assign the component and add this neighbor to the queue
            components[neighbor] = component;
            queue.push(neighbor);
          }
        }
      }
      component++;
    }
    ASSERT(component > 0);
    max_components = std::max(component, max_components);
    if (component > 1) n_multiple_components++;

    // create polygons for each component
    for (int c = 0; c < component; c++) {
      // find an edge to start with that has this component on the left
      // and is on the boundary
      index_t p, q;
      bool found = false;
      for (const auto& [edge, faces] : edges) {
        auto fl = faces[0];
        auto fr = faces[1];
        if (fr >= 0) continue;  // not on the boundary
        auto cl = components[fl];
        if (cl != c) continue;  // not in this component
        p = edge[0];
        q = edge[1];
        found = true;
      }
      ASSERT(found);

      polygon.clear();
      polygon.push_back(p);

      index_t root = p;
      do {
        found = false;
        index_t next;
        for (const auto& [edge, faces] : edges) {
          auto fl = faces[0];
          auto fr = faces[1];
          if (fr >= 0) continue;  // not on the boundary
          auto cl = components[fl];
          if (cl != c) continue;  // not in this component
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
      } while (q != root);
      mesh.polygons().add(polygon.data(), polygon.size());
    }
  }
  LOG << fmt::format("maximum # components = {}, # multiple components = {}",
                     max_components, n_multiple_components);

  PolygonTriangulationThread triangulation(mesh.vertices(), mesh.polygons());
  triangulation.triangulate(TangentSpaceType::kSphere, 0, mesh.polygons().n());

  for (size_t k = 0; k < triangulation.n(); k++) {
    mesh.triangles().add(triangulation.triangle(k));
    mesh.triangles().set_group(k, triangulation.group(k));
  }

  oceans.triangles().clear();
  meshb::write(oceans, "test2.meshb");
  meshb::write(mesh, "earclip.meshb");
}
UT_TEST_CASE_END(polygon_triangulation_test)

UT_TEST_SUITE_END(triangulate_tests)