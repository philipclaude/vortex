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
#include "library.h"

#include <math.h>

#include <limits>
#include <type_traits>

#include "elements.h"
#include "io.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(library_test_suite)

UT_TEST_CASE(grid_triangle_test) {
  // testing number of vertices and triangles in a triange mesh
  double tol = 1e-10;
  for (int i = 1; i < 6; i++) {
    for (int j = 1; j < 6; j++) {
      int n1t = i * 100;
      int n2t = j * 100;
      Grid<Triangle> mesh({n1t, n2t});
      int vertst = mesh.vertices().n();
      int trit = mesh.triangles().n();
      UT_ASSERT_EQUALS(vertst, ((n1t + 1) * (n2t + 1)));
      UT_ASSERT_EQUALS(trit, (2 * n1t * n2t));
      // testing sum area of each triangle
      double tot_areat = 0;
      for (int k = 0; k < trit; k++) {
        auto* t = mesh.triangles()[k];
        const auto* p1t = mesh.vertices()[t[0]];
        const auto* p2t = mesh.vertices()[t[1]];
        const auto* p3t = mesh.vertices()[t[2]];
        tot_areat += Triangle::area(p1t, p2t, p3t);
      }
      // total area should be very close to 1
      UT_ASSERT_NEAR(1., tot_areat, tol);
    }
  }
}
UT_TEST_CASE_END(grid_triangle_test)

UT_TEST_CASE(grid_quad_test) {
  // testing the number of vertices and quads in a quad mesh
  double tol = 1e-10;
  for (int i = 1; i < 6; i++) {
    for (int j = 1; j < 6; j++) {
      int n1q = i * 100;
      int n2q = j * 100;
      Grid<Quad> mesh({n1q, n2q});
      int vertsq = mesh.vertices().n();
      int quads = mesh.quads().n();
      UT_ASSERT_EQUALS(vertsq, ((n1q + 1) * (n2q + 1)));
      UT_ASSERT_EQUALS(quads, (n1q * n2q));
      // testing sum area of each quad
      double tot_areaq = 0;
      for (int k = 0; k < quads; k++) {
        auto* t = mesh.quads()[k];
        const auto* p1q = mesh.vertices()[t[0]];
        const auto* p2q = mesh.vertices()[t[1]];
        const auto* p3q = mesh.vertices()[t[2]];
        const auto* p4q = mesh.vertices()[t[0]];
        const auto* p5q = mesh.vertices()[t[2]];
        const auto* p6q = mesh.vertices()[t[3]];
        tot_areaq += Triangle::area(p1q, p2q, p3q);
        tot_areaq += Triangle::area(p4q, p5q, p6q);
      }
      // total area should be very close to 1
      UT_ASSERT_NEAR(tot_areaq, 1., tol);
    }
  }
}
UT_TEST_CASE_END(grid_quad_test)

UT_TEST_CASE(grid_polygon_test) {
  // testing the number of vertices and polygons in a polygon mesh
  double tol = 1e-10;
  for (int i = 1; i < 6; i++) {
    for (int j = 1; j < 6; j++) {
      int nx = i * 100;
      int ny = j * 100;
      Grid<Polygon> mesh({nx, ny});
      auto& polygons = mesh.polygons();
      auto& vertices = mesh.vertices();
      UT_ASSERT_EQUALS(vertices.n(), (nx + 1) * (ny + 1));
      UT_ASSERT_EQUALS(polygons.n(), nx * ny);
      // testing sum area of each polygon (actually quad here)
      double total_area = 0.0;
      for (int k = 0; k < polygons.n(); k++) {
        const auto& polygon = polygons[k];
        const auto& vertex_indices = polygons.length(k);
        for (size_t v = 1; v < vertex_indices - 1; ++v) {
          const coord_t* p0 = vertices[polygon[0]];
          const coord_t* p1 = vertices[polygon[v]];
          const coord_t* p2 = vertices[polygon[v + 1]];
          total_area += Triangle::area(p0, p1, p2);
        }
      }
      // total area should be very close to 1
      UT_ASSERT_NEAR(total_area, 1., tol);
    }
  }
}

UT_TEST_CASE_END(grid_polygon_test)

UT_TEST_CASE(sphere_test) {
  // testing number of triangles in a subdivided sphere mesh
  double tol1 = 1e-10;
  double tol2 = 5e-3;
  auto test_shape = [&](auto& shape) {
    using Shape_t = typename std::remove_reference<decltype(shape)>::type;
    // vectors for area check
    std::vector<double> error;
    std::vector<double> meshsize;
    for (int i = 0; i < 7; i++) {
      SubdividedSphere<Shape_t> mesh(i);
      int tris = mesh.triangles().n();
      int tris_check = std::pow(4, i) * Shape_t::n_faces;
      UT_ASSERT_EQUALS(tris, tris_check);
      double tot_areas = 0;
      double strt_areas = 0;
      for (int k = 0; k < tris; k++) {
        const auto* t = mesh.triangles()[k];
        const coord_t* p1 = mesh.vertices()[t[0]];
        const coord_t* p2 = mesh.vertices()[t[1]];
        const coord_t* p3 = mesh.vertices()[t[2]];
        tot_areas += SphericalTriangle::area(p1, p2, p3);
        strt_areas += Triangle::area(p1, p2, p3);
      }
      // recording straight-sided area error
      error.push_back(fabs(strt_areas - (4 * M_PI)));
      meshsize.push_back(sqrt(tris));
      //  abs. val. of straight-sided area against meshsize (approx
      //  sqrt(mesh.triangles().n())) should be very close to 2 in the
      //  asymptotic range
      if (i > 4) {
        double slope = fabs(
            log(error[error.size() - 2] / error[error.size() - 1]) /
            log(meshsize[meshsize.size() - 2] / meshsize[meshsize.size() - 1]));
        UT_ASSERT_NEAR(slope, 2., tol2);
      }
      // spherical triangle area should converge to 4pi
      UT_ASSERT_NEAR(tot_areas, 4 * M_PI, tol1);
    }
  };
  Icosahedron icosahedron;
  test_shape(icosahedron);
  Octahedron octahedron;
  test_shape(octahedron);
}
UT_TEST_CASE_END(sphere_test)

UT_TEST_CASE(squircle_test) {
  auto get_area = [](const Mesh& m) -> double {
    double area = 0;
    for (int k = 0; k < m.triangles().n(); k++) {
      const auto* t = m.triangles()[k];
      const auto* p1 = m.vertices()[t[0]];
      const auto* p2 = m.vertices()[t[1]];
      const auto* p3 = m.vertices()[t[2]];
      area += Triangle::area(p1, p2, p3);
    }
    return area;
  };

  int n = 50;
  double r = 0.1;
  int nr = n;
  int nt = 2 * n;
  Squircle mesh(r, nr, nt);
  meshb::write(mesh, "circle_square.meshb");
  UT_ASSERT_NEAR(get_area(mesh), 4.0 - M_PI * r * r, 1e-4);

  Squircle mesh_half(r, nr, nt, true);
  meshb::write(mesh_half, "circle_square_half.meshb");
  UT_ASSERT_NEAR(get_area(mesh_half), 2.0 - 0.5 * M_PI * r * r, 1e-4);

  UT_ASSERT_EQUALS(mesh.triangles().n(), 4 * nr * nt);
  UT_ASSERT_EQUALS(mesh.triangles().n(), 2 * mesh_half.triangles().n());
  UT_ASSERT_EQUALS(mesh.lines().n(), 4 * nt);
  UT_ASSERT_EQUALS(mesh_half.lines().n(), 2 * nt + 2 * nr);
  UT_CATCH_EXCEPTION(Squircle(1.1, 10, 10));
}
UT_TEST_CASE_END(squircle_test)

UT_TEST_SUITE_END(library_test_suite)
