#include "library.h"

#include "io.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(library_test_suite)

UT_TEST_CASE(grid_triangle_test) {}
UT_TEST_CASE_END(grid_triangle_test)

UT_TEST_CASE(polygon_test) {
  int n = 10;
  Grid<Polygon> mesh({n, n});
  meshb::write(mesh, "polygrid.meshb");
}
UT_TEST_CASE_END(polygon_test)

UT_TEST_CASE(sphere_test) {
  Sphere mesh(2);
  meshb::write(mesh, "sphere.meshb");
}
UT_TEST_CASE_END(sphere_test)

UT_TEST_SUITE_END(library_test_suite)
