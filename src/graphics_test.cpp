#include "graphics.h"

#include "io.h"
#include "library.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(graphics_test_suite)

UT_TEST_CASE(test1) {
  // Sphere mesh(3);
  //  mesh.vertices().print();
  //  mesh.triangles().print();
  // Grid<Polygon> mesh({10, 10}, 3);
  Mesh mesh(3);
  meshb::read("water.meshb", mesh);
  // Grid<Quad> mesh({10, 10}, 3);
  mesh.fields().set_defaults(mesh);
  Viewer viewer(mesh, 7681, true);
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(graphics_test_suite)