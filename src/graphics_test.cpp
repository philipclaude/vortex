#include "graphics.h"

#include "library.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(graphics_test_suite)

UT_TEST_CASE(test1) {
  int ws_port = 7681;

  // Sphere mesh(4);
  // mesh.vertices().print();
  // mesh.triangles().print();
  Grid<Polygon> mesh({10, 10}, 3);
  mesh.fields().set_defaults(mesh);
  LOG << "# triangles = " << mesh.triangles().n();
  // Viewer viewer(mesh, ws_port);
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(graphics_test_suite)