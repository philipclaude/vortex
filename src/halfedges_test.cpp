#include "halfedges.h"

#include "io.h"
#include "library.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(halfedges_test_suite)

UT_TEST_CASE(test_closed) {
  Sphere mesh(3);
  HalfMesh hmesh(mesh);
  UT_ASSERT(hmesh.check());
}
UT_TEST_CASE_END(test_closed)

UT_TEST_CASE(test_open) {
  Grid<Triangle> mesh({10, 10}, 3);
  HalfMesh hmesh(mesh);
  UT_ASSERT(hmesh.check());

  Mesh m(3);
  hmesh.extract(m);
  meshb::write(m, "grid.meshb");
}
UT_TEST_CASE_END(test_open)

UT_TEST_SUITE_END(halfedges_test_suite)
