#include "mesh.h"

#include "io.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(mesh_test_suite)

UT_TEST_CASE(test1) {
  Mesh mesh(3);
  index_t triangle[3] = {1, 2, 3};
  mesh.triangles().add(triangle);
  mesh.triangles().print();
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(mesh_test_suite)