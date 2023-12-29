#include "io.h"

#include "library.h"
#include "mesh.h"
#include "numerics.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(io_tests)

UT_TEST_CASE(obj_test) {
  Sphere sphere(4);
  obj::write(sphere, "sphere.obj");

  Mesh mesh(3);
  read_mesh("sphere.obj", mesh);
  UT_ASSERT_EQUALS(mesh.vertices().dim(), 3);
  UT_ASSERT_EQUALS(mesh.vertices().n(), sphere.vertices().n());
  UT_ASSERT_EQUALS(mesh.triangles().n(), sphere.triangles().n());
}
UT_TEST_CASE_END(obj_test)

UT_TEST_CASE(meshb_test) {
  Sphere sphere(4);
  meshb::write(sphere, "sphere.meshb");

  Mesh mesh(3);
  read_mesh("sphere.meshb", mesh);
  UT_ASSERT_EQUALS(mesh.vertices().dim(), 3);
  UT_ASSERT_EQUALS(mesh.vertices().n(), sphere.vertices().n());
  UT_ASSERT_EQUALS(mesh.triangles().n(), sphere.triangles().n());
}
UT_TEST_CASE_END(meshb_test)

UT_TEST_SUITE_END(io_tests)