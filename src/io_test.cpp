#include "io.h"

#include "library.h"
#include "mesh.h"
#include "numerics.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(IO_TESTS)

UT_TEST_CASE(test1) {
  Mesh mesh(3);
  obj::read("/Users/philip/Desktop/earth-r50.obj", mesh);
  for (size_t i = 0; i < mesh.vertices().n(); i++) {
    vec3d p(mesh.vertices()[i]);
    p = normalize(p);
    for (int d = 0; d < 3; d++) mesh.vertices()[i][d] = p[d];
  }
  obj::write(mesh, "../data/earth.obj");
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(IO_TESTS)