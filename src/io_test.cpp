#include "io.h"

#include "mesh.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(IO_TESTS)

UT_TEST_CASE(shp_test) {
  Mesh mesh(3);
  int res = 110;
  std::string name = "coastline";
  std::string base = fmt::format(
      "/Users/philip/Codes/external/natural-earth-vector/{}m_physical/"
      "ne_{}m_{}",
      res, res, name);
  shp::read(base, mesh);
  meshb::write(mesh, "../data/test.meshb");
}
UT_TEST_CASE_END(shp_test)

UT_TEST_SUITE_END(IO_TESTS)