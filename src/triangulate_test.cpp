#include "triangulate.h"

#include "io.h"
#include "library.h"
#include "mesh.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(TRIANGULATE_TESTS)

UT_TEST_CASE(test1) {
  Mesh mesh(3);
  int res = 110;
  std::string name = "coastline";
  std::string base = fmt::format(
      "/Users/philip/Codes/external/natural-earth-vector/{}m_physical/"
      "ne_{}m_{}",
      res, res, name);
  shp::read(base, mesh);

  Sphere oceans(4);
  OceanTriangulator triangulator(oceans, mesh);
  triangulator.triangulate();

  meshb::write(oceans, "../data/oceans.meshb");
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(TRIANGULATE_TESTS)