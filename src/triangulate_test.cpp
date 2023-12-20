#include "triangulate.h"

#include "io.h"
#include "library.h"
#include "mesh.h"
#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(TRIANGULATE_TESTS)

UT_TEST_CASE(test1) {
  Mesh input_mesh(3);
  read_mesh("../build/release/test.meshb", input_mesh);

  // only keep vertices on the lines
  Mesh coast(3);
  std::unordered_map<index_t, index_t> vertex_map;
  vertex_map.reserve(input_mesh.vertices().n());
  for (size_t k = 0; k < input_mesh.lines().n(); k++) {
    auto* e = input_mesh.lines()[k];
    for (int j = 0; j < 2; j++) {
      if (vertex_map.find(e[j]) == vertex_map.end()) {
        vertex_map.insert({e[j], vertex_map.size()});
        coast.vertices().add(input_mesh.vertices()[e[j]]);
      }
      e[j] = vertex_map.at(e[j]);
    }
  }

  Sphere oceans(4);
  OceanTriangulator triangulator(oceans, coast);
  triangulator.triangulate();

  meshb::write(oceans, "../data/oceans.meshb");
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(TRIANGULATE_TESTS)