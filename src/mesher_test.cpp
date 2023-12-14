#include "mesher.h"

#include "io.h"
#include "mesh.h"
#include "tester.h"
#include "texture.h"

using namespace vortex;

UT_TEST_SUITE(mesher_test_suite)

UT_TEST_CASE(test1) {
  TextureOptions tex_opts{.format = TextureFormat::kGrayscale};
  std::string filename = "land-ocean-ice.png";
  Texture texture(filename, tex_opts);

  Mesh mesh(3);
  MeshingParameters msh_opts;
  EarthMesher mesher(texture);
  mesher.generate(msh_opts, mesh);

  meshb::write(mesh, "earth.meshb");
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(mesher_test_suite)