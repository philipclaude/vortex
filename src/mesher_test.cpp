#include "mesher.h"

#include "io.h"
#include "mesh.h"
#include "tester.h"
#include "texture.h"

using namespace vortex;

UT_TEST_SUITE(mesher_test_suite)

UT_TEST_CASE(test1) {
  TextureOptions tex_opts{.format = TextureFormat::kGrayscale};
  std::string filename = "lcc_global_2048.png";
  Texture texture(filename, tex_opts);
  texture.make_binary(10, 10, 255);
  texture.smooth(10);
  texture.make_periodic();
  texture.write("texture.jpg");

  Mesh mesh(3);
  MeshingParameters msh_opts{.max_iter = 10, .h_min = 0.005, .h_max = 0.01};
  EarthMesher mesher(texture);
  mesher.generate(msh_opts, mesh);

  meshb::write(mesh, "earth.meshb");
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(mesher_test_suite)