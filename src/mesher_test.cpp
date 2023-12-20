#include "mesher.h"

#include "halfedges.h"
#include "io.h"
#include "mesh.h"
#include "numerics.h"
#include "tester.h"
#include "texture.h"

using namespace vortex;

UT_TEST_SUITE(mesher_test_suite)

UT_TEST_CASE(test1) {
  TextureOptions tex_opts;
  tex_opts.format = TextureFormat::kGrayscale;
  std::string filename = "oceans_2048.png";
  Texture texture(filename, tex_opts);
  texture.make_binary(10, 10, 255);
  texture.smooth(10);
  texture.make_periodic();
  texture.write("texture.jpg");

  MeshingParameters msh_opts;
  msh_opts.max_iter = 10;
  msh_opts.h_min = 0.05;
  msh_opts.h_max = 0.1;
  EarthMesher mesher(texture);
  mesher.generate(msh_opts);

  Texture mask("oceans_2048.png", tex_opts);
  mask.make_binary(10, 10, 255);

  auto& mesh = mesher.mesh();
  mesh.deactivate_by([&mask](const auto& n) {
    vec3d p(n.point());
    vec3d u;
    sphere_params(p, u);
    double t;
    mask.sample(u[0], u[1], &t);
    if (t < 11) return true;
    return false;
  });

  Mesh output_mesh(3);
  mesh.extract(output_mesh);

  meshb::write(output_mesh, "earth3.meshb");
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(mesher_test_suite)