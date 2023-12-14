#include "mesher.h"

#include <algorithm>
#include <cmath>

#include "halfedges.h"
#include "library.h"

namespace vortex {

namespace {

void collapse(HalfMesh& mesh, const Texture& texture) {
  //
}

void split(HalfMesh& mesh, const Texture& texture) {
  //
}

void flips(HalfMesh& mesh, const Texture& texture) {
  //
}

void smooth(HalfMesh& mesh, const Texture& texture) {
  //
}

}  // namespace

EarthMesher::EarthMesher(const Texture& texture) : texture_(texture) {}

void EarthMesher::generate(MeshingParameters params, Mesh& output_mesh) const {
  // estimate the number of levels from the mean of the min and max sizes
  double h_avg = 0.5 * (params.h_min + params.h_max);
  double at = std::sqrt(3.0) * h_avg * h_avg / 4.0;
  int n_triangles = std::floor(4.0 * M_PI / at);
  int n_level = int(std::log(n_triangles / 20) / std::log(4.0));
  LOG << fmt::format("n_level = {}", n_level);
  Sphere sphere(n_level);

  HalfMesh mesh(sphere);

  for (int iter = 0; iter < params.max_iter; iter++) {
    // collapse short edges
    collapse(mesh, texture_);

    // split long edges without creating short edges
    split(mesh, texture_);

    // optimize by flips
    flips(mesh, texture_);

    // optimize by smoothing
    smooth(mesh, texture_);
  }

  mesh.extract(output_mesh);
}

}  // namespace vortex