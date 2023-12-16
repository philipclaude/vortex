#pragma once

namespace vortex {

class Texture;
class Mesh;

struct MeshingParameters {
  double h_min{0.025};
  double h_max{0.05};
  int max_iter{5};
};

class EarthMesher {
 public:
  EarthMesher(const Texture& texture);
  void generate(MeshingParameters params, Mesh& output_mesh) const;

 private:
  const Texture& texture_;
};

}  // namespace vortex