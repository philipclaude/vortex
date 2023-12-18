#pragma once

#include <memory>

namespace vortex {

class Texture;
class Mesh;
class HalfMesh;

struct MeshingParameters {
  double h_min{0.025};
  double h_max{0.05};
  int max_iter{5};
};

class EarthMesher {
 public:
  EarthMesher(const Texture& texture);
  ~EarthMesher();
  void generate(MeshingParameters params);
  auto& mesh() { return *mesh_; }
  void extract(Mesh& mesh) const;

 private:
  std::unique_ptr<HalfMesh> mesh_;
  const Texture& texture_;
};

}  // namespace vortex