#pragma once

namespace vortex {

class Mesh;

class OceanTriangulator {
 public:
  struct Options {};
  OceanTriangulator(Mesh& mesh, const Mesh& coast);
  void triangulate();

 private:
  void insert_points();
  void recover_edges();

  Mesh& mesh_;
  const Mesh& coast_;
};

}  // namespace vortex