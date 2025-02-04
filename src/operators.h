#pragma once

#include <vector>

#include "defs.h"

namespace vortex {

class VoronoiDiagram;

template <typename Domain_t>
class VoronoiOperators {
 public:
  VoronoiOperators(const VoronoiDiagram& voronoi);

  void set_boundary_value(double x) { boundary_value_ = x; }
  void calculate_gradient(const coord_t* f, coord_t* grad_f);
  void calculate_divergence(const coord_t* u, coord_t* div_u);

 private:
  const VoronoiDiagram& voronoi_;
  double boundary_value_{1e20};
};

}  // namespace vortex