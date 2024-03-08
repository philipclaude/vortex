#pragma once

#include <vector>

#include "defs.h"
#include "elements.h"
#include "math/vec.h"

namespace vortex {

template <typename T>
class TriangleQuadrature {
 public:
  static const int dim = T::dimension;
  using point_t = vecs<dim, coord_t>;

  TriangleQuadrature(int order) { define(order); }

  template <typename Integrand>
  coord_t integrate(const Integrand& fn, const coord_t* pa, const coord_t* pb,
                    const coord_t* pc) const;

 private:
  void define(int order);
  std::vector<point_t> points_;
  std::vector<coord_t> weights_;
};

}  // namespace vortex