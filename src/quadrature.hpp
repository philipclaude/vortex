#include "quadrature.h"

namespace vortex {

template <typename T>
template <typename Integrand>
coord_t TriangleQuadrature<T>::integrate(const Integrand& fn, const coord_t* pa,
                                         const coord_t* pb,
                                         const coord_t* pc) const {
  coord_t integral = 0.0;
  for (size_t k = 0; k < points_.size(); k++) {
    // evaluate the point in physical space
    vec3d point = T::get_physical_coordinates(pa, pb, pc, &points_[k][0]);

    // add the contribution to the integral
    integral += weights_[k] * fn(point);
  }
  // technically the jacobian determinant should be used at each quadrature
  // point but maybe this will work fine
  return integral * 2 * T::area(pa, pb, pc);
}

}  // namespace vortex