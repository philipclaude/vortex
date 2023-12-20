#include "numerics.h"

#include "predicates.h"

#if ROBUST_PREDICATES == 0

void exactinit() {}
double orient2d(double* a, double* b, double* c) {
  const double acx = a[0] - c[0];
  const double bcx = b[0] - c[0];
  const double acy = a[1] - c[1];
  const double bcy = b[1] - c[1];
  return acx * bcy - acy * bcx;
}

#endif