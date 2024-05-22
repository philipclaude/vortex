#include "spmat.h"

#include <OpenNL_psm.h>

#include "defs.h"
#include "vec.h"

namespace vortex {

template <>
void spmat<double>::solve_nl(const vecd<double>& b, vecd<double>& x, double tol,
                             bool symmetric) const {
  ASSERT(b.m() == n_rows());

  nlNewContext();
  nlSolverParameteri(NL_NB_VARIABLES, NLint(n_rows()));
  nlSolverParameteri(NL_NB_SYSTEMS, 1);
  if (symmetric) {
    nlSolverParameteri(NL_SOLVER, NL_CG);
    nlSolverParameteri(NL_SYMMETRIC, NL_TRUE);
  } else {
    nlSolverParameteri(NL_SOLVER, NL_BICGSTAB);
    nlSolverParameteri(NL_SYMMETRIC, NL_FALSE);
  }
  nlSolverParameteri(NL_PRECONDITIONER, NL_PRECOND_JACOBI);
  nlSolverParameterd(NL_THRESHOLD, tol);
  nlSolverParameteri(NL_MAX_ITERATIONS, 100);

  nlBegin(NL_SYSTEM);
  nlBegin(NL_MATRIX);

  for (size_t k = 0; k < rows_.size(); k++) {
    const auto& r = rows_[k];
    nlBegin(NL_ROW);
    for (const auto& [col, value] : r) {
      nlCoefficient(col, value);
      nlRightHandSide(b[k]);
    }
    nlEnd(NL_ROW);
  }

  nlEnd(NL_MATRIX);
  nlEnd(NL_SYSTEM);

  nlSolve();

  ASSERT(x.m() == b.m());
  for (auto k = 0; k < n_rows(); k++) x[k] = nlGetVariable(k);

  nlDeleteContext(nlGetCurrent());
}

template <typename T>
double spmat<T>::solve_jacobi(const vecd<T>& b, vecd<T>& x, double tol,
                              int max_iter, bool verbose) const {
  int n = b.m();
  if (max_iter < 0) max_iter = n * n;

  // dereference the matrix for readability
  const spmat<T>& A = *this;

  vecd<T> r(n);
  int iter = 0;
  double e = norm(A * x - b);
  while (e > tol && iter++ < max_iter) {
    r = b - A * x;
    for (int i = 0; i < n; i++)
      x(i) =
          (r(i) + A(i, i) * x(i)) / A(i, i);  // re-add the diagonal component

    e = norm(r);
    if (verbose) LOG << fmt::format("iter {}, e = {}", iter, e);
  }

  return e;
}

template class spmat<double>;

}  // namespace vortex
