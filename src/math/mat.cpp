#include "mat.hpp"

#include "sym.h"
#include "vec.h"

namespace vortex {

template <typename T>
matd<T> operator+(const matd<T>& A) {
  return A;
}

template <typename T>
matd<T> operator-(const matd<T>& A) {
  matd<T> B(A.m(), A.n());
  for (int i = 0; i < A.m(); i++)
    for (int j = 0; j < A.n(); j++) B(i, j) = -A(i, j);
  return B;
}

template <typename T>
vecd<T> operator*(const matd<T>& A, const vecd<T>& x) {
  ASSERT(A.n() == x.m()) << "bad sizes";
  vecd<T> b(A.m());
  for (int i = 0; i < A.m(); i++) {
    T sum = 0;
    for (int k = 0; k < A.n(); k++) sum += A(i, k) * x(k);
    b(i) = sum;
  }
  return b;
}

template vecd<double> operator*(const matd<double>&, const vecd<double>&);

#define INSTANTIATE_MATD(T) template class matd<T>;
INSTANTIATE_MATD(double)
#undef INSTANTIATE_MATD

/*
#define INSTANTIATE_MUL(R,S) template symd< typename result_of<R,S>::type >
operator*( const symd<R>& , const symd<S>& ); INSTANTIATE_MUL( double , double )
#undef INSTANTIATE_MUL
*/

#define INSTANTIATE_MULD(R, S)                                            \
  template matd<typename result_of<R, S>::type> operator*(const matd<R>&, \
                                                          const matd<S>&);
INSTANTIATE_MULD(double, double)
INSTANTIATE_MULD(float, float)
#undef INSTANTIATE_MULD

#define INSTANTIATE_OPD(R, S)                                              \
  template matd<typename result_of<R, S>::type> operator+(const matd<R>&,  \
                                                          const matd<S>&); \
  template matd<typename result_of<R, S>::type> operator-(const matd<R>&,  \
                                                          const matd<S>&);
INSTANTIATE_OPD(double, double)
#undef INSTANTIATE_OPD

#define INSTANTIATE_ASSIGN(T, S) \
  template matd<T>& matd<T>::operator=(const symd<S>& M);

// INSTANTIATE_ASSIGN(double,double)
#undef INSTANTIATE_ASSIGN

template matd<double> operator-(const matd<double>&);
template matd<double> operator+(const matd<double>&);

}  // namespace vortex
