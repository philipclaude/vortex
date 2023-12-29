#include "vec.hpp"

namespace vortex {

#define INSTANTIATE_VECD(R, S)                                             \
  template vecd<typename result_of<R, S>::type> operator+(const vecd<R>&,  \
                                                          const vecd<S>&); \
  template vecd<typename result_of<R, S>::type> operator-(const vecd<R>&,  \
                                                          const vecd<S>&); \
  template vecd<typename result_of<R, S>::type> operator*(const R&,        \
                                                          const vecd<S>&); \
  template vecd<typename result_of<R, S>::type> operator*(const vecd<R>&,  \
                                                          const S&);
INSTANTIATE_VECD(double, double)

#undef INSTANTIATE_VECD

template <typename T>
vecs<3, T> cross(const vecs<3, T>& u, const vecs<3, T>& v) {
  vecs<3, T> w;
  w(0) = u(1) * v(2) - u(2) * v(1);
  w(1) = -(u(0) * v(2) - u(2) * v(0));
  w(2) = u(0) * v(1) - u(1) * v(0);
  return w;
}

// template<typename R,typename S,int M>
// typename result_of<R,S>::type
// dot( const vecs<M,R>& u , const vecs<M,S>& v ) {
//   typename result_of<R,S>::type d(0);
//   for (int i = 0; i < M; i++)
//     d += u[i]*v[i];
//   return d;
// }

template <typename R, typename S>
typename result_of<R, S>::type dot(const vecd<R>& u, const vecd<S>& v) {
  ASSERT(u.m() == v.m());
  typename result_of<R, S>::type d(0);
  for (int i = 0; i < u.m(); i++) d += u[i] * v[i];
  return d;
}

template <typename T>
T length(const vecd<T>& u) {
  T m(0);
  for (int i = 0; i < u.m(); i++) m += u[i] * u[i];
  return std::sqrt(m);
}

template vecs<3, double> cross(const vecs<3, double>&, const vecs<3, double>&);
template vecs<3, float> cross(const vecs<3, float>&, const vecs<3, float>&);

// template typename result_of<double,double>::type dot<double,double,3>( const
// vecs<3,double>& , const vecs<3,double>& );
template typename result_of<double, double>::type dot<double, double>(
    const vecd<double>&, const vecd<double>&);

template double length(const vecd<double>&);
template float length(const vecd<float>&);

}  // namespace vortex
