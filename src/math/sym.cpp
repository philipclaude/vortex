#include "sym.h"

#include <math.h>
#include <stdio.h>

#include "linalg.h"
#include "mat.h"
#include "mat.hpp"
#include "sym.hpp"

namespace vortex {

template <typename type>
symd<type>::symd(const int _n) : n_(_n) {
  // constructor from size
  data_.resize(nb());
  std::fill(data_.begin(), data_.end(), 0.);
}

template <typename type>
symd<type>::symd(const vecd<type>& lambda, const matd<type>& q)
    : n_(lambda.m()) {
  from_eig(lambda, q);
}

template <typename type>
symd<type>::symd(const std::pair<vecd<type>, matd<type> >& decomp)
    : symd(decomp.first, decomp.second) {}

template <typename type>
symd<type>::symd(const matd<type>& M) : n_(M.n()) {
  ASSERT(M.m() == M.n());
  data_.resize(nb());
  for (int i = 0; i < n_; i++)
    for (int j = 0; j < n_; j++) operator()(i, j) = M(i, j);
}

template <typename type>
void symd<type>::set(const matd<type>& M) {
  n_ = M.m();
  ASSERT(n_ == M.n());
  for (int i = 0; i < n_; i++)
    for (int j = 0; j < n_; j++) (*this)(i, j) = M(i, j);
}

template <typename type>
void symd<type>::from_eig(const vecd<type>& lambda, const matd<type>& q) {
  // constructor from eigendecomposition
  n_ = lambda.m();
  data_.resize(nb());
  matd<type> M = q * (diag(lambda) * transpose(q));
  set(M);
}

// eigenvalues and eigenvectors
template <typename type>
std::pair<vecd<type>, matd<type> > symd<type>::eig() const {
  std::pair<vecd<type>, matd<type> > decomp = __eig__();
  return decomp;
}

template <typename type>
std::pair<vecd<type>, matd<type> > symd<type>::__eig__() const {
  // please see the following paper:
  // "Efficient numerical diagonalization of hermitian 3x3 matrices"
  // by Joachim Kopp
  // Int. J. Mod. Phys. C 19 (2008) 523-548
  // arXiv.org: physics/0610206
  // (https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/)
  vecd<type> L(n_);
  matd<type> E(n_, n_);

  type sd, so;
  type s, c, t;
  type g, h, z, theta;
  type thresh;

  E.eye();
  symd A(n_);
  A.copy(*this);

  for (int i = 0; i < n_; i++) L[i] = A(i, i);

  // calculate the square of the trace
  sd = 0.;
  for (int i = 0; i < n_; i++) sd += fabs(L[i]);
  sd = sd * sd;

  for (int iter = 0; iter < 50; iter++) {
    // test for convergence
    so = 0.;
    for (int p = 0; p < n_; p++)
      for (int q = p + 1; q < n_; q++) so += fabs(A(p, q));

    if (so == 0.0) return std::make_pair(L, E);

    if (iter < 4)
      thresh = 0.2 * so / type(n_ * n_);
    else
      thresh = 0.;

    // sweep
    for (int p = 0; p < n_; p++) {
      for (int q = p + 1; q < n_; q++) {
        g = 100. * fabs(A(p, q));
        if (iter > 4 && fabs(L[p]) + g == fabs(L[p]) &&
            fabs(L[q]) + g == fabs(L[q]))
          A(p, q) = 0.;
        else if (fabs(A(p, q)) > thresh) {
          // calculate Jacobi transformation
          h = L[q] - L[p];
          if (fabs(h) + g == fabs(h)) {
            t = A(p, q) / h;
          } else {
            theta = 0.5 * h / A(p, q);
            if (theta < 0.0)
              t = -1. / (::sqrt(1. + theta * theta) - theta);
            else
              t = 1. / (::sqrt(1. + theta * theta) + theta);
          }
          c = 1. / ::sqrt(1. + t * t);
          s = t * c;
          z = t * A(p, q);

          // apply Jacobi transformation
          A(p, q) = 0.;
          L[p] -= z;
          L[q] += z;
          for (int r = 0; r < p; r++) {
            t = A(r, p);
            A(r, p) = c * t - s * A(r, q);
            A(r, q) = s * t + c * A(r, q);
          }

          for (int r = p + 1; r < q; r++) {
            t = A(p, r);
            A(p, r) = c * t - s * A(r, q);
            A(r, q) = s * t + c * A(r, q);
          }

          for (int r = q + 1; r < n_; r++) {
            t = A(p, r);
            A(p, r) = c * t - s * A(q, r);
            A(q, r) = s * t + c * A(q, r);
          }

          // update eigenvectors
          for (int r = 0; r < n_; r++) {
            t = E(r, p);
            E(r, p) = c * t - s * E(r, q);
            E(r, q) = s * t + c * E(r, q);
          }
        }
      }
    }
  }
  print();
  NOT_POSSIBLE;
  return std::make_pair(L, E);
}

template <typename T>
symd<T> operator+(const symd<T>& A) {
  return A;
}

template <typename T>
symd<T> operator-(const symd<T>& A) {
  symd<T> B(A.m());
  for (int i = 0; i < A.m(); i++)
    for (int j = 0; j < A.m(); j++) B(i, j) = -A(i, j);
  return B;
}

template <typename T>
vecd<T> operator*(const symd<T>& A, const vecd<T>& x) {
  LLAMA_ASSERT_msg(A.m() == x.m(), "bad sizes");
  vecd<T> b(A.m());
  for (int i = 0; i < A.m(); i++) {
    T sum = 0;
    for (int k = 0; k < A.m(); k++) sum += A(i, k) * x(k);
    b(i) = sum;
  }
  return b;
}

template <typename type>
void symd<type>::print(const std::string& title) const {
  if (!title.empty()) std::cout << title << std::endl;
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  for (int i = 0; i < n_; i++)
    for (int j = 0; j < n_; j++)
      std::cout << "(" + std::to_string(i) + "," + std::to_string(j) + "): "
                << (*this)(i, j) << std::endl;
}

// ======== syms ======

template <int N, typename T>
void syms<N, T>::from_eig(const vecs<N, T>& lambda, const mats<N, N, T>& q) {
  mats<N, N, T> M = q * (diag(lambda) * transpose(q));
  set(M);
}

// eigenvalues and eigenvectors
template <int N, typename T>
std::pair<vecs<N, T>, mats<N, N, T> > syms<N, T>::eig() const {
  std::pair<vecs<N, T>, mats<N, N, T> > decomp = __eig__();
  return decomp;
}

template <int N, typename T>
std::pair<vecs<N, T>, mats<N, N, T> > syms<N, T>::__eig__() const {
  // please see the following paper:
  // "Efficient numerical diagonalization of hermitian 3x3 matrices"
  // by Joachim Kopp
  // Int. J. Mod. Phys. C 19 (2008) 523-548
  // arXiv.org: physics/0610206
  // (https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/)
  vecs<N, T> L;
  mats<N, N, T> E;

  T sd, so;
  T s, c, t;
  T g, h, z, theta;
  T thresh;

  E.eye();
  syms<N, T> A;
  A.copy(*this);

  for (int i = 0; i < N; i++) L[i] = A(i, i);

  // calculate the square of the trace
  sd = 0.;
  for (int i = 0; i < N; i++) sd += fabs(L[i]);
  sd = sd * sd;

  for (int iter = 0; iter < 50; iter++) {
    // test for convergence
    so = 0.;
    for (int p = 0; p < N; p++)
      for (int q = p + 1; q < N; q++) so += fabs(A(p, q));

    if (so == 0.0) return std::make_pair(L, E);

    if (iter < 4)
      thresh = 0.2 * so / T(N * N);
    else
      thresh = 0.;

    // sweep
    for (int p = 0; p < N; p++) {
      for (int q = p + 1; q < N; q++) {
        g = 100. * fabs(A(p, q));
        if (iter > 4 && fabs(L[p]) + g == fabs(L[p]) &&
            fabs(L[q]) + g == fabs(L[q]))
          A(p, q) = 0.;
        else if (fabs(A(p, q)) > thresh) {
          // calculate Jacobi transformation
          h = L[q] - L[p];
          if (fabs(h) + g == fabs(h)) {
            t = A(p, q) / h;
          } else {
            theta = 0.5 * h / A(p, q);
            if (theta < 0.0)
              t = -1. / (::sqrt(1. + theta * theta) - theta);
            else
              t = 1. / (::sqrt(1. + theta * theta) + theta);
          }
          c = 1. / ::sqrt(1. + t * t);
          s = t * c;
          z = t * A(p, q);

          // apply Jacobi transformation
          A(p, q) = 0.;
          L[p] -= z;
          L[q] += z;
          for (int r = 0; r < p; r++) {
            t = A(r, p);
            A(r, p) = c * t - s * A(r, q);
            A(r, q) = s * t + c * A(r, q);
          }

          for (int r = p + 1; r < q; r++) {
            t = A(p, r);
            A(p, r) = c * t - s * A(r, q);
            A(r, q) = s * t + c * A(r, q);
          }

          for (int r = q + 1; r < N; r++) {
            t = A(p, r);
            A(p, r) = c * t - s * A(q, r);
            A(q, r) = s * t + c * A(q, r);
          }

          // update eigenvectors
          for (int r = 0; r < N; r++) {
            t = E(r, p);
            E(r, p) = c * t - s * E(r, q);
            E(r, q) = s * t + c * E(r, q);
          }
        }
      }
    }
  }
  print();
  NOT_POSSIBLE;
  return std::make_pair(L, E);
}

template <int N, typename T>
syms<N, T> operator-(const syms<N, T>& A) {
  syms<N, T> B;
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) B(i, j) = -A(i, j);
  return B;
}

template <int N, typename T>
vecs<N, T> operator*(const syms<N, T>& A, const vecs<N, T>& x) {
  vecs<N, T> b;
  for (int i = 0; i < N; i++) {
    T sum = 0;
    for (int k = 0; k < N; k++) sum += A(i, k) * x(k);
    b(i) = sum;
  }
  return b;
}

template <int N, typename T>
void syms<N, T>::print(const std::string& title) const {
  if (!title.empty()) std::cout << title << std::endl;
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      std::cout << "(" + std::to_string(i) + "," + std::to_string(j) + "): "
                << (*this)(i, j) << std::endl;
}

template class symd<double>;
template class symd<float>;

template class syms<2, double>;
template class syms<3, double>;

}  // namespace vortex
