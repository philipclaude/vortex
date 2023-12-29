#pragma once

#include <string>
#include <utility>
#include <vector>

#include "log.h"
#include "mat.h"
#include "vec.h"

namespace vortex {

/**
 * \brief Represents a symmetric matrix, thus reducing storage by only saving
 * the unique matrix entries. This is a square matrix, so we only need to store
 * a single value for the number of rows or columns (n). Only dynamic storage is
 * currently supported (unlike the difference between matd and mats).
 */
template <typename type>
class symd {
 public:
  /**
   * \brief Constructs a 0x0 matrix.
   */
  symd() : n_(0) {}

  /**
   * \brief Constructs a n x n matrix.
   */
  symd(const int _n);

  /**
   * \brief Constructs this matrix from the eigendecomposition Q * D * Q^T
   *
   * \param[in] lambda - vector of n eigenvalues
   * \param[in] q -  n x n matrix of eigenvectors
   */
  symd(const vecd<type>& lambda, const matd<type>& q);

  /**
   * \brief Constructs this matrix from the eigendecomposition Q * D * Q^T
   *
   * \param[in] decomp (first = vector of n eigenvalues, second = n x n matrix
   * of eigenvectors)
   */
  symd(const std::pair<vecd<type>, matd<type> >& decomp);

  /**
   * \brief Constructs a symmetric matrix from a general matrix (which is
   * hopefully symmetric).
   *
   * \param[in] A - n x n matrix which is hopefully symmetric
   */
  symd(const matd<type>& A);

  /**
   * \brief Constructs a matrix by copying values from some other symmetric
   * matrix.
   *
   * \param[in] A - n x n symmetric matrix with a potentially different type S
   */
  template <typename S>
  symd(const symd<S>& A) : symd(A.n()) {
    for (int i = 0; i < n_; i++)
      for (int j = 0; j < n_; j++) (*this)(i, j) = A(i, j);
  }

  /**
   * \brief Constructs a matrix by copying values from some other matrix
   * (hopefully symemtric)
   *
   * \param[in] A - n x n matrix (hopefully symmetric) with a potentially
   * different type S
   */
  template <typename S>
  symd(const matd<S>& A) : symd(A.n()) {
    ASSERT(A.m() == A.n());
    for (int i = 0; i < n_; i++)
      for (int j = 0; j < n_; j++) (*this)(i, j) = A(i, j);
  }

  /**
   * \brief Reallocates storage for this symmetric matrix.
   *
   * \param[in] n - the number of rows and columns of this matrix.
   */
  void resize(const int n) {
    n_ = n;
    data_.resize(nb());
  }

  /**
   * \brief Constructs this matrix from the eigendecomposition Q * D * Q^T
   *
   * \param[in] lambda - vector of n eigenvalues
   * \param[in] q -  n x n matrix of eigenvectors
   */
  void from_eig(const vecd<type>& lambda, const matd<type>& q);

  /**
   * \brief Returns the number of unique entries stored in this matrix.
   *
   * \return number of unique entries
   */
  int nb() const { return n_ * (n_ + 1) / 2; }

  /**
   * \brief Returns the number of columns.
   */
  int n() const { return n_; }

  /**
   * \brief Returns the number of rows.
   */
  int m() const { return n_; }

  /**
   * \brief Copies the entries from some other symmetric matrix into this
   * matrix.
   *
   * \param[in] T - symmetric matrix to copy
   */
  void copy(const symd& T) {
    ASSERT(n_ == T.n());
    for (int k = 0; k < nb(); k++) data_[k] = T.data(k);
  }

  /**
   * \brief Copies the entries from some other matrix into this matrix.
   *
   * \param[in] T - matrix to copy (hopefully symmetric)
   */
  void set(const matd<type>& M);

  /**
   * \brief Copies the entries from some other symmetric matrix into this
   * matrix.
   *
   * \param[in] T - symmetric matrix to copy
   */
  void set(const symd<type>& S) {
    for (int i = 0; i < nb(); i++) data_[i] = S.data(i);
  }

  /**
   * \brief Sets all entries in this symmetric matrix to zero.
   */
  void zero() {
    for (int i = 0; i < nb(); i++) data_[i] = 0.;
  }

  /**
   * \brief Sets this matrix to the identity matrix (diagonal entries = 1,
   * off-diagonals = 0)
   */
  void eye() {
    zero();
    for (int i = 0; i < n_; i++) operator()(i, i) = 1.;
  }

  /**
   * \brief Read/write access to the entry into the flattened 1d array of
   * entries. This function should not be used often.
   *
   * \param[in] k - index of the item (in the flattened array) to retrieve.
   *
   * \return entry at index k
   */
  type& data(const int k) { return data_[k]; }

  /**
   * \brief Read access to the entry into the flattened 1d array of entries.
   *        This function should not be used often.
   *
   * \param[in] k - index of the item (in the flattened array) to retrieve.
   *
   * \return entry at index k
   */
  type data(const int k) const { return data_[k]; }

  /**
   * \brief Read/write access to entry (i,j)
   *
   * \param[in] i - row index
   * \param[in] j - column index
   *
   * \return reference to entry (i,j)
   */
  inline type& operator()(const int i, const int j) {
    ASSERT(i < n_ && j < n_);
    return (i > j) ? data_[i * (i + 1) / 2 + j] : data_[j * (j + 1) / 2 + i];
  }

  /**
   * \brief Read access to entry (i,j)
   *
   * \param[in] i - row index
   * \param[in] j - column index
   *
   * \return reference to entry (i,j)
   */
  inline type operator()(const int i, const int j) const {
    ASSERT(i < n_ && j < n_);
    return (i > j) ? data_[i * (i + 1) / 2 + j] : data_[j * (j + 1) / 2 + i];
  }

  /**
   * \brief Unary operator, so that symmetric matrices can be included in
   * expressions such as +A
   */
  symd operator+() const { return *this; }

  /**
   * \brief Matrix-matrix addition: (*this) + U
   *
   * \param[in] this - left operand in addition
   * \param[in] U - right operand in addition
   *
   * \return new symmetric matrix (*this) + U
   */
  symd operator+(const symd& U) {
    symd V(U.n());
    for (int i = 0; i < U.nb(); i++) V.data(i) = data_[i] + U.data(i);
    return V;
  }

  /**
   * \brief Matrix-scalar multiplication.
   *
   * \param[in] this - matrix to multiply
   * \param[in] a - scalar
   *
   * \return (*this) * a
   */
  symd operator*(const type a) const {
    symd C(n());
    for (int i = 0; i < nb(); i++) C.data(i) = a * data_[i];
    return C;
  }

  /**
   * \brief Assignment operator from a matrix - saves entries of A to this
   * matrix.
   *
   * \param[in] A - n x n (hopefully symmetric) matrix
   *
   * \return (*this)
   */
  symd& operator=(const matd<type>& A) {
    // assumes A is already symmetric
    ASSERT(A.m() == A.n());
    resize(A.m());
    for (int i = 0; i < n_; i++)
      for (int j = 0; j < n_; j++) (*this)(i, j) = A(i, j);
    return *this;
  }

  /**
   * \brief Assignment operator from a scalar - sets all entries to a constant.
   *
   * \param a - scalar all entries will be set to
   *
   * \return (*this) where all entries are equal to a
   */
  symd& operator=(const double& a) {
    std::fill(data_.begin(), data_.end(), a);
    return *this;
  }

  /**
   * \brief Computes the eigendecomposition of a symmetric matrix,
   *        i.e. M = Q * D * Q^T, returns Q * exp(D) * Q^T,
   *        where D is a diagonal matrix with the eigenvalues
   *        and Q are the eigenvectors.
   *
   * \param[in] M - the symmetric d x d matrix
   *
   * \return pair {D,Q} where D is a d-vector of diagonal entries and Q  is a d
   * x d matrix of eigenvectors
   */
  std::pair<vecd<type>, matd<type> > eig() const;

  /**
   * \brief Prints this matrix to the console.
   *
   * \param[in] title (optional) prefix for printing entries
   *
   */
  void print(const std::string& title = std::string()) const;

 private:
  int n_;                   // number of rows and also number of columns
  std::vector<type> data_;  // flattened 1d array used to store matrix entries

  /**
   * \brief Computes the eigendecomposition of a symmetric matrix using the
   * method described here: https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/
   *        i.e. M = Q * D * Q^T, returns Q * exp(D) * Q^T,
   *        where D is a diagonal matrix with the eigenvalues
   *        and Q are the eigenvectors.
   *
   * \param[in] M - the symmetric d x d matrix
   *
   * \return pair {D,Q} where D is a d-vector of diagonal entries and Q  is a d
   * x d matrix of eigenvectors
   */
  std::pair<vecd<type>, matd<type> > __eig__() const;
};

template <int N, typename T>
class syms {
 private:
  static const int NB = N * (N + 1) / 2;

 public:
  /**
   * \brief Constructs a NxN symmetric matrix.
   */
  syms() { zero(); }

  /**
   * \brief Constructs this matrix from the eigendecomposition Q * D * Q^T
   *
   * \param[in] lambda - vector of n eigenvalues
   * \param[in] q -  n x n matrix of eigenvectors
   */
  syms(const vecs<N, T>& lambda, const mats<N, N, T>& q) {
    from_eig(lambda, q);
  }

  /**
   * \brief Constructs this matrix from the eigendecomposition Q * D * Q^T
   *
   * \param[in] decomp (first = vector of n eigenvalues, second = n x n matrix
   * of eigenvectors)
   */
  syms(const std::pair<vecs<N, T>, mats<N, N, T> >& decomp)
      : syms(decomp.first, decomp.second) {}

  /**
   * \brief Constructs a symmetric matrix from a general matrix (which is
   * hopefully symmetric).
   *
   * \param[in] A - n x n matrix which is hopefully symmetric
   */
  syms(const mats<N, N, T>& A) {
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) (*this)(i, j) = A(i, j);
  }

  /**
   * \brief Constructs a matrix by copying values from some other symmetric
   * matrix.
   *
   * \param[in] A - n x n symmetric matrix with a potentially different type S
   */
  template <typename S>
  syms(const syms<N, S>& A) {
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) (*this)(i, j) = A(i, j);
  }

  /**
   * \brief Constructs a matrix by copying values from some other matrix
   * (hopefully symemtric)
   *
   * \param[in] A - n x n matrix (hopefully symmetric) with a potentially
   * different type S
   */
  template <typename S>
  syms(const mats<N, N, S>& A) {
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) (*this)(i, j) = A(i, j);
  }

  /**
   * \brief Constructs this matrix from the eigendecomposition Q * D * Q^T
   *
   * \param[in] lambda - vector of n eigenvalues
   * \param[in] q -  n x n matrix of eigenvectors
   */
  void from_eig(const vecs<N, T>& lambda, const mats<N, N, T>& q);

  /**
   * \brief Returns the number of unique entries stored in this matrix.
   *
   * \return number of unique entries
   */
  int nb() const { return NB; }

  /**
   * \brief Copies the entries from some other symmetric matrix into this
   * matrix.
   *
   * \param[in] T - symmetric matrix to copy
   */
  void copy(const syms<N, T>& A) {
    for (int k = 0; k < nb(); k++) data_[k] = A.data(k);
  }

  /**
   * \brief Copies the entries from some other matrix into this matrix.
   *
   * \param[in] T - matrix to copy (hopefully symmetric)
   */
  void set(const mats<N, N, T>& M) {
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) (*this)(i, j) = M(i, j);
  }

  /**
   * \brief Copies the entries from some other symmetric matrix into this
   * matrix.
   *
   * \param[in] T - symmetric matrix to copy
   */
  void set(const syms<N, T>& S) {
    for (int i = 0; i < nb(); i++) data_[i] = S.data(i);
  }

  /**
   * \brief Sets all entries in this symmetric matrix to zero.
   */
  void zero() {
    for (int i = 0; i < nb(); i++) data_[i] = 0.;
  }

  /**
   * \brief Sets this matrix to the identity matrix (diagonal entries = 1,
   * off-diagonals = 0)
   */
  void eye() {
    zero();
    for (int i = 0; i < N; i++) operator()(i, i) = 1.;
  }

  /**
   * \brief Read/write access to the entry into the flattened 1d array of
   * entries. This function should not be used often.
   *
   * \param[in] k - index of the item (in the flattened array) to retrieve.
   *
   * \return entry at index k
   */
  T& data(const int k) { return data_[k]; }

  /**
   * \brief Read access to the entry into the flattened 1d array of entries.
   *        This function should not be used often.
   *
   * \param[in] k - index of the item (in the flattened array) to retrieve.
   *
   * \return entry at index k
   */
  T data(const int k) const { return data_[k]; }

  /**
   * \brief Read/write access to entry (i,j)
   *
   * \param[in] i - row index
   * \param[in] j - column index
   *
   * \return reference to entry (i,j)
   */
  inline T& operator()(const uint8_t i, const uint8_t j) {
    ASSERT(i < N && j < N);
    return (i > j) ? data_[i * (i + 1) / 2 + j] : data_[j * (j + 1) / 2 + i];
  }

  /**
   * \brief Read access to entry (i,j)
   *
   * \param[in] i - row index
   * \param[in] j - column index
   *
   * \return reference to entry (i,j)
   */
  inline T operator()(const uint8_t i, const uint8_t j) const {
    ASSERT(i < N && j < N);
    return (i > j) ? data_[i * (i + 1) / 2 + j] : data_[j * (j + 1) / 2 + i];
  }

  /**
   * \brief Unary operator, so that symmetric matrices can be included in
   * expressions such as +A
   */
  syms operator+() const { return *this; }

  /**
   * \brief Matrix-matrix addition: (*this) + U
   *
   * \param[in] this - left operand in addition
   * \param[in] U - right operand in addition
   *
   * \return new symmetric matrix (*this) + U
   */
  syms<N, T> operator+(const syms<N, T>& U) {
    syms<N, T> V;
    for (int i = 0; i < U.nb(); i++) V.data(i) = data_[i] + U.data(i);
    return V;
  }

  /**
   * \brief Matrix-scalar multiplication.
   *
   * \param[in] this - matrix to multiply
   * \param[in] a - scalar
   *
   * \return (*this) * a
   */
  syms<N, T> operator*(const T a) const {
    syms<N, T> C;
    for (int i = 0; i < nb(); i++) C.data(i) = a * data_[i];
    return C;
  }

  /**
   * \brief Assignment operator from a scalar - sets all entries to a constant.
   *
   * \param a - scalar all entries will be set to
   *
   * \return (*this) where all entries are equal to a
   */
  syms<N, T>& operator=(const double& a) {
    for (int i = 0; i < nb(); i++) data_[i] = a;
    return *this;
  }

  /**
   * \brief Computes the eigendecomposition of a symmetric matrix,
   *        i.e. M = Q * D * Q^T, returns Q * exp(D) * Q^T,
   *        where D is a diagonal matrix with the eigenvalues
   *        and Q are the eigenvectors.
   *
   * \param[in] M - the symmetric d x d matrix
   *
   * \return pair {D,Q} where D is a d-vector of diagonal entries and Q  is a d
   * x d matrix of eigenvectors
   */
  std::pair<vecs<N, T>, mats<N, N, T> > eig() const;

  /**
   * \brief Prints this matrix to the console.
   *
   * \param[in] title (optional) prefix for printing entries
   *
   */
  void print(const std::string& title = std::string()) const;

 private:
  T data_[NB];  // flattened 1d array used to store matrix entries

  /**
   * \brief Computes the eigendecomposition of a symmetric matrix using the
   * method described here: https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/
   *        i.e. M = Q * D * Q^T, returns Q * exp(D) * Q^T,
   *        where D is a diagonal matrix with the eigenvalues
   *        and Q are the eigenvectors.
   *
   * \param[in] M - the symmetric d x d matrix
   *
   * \return pair {D,Q} where D is a d-vector of diagonal entries and Q  is a d
   * x d matrix of eigenvectors
   */
  std::pair<vecs<N, T>, mats<N, N, T> > __eig__() const;
};

/**
 * \brief Computes symmetric matrix-matrix multiplication A * B
 */
template <typename R, typename S>
matd<typename result_of<R, S>::type> operator*(const symd<R>& A,
                                               const symd<S>& B);

/**
 * \brief Computes symmetric matrix addition A + B
 */
template <typename R, typename S>
symd<typename result_of<R, S>::type> operator+(const symd<R>& A,
                                               const symd<S>& B);

/**
 * \brief Computes symmetric matrix subtraction A - B
 */
template <typename R, typename S>
symd<typename result_of<R, S>::type> operator-(const symd<R>& A,
                                               const symd<S>& B);

/**
 * \brief Computes matrix addition A - B when A is a symmetric matrix and B is a
 * general matrix
 */
template <typename R, typename S>
matd<typename result_of<R, S>::type> operator-(const symd<R>& A,
                                               const matd<S>& B);

/**
 * \brief Computes matrix addition A - B when A is a general matrix and B is a
 * symmetric matrix
 */
template <typename R, typename S>
matd<typename result_of<R, S>::type> operator-(const matd<R>& A,
                                               const symd<S>& B);

/**
 * \brief Computes symmetric matrix-matrix multiplication A * B
 */
template <int N, typename R, typename S>
mats<N, N, typename result_of<R, S>::type> operator*(const syms<N, R>& A,
                                                     const syms<N, S>& B);

/**
 * \brief Computes symmetric matrix addition A + B
 */
template <int N, typename R, typename S>
syms<N, typename result_of<R, S>::type> operator+(const syms<N, R>& A,
                                                  const syms<N, S>& B);

/**
 * \brief Computes symmetric matrix subtraction A - B
 */
template <int N, typename R, typename S>
syms<N, typename result_of<R, S>::type> operator-(const syms<N, R>& A,
                                                  const syms<N, S>& B);

/**
 * \brief Computes matrix addition A - B when A is a symmetric matrix and B is a
 * general matrix
 */
template <int N, typename R, typename S>
mats<N, N, typename result_of<R, S>::type> operator-(const syms<N, R>& A,
                                                     const mats<N, N, S>& B);

/**
 * \brief Computes matrix addition A - B when A is a general matrix and B is a
 * symmetric matrix
 */
template <int N, typename R, typename S>
mats<N, N, typename result_of<R, S>::type> operator-(const mats<N, N, R>& A,
                                                     const syms<N, S>& B);

}  // namespace vortex
