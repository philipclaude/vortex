#pragma once

#include <fmt/format.h>

#include <iostream>
#include <vector>

#include "log.h"
#include "result_of.h"

namespace vortex {

template <typename T>
class vecd;
template <typename T>
class symd;
template <int N, typename T>
class syms;
template <int N, typename T>
class vecs;

/**
 * \brief Represents a dynamically-allocated m x n matrix with m rows and n
 * columns, where the entries all have some type 'T'. Entries are stored in
 * column-major order.
 *
 *  Supported operators (see at the bottom of this file and in mat.hpp):
 *  matrix addition: A + B
 *  matrix subtraction: A - B
 *  multiplication (upper-case for matrices, lower-case for scalars): A * B , A
 * * b , a * B matrix-vector multiplication: A * x where x is a vector unary
 * operators: +A, -A
 */
template <typename T>
class matd {
 public:
  /**
   * \brief Constructs a square n x n matrix.
   *
   * \param[in] n - number of rows/columns.
   */
  matd(int n) : m_(n), n_(n), data_(n * n, 0) {}

  /**
   * \brief Constructs an empty matrix. The function 'resize' can be used to set
   * the sizes later.
   */
  matd() : m_(0), n_(0) {}

  /**
   * \brief Constructs a rectangular m x n matrix.
   *
   * \param[in] m - number of rows
   * \param[in] n - number of columns
   */
  matd(int m, int n) : m_(m), n_(n), data_(m * n, 0) {}

  /**
   * \brief Constructs a rectangular matrix from another matrix with a
   *        potentially different type than T.
   *
   * \param[in] A - the m x n matrix to copy where the entries have a type S
   */
  template <typename S>
  matd(const matd<S>& A) : m_(A.m()), n_(A.n()), data_(m_ * n_) {
    for (int i = 0; i < m_; i++)
      for (int j = 0; j < n_; j++) (*this)(i, j) = A(i, j);
  }

  /**
   * \brief Constructs a rectangular matrix from a symmetric matrix with a
   *        potentially different type than T.
   *
   * \param[in] A - the m x n symmetric matrix to copy where the entries have a
   * type S
   */
  template <typename S>
  matd(const symd<S>& A);

  /**
   * \brief Read/write access to entry (i,j)
   *
   * \param[in] i - row
   * \param[in] j - column
   *
   * \return Aij
   */
  inline T& operator()(int i, int j) {
    ASSERT(i < m_ && j < n_)
        << fmt::format("i = {}, j = {} of {} x {} matrix", i, j, m_, n_);
    return data_[j * m_ + i];
  }

  /**
   * \brief Read access to entry (i,j)
   *
   * \param[in] i - row
   * \param[in] j - column
   *
   * \return Aij
   */
  inline const T& operator()(int i, int j) const {
    ASSERT(i < m_ && j < n_)
        << fmt::format("i = {}, j = {} of {} x {} matrix", i, j, m_, n_);
    return data_[j * m_ + i];
  }

  /**
   * \brief Returns the number of rows in the matrix.
   */
  int m() const { return m_; }

  /**
   * \brief Returns the number of columns in the matrix.
   */
  int n() const { return n_; }

  /**
   * \brief Converts this matrix to a scalar double.
   */
  operator double() const {
    ASSERT(m_ == 1 && n_ == 1) << "matrix must be 1x1 to convert to scalar";
    return data_[0];
  }

  /**
   * \brief Converts this matrix to a scalar float.
   */
  operator float() const {
    ASSERT(m_ == 1 && n_ == 1) << "matrix must be 1x1 to convert to scalar";
    return data_[0];
  }

  /**
   * \brief Sets all matrix entries to zero.
   */
  void zero() {
    for (int i = 0; i < m_ * n_; i++) data_[i] = 0;
  }

  /**
   * \brief Sets this (square) matrix to the identity matrix, with ones on the
   *        diagonal, and zeros on the off-diagonals.
   */
  void eye() {
    ASSERT(m_ == n_);
    zero();
    for (int i = 0; i < m_; i++) (*this)(i, i) = 1;
  }

  /**
   * \brief Sets all entries to a constant.
   *
   * \param[in] x - the value each entry should be set to.
   */
  void set(const double& x) {
    for (size_t j = 0; j < data_.size(); j++) data_[j] = x;
  }

  /**
   * \brief Copies value from another matrix.
   *
   * \param[in] A - m x n matrix to copy (should have the same size as this
   * matrix).
   */
  void set(const matd<T>& A) {
    ASSERT(A.m() == m_ && A.n() == n_);
    for (int i = 0; i < m_; i++)
      for (int j = 0; j < n_; j++) (*this)(i, j) = A(i, j);
  }

  /**
   * \brief Copies value from a symmetric matrix.
   *
   * \param[in] S - m x m symmetric matrix to copy.
   */
  void set(const symd<T>& S);

  /**
   * \brief Sets the sizes of the m x n marix, and re-allocates storage for the
   * entries.
   *
   * \param[in] m - number of rows
   * \param[in] n - number of columns
   */
  void resize(int m, int n) {
    m_ = m;
    n_ = n;
    data_.resize(m * n);
  }

  /**
   * \brief Sets a particular row to a vector of values.
   *
   * \param[in] i - the row to set
   * \param[in] row - the values to set into row i
   */
  void set_row(int i, const vecd<T>& row);

  /**
   * \brief Retrieves a particular row of the matrix.
   *
   * \param[in] i - the row to retrieve
   * \param[out] row - where the values will be written to (the size should be
   * the number of columns in this matrix).
   */
  void get_row(int i, vecd<T>& row) const;

  /**
   * \brief Prints out the values in the matrix.
   *
   * \param[in] v0 (optional) - label for each entry
   */
  void print(const std::string& v0 = std::string()) const {
    std::string v = v0.empty() ? "" : v0;
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    for (int i = 0; i < m_; i++)
      for (int j = 0; j < n_; j++)
        std::cout << v
                  << "(" + std::to_string(i) + "," + std::to_string(j) + "): "
                  << (*this)(i, j) << std::endl;
  }

 protected:
  int m_;                // number of rows
  int n_;                // number of columns
  std::vector<T> data_;  // matrix entries, stored in column-major order
};

/**
 * \brief Represents a statically-allocated M x N matrix with M rows and N
 * columns, where the entries all have some type 'T'. M and N must be known at
 * compile time. Entries are stored in column-major order.
 *
 * Supported operators (see at the bottom of this file and in mat.hpp):
 *  matrix addition: A + B
 *  matrix subtraction: A - B
 *  multiplication (upper-case for matrices, lower-case for scalars): A * B , A
 * * b , a * B matrix-vector multiplication: A * x where x is a vector unary
 * operators: +A, -A
 */
template <int M, int N, typename T>
class mats {
 public:
  /**
   * \brief Constructs a rectangular matrix and sets all entries to zero.
   */
  mats() { zero(); }

  /**
   * \brief Copies an M x N matrix A of a potentially different type.
   *
   * \param[in] A - M x N matrix whose entries have a type S
   */
  template <typename S>
  mats(const mats<M, N, S>& A) {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++) (*this)(i, j) = A(i, j);
  }

  /**
   * \brief Constructs a rectangular matrix from a symmetric matrix with a
   *        potentially different type than T.
   *
   * \param[in] A - the m x n symmetric matrix to copy where the entries have a
   * type S
   */
  template <typename S>
  mats(const syms<M, S>& A);

  /**
   * \brief Sets all matrix entries to zero.
   */
  void zero() {
    for (int i = 0; i < M * N; i++) data_[i] = 0;
  }

  /**
   * \brief Sets this (square) matrix to the identity matrix, with ones on the
   *        diagonal, and zeros on the off-diagonals.
   */
  void eye() {
    ASSERT(M == N);
    zero();
    for (int i = 0; i < M; i++) (*this)(i, i) = 1;
  }

  /**
   * \brief Copies the values from some matrix to this matrix.
   *
   * \param[in] b - M x N matrix to copy whose values have a tyep S
   */
  template <typename S>
  mats<M, N, T>& operator=(const mats<M, N, S>& b) {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++) (*this)(i, j) = b(i, j);
    return *this;
  }

  /**
   * \brief Sets all entries to some scalar (double).
   *
   * \param[in] b - scalar value to assign to each entry.
   */
  mats<M, N, T>& operator=(const double& b) {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++) (*this)(i, j) = b;
    return *this;
  }

  /**
   * \brief Read/write access to entry (i,j)
   *
   * \param[in] i - row
   * \param[in] j - column
   *
   * \return Aij
   */
  T& operator()(int i, int j) {
    ASSERT(i < M && j < N);
    return data_[j * M + i];
  }

  /**
   * \brief Read access to entry (i,j)
   *
   * \param[in] i - row
   * \param[in] j - column
   *
   * \return Aij
   */
  const T& operator()(int i, int j) const {
    ASSERT(i < M && j < N);
    return data_[j * M + i];
  }

  /**
   * \brief Prints out the values in the matrix.
   */
  void print() const {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        std::cout << "(" + std::to_string(i) + "," + std::to_string(j) + "): "
                  << (*this)(i, j) << std::endl;
  }

  /**
   * \brief Returns a read-only pointer to the matrix data.
   *
   * \return const pointer to data
   */
  const T* data() const { return data_; }

  /**
   * \brief Converts this matrix to a scalar (as long as it is 1x1).
   *        Useful when computing something like x^T * A * x
   */
  operator T() const {
    ASSERT(M == 1 && N == 1);
    return data_[0];
  }

 protected:
  T data_[M * N];
};

/**
 * \brief Allows to print the matrix to an output stream, such as std::cout
 *
 * \param[in] os - output stream to write to
 * \param[in] z - the matrix to print
 */
template <int M, int N, class T>
std::ostream& operator<<(std::ostream& os, const mats<M, N, T>& z) {
  z.print();
  return os;
}

/**
 * \brief Computes and returns C = A + B
 */
template <int M, int N, class R>
mats<M, N, R> operator+(const mats<M, N, R>& A, const mats<M, N, R>& B) {
  mats<M, N, R> C;  // initializes to zero
  for (int i = 0; i < M; i++)
    for (int j = 0; j < N; j++) C(i, j) = A(i, j) + B(i, j);
  return C;
}

/**
 * \brief Computes the matrix-vector multiplication A * x.
 */
template <typename T>
vecd<T> operator*(const matd<T>& A, const vecd<T>& x);

/**
 * \brief Computes the matrix-matrix multiplication A * B.
 */
template <typename R, typename S>
matd<typename result_of<R, S>::type> operator*(const matd<R>& A,
                                               const matd<S>& B);

/**
 * \brief Computes the matrix-scalar multiplication A * b.
 */
template <typename R, typename S>
matd<typename result_of<R, S>::type> operator*(const matd<R>& A, const S& b);

/**
 * \brief Compute the scalar-matrix multiplication a * B.
 */
template <typename R, typename S>
matd<typename result_of<R, S>::type> operator*(const S& a, const matd<R>& B);

/**
 * \brief Computes the matrix addition A + B.
 */
template <typename R, typename S>
matd<typename result_of<R, S>::type> operator+(const matd<R>& A,
                                               const matd<S>& B);

/**
 * \brief Computes the matrix subtraction A - B.
 */
template <typename R, typename S>
matd<typename result_of<R, S>::type> operator-(const matd<R>& A,
                                               const matd<S>& B);

/**
 * \brief Unary + operator, for including A in expressions such as +A + B.
 */
template <typename R>
matd<R> operator+(const matd<R>& A);

/**
 * \brief Unary - operator, for including A in expressions such as -A + B.
 */
template <typename R>
matd<R> operator-(const matd<R>& A);

/**
 * \brief Unary + operator, for including A in expressions such as +A + B.
 */
template <int M, int N, typename R>
mats<M, N, R> operator+(const mats<M, N, R>& A);

/**
 * \brief Unary - operator, for including A in expressions such as -A + B.
 */
template <int M, int N, typename R>
mats<M, N, R> operator-(const mats<M, N, R>& A);

/**
 * \brief Specialized definition for 4x4 matrices.
 */
typedef mats<4, 4, double> mat4d;
typedef mats<4, 4, float> mat4f;
typedef mats<3, 3, float> mat3f;

}  // namespace vortex
