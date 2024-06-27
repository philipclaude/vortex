#pragma once
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include "log.h"
#include "stlext.h"
#include "vec.h"

namespace vortex {

/*
 * Represents a sparse matrix in which the values stored are any type T.
 */
template <typename T>
class spmat {
 public:
  /**
   * \brief Initializes a m x n sparse matrix, reserving space for the non-zero
   * entries.
   *
   * \param[in] m - number of rows
   * \param[in] n - number of columns
   * \param[in] bandwidth (optional) - estimate on the number of non-zero
   * entries per row or column (default: 10)
   */
  spmat(int m, int n, int bandwidth = 10) : n_rows_(m), n_cols_(n) {
    rows_.resize(m);
    for (auto& row : rows_) row.reserve(bandwidth);
  }

  /**
   * \brief Read/write access to an entry in the matrix. Adds a triplet if
   * necessary.
   *
   * \param[in] i - row index
   * \param[in] j - column index
   *
   * \return reference to the entry at (i,j)
   */
  T& operator()(int i, int j) {
    auto it = rows_[i].find(j);
    if (it == rows_[i].end()) {
      rows_[i].insert({j, 0});  // initialize the entry to zero
    }
    it = rows_[i].find(j);
    ASSERT(it != rows_[i].end());
    return it->second;
  }

  /**
   * \brief Read access to an entry in the matrix.
   *
   * \param[in] i - row index
   * \param[in] j - column index
   *
   * \return const reference to the entry at (i,j)
   */
  const T& operator()(int i, int j) const {
    ASSERT(rows_[i].find(j) != rows_[i].end());
    return rows_[i].at(j);
  }

  /**
   * \brief Solves A * x = b using OpenNL (recommended).
   *
   * \param[in]    b - right-hand side vector
   * \param[inout] x - solution vector
   * \param[in]    symmetric - option to specify if the matrix is symmetric
   */
  void solve_nl(const vecd<T>& b, vecd<T>& x, double tol,
                bool symmetric = true) const;

  /**
   * \brief Solves A * x = b using the Jacobi method.
   *
   * \param[in]    b - right-hand side vector
   * \param[inout] x - solution vector
   * \param[in]    tol - tolerance to use in convergence test
   * \param[in]    max_iter - maximum number of iterations
   * \param[in]    verbose - whether to print convergence information at each
   * iteration
   */
  double solve_jacobi(const vecd<T>& b, vecd<T>& x, double tol,
                      int max_iter = -1, bool verbose = false) const;

  /**
   * \brief Returns the number of rows in the matrix.
   */
  int n_rows() const { return rows_.size(); }

  /**
   * \brief Returns the number of columns in the matrix.
   */
  int n_cols() const { return n_cols_; }

  /**
   * \brief Returns the number of nonzero entries in the matrix.
   */
  int nnz() const {
    size_t count = 0;
    for (const auto& row : rows_) count += row.size();
    return count;
  }

  /**
   * \brief Prints all the triplets (row,col,value) stored in the matrix.
   */
  void print() const {
    for (size_t row = 0; row < rows_.size(); row++) {
      for (const auto& [col, value] : rows_[row])
        std::cout << "(" << row << ", " << col << "):" << value << std::endl;
    }
  }

  /**
   * \brief Prints the dense representation of the matrix (not recommended for
   * large matrices).
   */
  void print_full() const {
    if (rows_.size() > 50) {
      std::cout << "matrix is too big" << std::endl;
      return;
    }
    for (size_t row = 0; row < rows_.size(); row++) {
      for (size_t col = 0; col < n_cols_; col++) {
        double value = 0.0;
        auto it = rows_[row].find(col);
        if (it != rows_[row].end()) value = it->second;
        std::cout << value << " ";
      }
      std::cout << std::endl;
    }
  }

  const auto& rows() const { return rows_; }
  auto& rows() { return rows_; }

  void clear() {
    for (auto& row : rows_) row.clear();
  }

 private:
  std::vector<std::unordered_map<uint32_t, T>> rows_;
  size_t n_rows_{0};
  size_t n_cols_{0};
};

/**
 * \brief Computes sparse matrix-vector multiplication b = A * x (only considers
 * non-zero entries).
 *
 * \param[in] A - sparse matrix
 * \param[in] x - vector
 *
 * \return b = A * x
 */
template <typename T>
vecd<T> operator*(const spmat<T>& A, const vecd<T>& x) {
  vecd<T> b(x.m());
  spmv(A, x, b);
  return b;
}

/**
 * \brief Computes sparse matrix-vector multiplication b = A * x.
 *
 * \param[in] A - sparse matrix
 * \param[in] x - vector
 *
 * \return b = A * x
 */
template <typename T>
void spmv(const spmat<T>& A, const vecd<T>& x, vecd<T>& b) {
  b.zero();
  const auto& rows = A.rows();
  for (size_t row = 0; row < rows.size(); row++) {
    for (const auto& [col, entry] : rows[row]) b(row) += entry * x(col);
  }
}

/**
 * \brief Computes sparse matrix-matrix multiplication C = A * B.
 *
 * \param[in] A - sparse matrix
 * \param[in] B - sparse matrix
 * \param[out] C - sparse matrix
 */
template <typename T>
void spmm(const spmat<T>& A, const spmat<T>& B, spmat<T>& C) {
  ASSERT(A.n_cols() == B.n_rows());
  ASSERT(C.n_rows() == A.n_rows());
  ASSERT(C.n_cols() == B.n_cols());
  const auto& a_rows = A.rows();
  const auto& b_rows = B.rows();
  for (size_t a_row = 0; a_row < a_rows.size(); a_row++) {
    for (const auto& [a_col, a_value] : a_rows[a_row]) {
      for (const auto& [b_col, b_value] : b_rows[a_col])
        C(a_row, b_col) += a_value * b_value;
    }
  }
}

/**
 * \brief Computes the transpose of a matrix, optionally multiplying factors by
 * some constant.
 *
 * \param[in] A - sparse matrix
 * \param[out] At - tranpose of A (sparse matrix)
 * \param[in] f - factor to multiply the entries by
 */
template <typename T>
void transpose(const spmat<T>& A, spmat<T>& At, double f = 1) {
  ASSERT(At.n_rows() == A.n_cols());
  ASSERT(At.n_cols() == A.n_rows());
  const auto& a_rows = A.rows();
  for (size_t i = 0; i < a_rows.size(); i++) {
    for (const auto& [j, value] : a_rows[i]) At(j, i) = f * value;
  }
}

/**
 * \brief Calculates the norm of a vector. This should probably be in vec.hpp.
 *
 * \param[in] x - vector
 *
 * \return || x ||
 */
template <typename T>
T norm(const vecd<T>& x) {
  T result = 0;
  for (int k = 0; k < x.m(); k++) result += x(k) * x(k);
  return std::sqrt(result);
}

}  // namespace vortex
