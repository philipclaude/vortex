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
  spmat(int m, int n, int bandwidth = 10) {
    rows_.reserve(m);
    cols_.reserve(n);
    triplets_.reserve(std::max(m, n) * bandwidth);
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
    typename std::unordered_map<std::pair<int, int>, T>::iterator it =
        triplets_.find({i, j});
    if (it == triplets_.end()) {
      triplets_.insert({{i, j}, 0});  // initialize the entry to zero
      rows_.insert(i);
      cols_.insert(j);
    }
    it = triplets_.find({i, j});
    ASSERT(it != triplets_.end());
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
    ASSERT(triplets_.find({i, j}) != triplets_.end());
    return triplets_.at({i, j});
  }

  /**
   * \brief Solves A * x = b using OpenNL (recommended).
   *
   * \param[in]    b - right-hand side vector
   * \param[inout] x - solution vector
   * \param[in]    symmetric - option to specify if the matrix is symmetric
   */
  void solve_nl(const vecd<T>& b, vecd<T>& x, bool symmetric = true) const;

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
   * \brief Returns the number of rows in the matrix, according to the triplets
   * added.
   */
  int nb_rows() const { return rows_.size(); }

  /**
   * \brief Returns the number of columns in the matrix, according to the
   * triplets added.
   */
  int nb_cols() const { return cols_.size(); }

  /**
   * \brief Returns the number of nonzero entries in the matrix.
   */
  int nb_nnz() const { return triplets_.size(); }

  /**
   * \brief Returns the triplets (row,col,value) stored in the matrix.
   */
  const std::unordered_map<std::pair<int, int>, T>& triplets() const {
    return triplets_;
  }

  /**
   * \brief Prints all the triplets (row,col,value) stored in the matrix.
   */
  void print() const {
    for (auto& t : triplets_) {
      std::cout << "A(" << t.first.first << "," << t.first.second
                << ") = " << t.second << std::endl;
    }
  }

  /**
   * \brief Prints the dense representation of the matrix (not recommended for
   * large matrices).
   */
  void print_full() const {
    for (int i = 0; i < this->nb_rows(); i++) {
      for (int j = 0; j < this->nb_cols(); j++) {
        double value = 0.0;
        auto it = triplets_.find({i, j});
        if (it != triplets_.end()) value = it->second;
        std::cout << value << " ";
      }
      std::cout << std::endl;
    }
  }

 private:
  std::unordered_set<int> rows_;
  std::unordered_set<int> cols_;
  std::unordered_map<std::pair<int, int>, T> triplets_;  // (row,col,value)
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
  b.zero();
  auto& triplets = A.triplets();
  for (auto& t : triplets) {
    int row = t.first.first;
    int col = t.first.second;
    const T& entry = t.second;
    b(row) += entry * x(col);
  }
  return b;
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
