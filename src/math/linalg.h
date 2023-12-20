#pragma once

#include <cmath>
#include <vector>

namespace vortex {

// forward declarations
template <int M, int N, typename T> class mats;
template <typename type> class matd;
template <typename T> class symd;
template <int N, typename T> class syms;
template <typename T> class vecd;
template <int N, typename T> class vecs;

/*\
 * =============================================================================
 *
 * mats
 *
 * =============================================================================
\*/
/**
 * \brief Computes the transpose of a matrix.
 *
 * \param[in] A - the M x N matrix to transpose
 *
 * \return A^T, a N x M matrix
 */
template <int M, int N, typename T>
mats<N, M, T> transpose(const mats<M, N, T>& A);

/**
 * \brief Computes the determinant of a square matrix.
 *
 * \param[in] A - the M x M matrix to compute the determinant of
 *
 * \return det(A), the determinant of A
 */
template <int M, typename T> T det(const mats<M, M, T>& A);

/**
 * \brief Computes the trace of a square matrix (sum of diagonal entries).
 *
 * \param[in] A - the M x M matrix to compute the trace of
 *
 * \return tr(A), the trace of A
 */
template <int M, typename T> T trace(const mats<M, M, T>& A);

/**
 * \brief Computes the inverse of a square matrix.
 *
 * \param[in] A - the M x M matrix to invert
 *
 * \return A^(-1), a M x M matrix
 */
template <int M, typename T> mats<M, M, T> inverse(const mats<M, M, T>& A);

/**
 * \brief Construct a diagonal matrix from a vector of diagonal entries.
 *        Off-diagonal entries will be zero.
 *
 * \param[in] d - the d.m() vector of diagonal entries.
 *
 * \return a N x N matrix
 */
template <int N, typename T> mats<N, N, T> diag(const vecs<N, T>& d);

/*\
 * =============================================================================
 *
 * matd
 *
 * =============================================================================
\*/
/**
 * \brief Computes the transpose of a matrix.
 *
 * \param[in] A - the A.m() x A.n() matrix to transpose
 *
 * \return a A^T A.n() x A.m() matrix
 */
template <typename T> matd<T> transpose(const matd<T>& A);

/**
 * \brief Construct a diagonal matrix from a vector of diagonal entries.
 *        Off-diagonal entries will be zero.
 *
 * \param[in] d - the d.m() vector of diagonal entries.
 *
 * \return a d.m() x d.m() matrix
 */
template <typename T> matd<T> diag(const vecd<T>& d);

/**
 * \brief Solve the system of equations A * x = b using LU factorization with
 * pivoting.
 *
 * \param[in] A - the m x m matrix with coefficients of the linear system
 * \param[in] b - the m-dimensional vector defining the right-hand side of the
 * system \param[out] x - the solution to the linear system
 */
template <typename T>
void solveLUP(const matd<T>& A, const vecd<T>& b, vecd<T>& x);

/**
 * \brief Computes the inverse of a square matrix using LU factorization with
 * pivoting.
 *
 * \param[in] A - the m x m matrix to invert
 * \param[out] Ainv - the resulting inverse A^(-1) (m x m)
 */
template <typename T> void inverseLUP(const matd<T>& A, matd<T>& Ainv);

/**
 * \brief Computes and returns the inverse of a square matrix.
 *
 * \param[in] A - the m x m matrix to invert
 *
 * \return A^(-1), the resulting inverse (m x m)
 */
template <typename T> matd<T> inverse(const matd<T>& M);

/**
 * \brief Computes the determinant of a square matrix.
 *
 * \param[in] A - the m x m matrix to compute the determinant of
 *
 * \return det(A), the determinant of A
 */
template <typename T> T det(const matd<T>& A);

/*\
 * =============================================================================
 *
 * symd
 *
 * =============================================================================
\*/
/**
 * \brief Computes the exponential of a symmetric matrix.
 *        From the eigendcomposition of M = Q * D * Q^T, returns Q * exp(D) *
 * Q^T, where exp(D) is the diagonal matrix whose entries are the exponent of
 * the eigenvalues and Q are the eigenvectors.
 *
 * \param[in] M - the symmetric d x d matrix
 *
 * \return exp(M) = Q * exp(D) * Q^T
 */
template <typename type> symd<type> expm(const symd<type>& M);

/**
 * \brief Computes the logarithm of a symmetric matrix.
 *        From the eigendcomposition of M = Q * D * Q^T, returns Q * log(D) *
 * Q^T, where log(D) is the diagonal matrix whose entries are the log of the
 * eigenvalues and Q are the eigenvectors.
 *
 * \param[in] M - the symmetric d x d matrix
 *
 * \return exp(M) = Q * log(D) * Q^T
 */
template <typename type> symd<type> logm(const symd<type>& M);

/**
 * \brief Computes the exponent of a symmetric matrix.
 *        From the eigendcomposition of M = Q * D * Q^T, returns Q * pow(D,p) *
 * Q^T, where pow(D,p) is the diagonal matrix whose entries are the eigenvalues
 * raised to the p'th power and Q are the eigenvectors.
 *
 * \param[in] M - the symmetric d x d matrix
 * \param[in] p - the exponent
 *
 * \return exp(M) = Q * pow(D,p) * Q^T
 */
template <typename type> symd<type> powm(const symd<type>& M, double p);

/**
 * \brief Computes the square-root of a symmetric matrix.
 *        From the eigendcomposition of M = Q * D * Q^T, returns Q * sqrt(D) *
 * Q^T, where sqrt(D) is the diagonal matrix whose entries are the square-root
 * of the eigenvalues and Q are the eigenvectors.
 *
 * \param[in] M - the symmetric d x d matrix
 *
 * \return exp(M) = Q * sqrt(D) * Q^T
 */
template <typename type> symd<type> sqrtm(const symd<type>& M);

/**
 * \brief Computes the determinant of a symmetric matrix.
 *
 * \param[in] M - the symmetric d x d matrix
 *
 * \return det(M)
 */
template <typename type> type det(const symd<type>& M);

/**
 * \brief Computes the inverse of a symmetric matrix.
 *
 * \param[in] M - the symmetric d x d matrix
 *
 * \return M^(-1), the d x d inverse of M
 */
template <typename T> symd<T> inverse(const symd<T>& M);

/**
 * \brief Computes the eigendecomposition of a symmetric matrix,
 *        i.e. M = Q * D * Q^T, returns Q * exp(D) * Q^T,
 *        where D is a diagonal matrix with the eigenvalues
 *        and Q are the eigenvectors.
 *
 * \param[in] M - the symmetric d x d matrix
 * \param[out] D - d-vector of diagonal entries
 * \param[out] Q - d x d matrix with eigenvectors
 */
template <typename T> void eig(const symd<T>& m, vecd<T>& D, matd<T>& Q);

/**
 * \brief Computes the eigendecomposition of a symmetric matrix,
 *        i.e. M = Q * D * Q^T, returns Q * exp(D) * Q^T,
 *        where D is a diagonal matrix with the eigenvalues
 *        and Q are the eigenvectors.
 *
 * \param[in] M - the symmetric d x d matrix

 * \return a pair {D,Q} with first (D) being the d-vector of diagonal entries
           and second (Q) being the d x d matrix of eigenvectors.
 */
template <typename T> std::pair<vecd<T>, matd<T>> eig(const symd<T>& m);

/**
 * \brief Computes the Log-Euclidean weighted average of a set of symmetric
 * matrices.
 *
 * \param[in] alpha - vector of weights
 * \param[in] tensors - the set of tensors
 *
 * \return exp( \sum_i log(tensors_i) * alpha_i )
 */
template <typename T>
symd<T> interp(const std::vector<double>& alpha,
               const std::vector<symd<T>>& tensors);

/*\
 * =============================================================================
 *
 * syms
 *
 * =============================================================================
\*/
/**
 * \brief Computes the exponential of a symmetric matrix.
 *        From the eigendcomposition of M = Q * D * Q^T, returns Q * exp(D) *
 * Q^T, where exp(D) is the diagonal matrix whose entries are the exponent of
 * the eigenvalues and Q are the eigenvectors.
 *
 * \param[in] M - the symmetric d x d matrix
 *
 * \return exp(M) = Q * exp(D) * Q^T
 */
template <int N, typename T> syms<N, T> expm(const syms<N, T>& M);

/**
 * \brief Computes the logarithm of a symmetric matrix.
 *        From the eigendcomposition of M = Q * D * Q^T, returns Q * log(D) *
 * Q^T, where log(D) is the diagonal matrix whose entries are the log of the
 * eigenvalues and Q are the eigenvectors.
 *
 * \param[in] M - the symmetric d x d matrix
 *
 * \return exp(M) = Q * log(D) * Q^T
 */
template <int N, typename T> syms<N, T> logm(const syms<N, T>& M);

/**
 * \brief Computes the exponent of a symmetric matrix.
 *        From the eigendcomposition of M = Q * D * Q^T, returns Q * pow(D,p) *
 * Q^T, where pow(D,p) is the diagonal matrix whose entries are the eigenvalues
 * raised to the p'th power and Q are the eigenvectors.
 *
 * \param[in] M - the symmetric d x d matrix
 * \param[in] p - the exponent
 *
 * \return exp(M) = Q * pow(D,p) * Q^T
 */
template <int N, typename T> syms<N, T> powm(const syms<N, T>& M, double p);

/**
 * \brief Computes the square-root of a symmetric matrix.
 *        From the eigendcomposition of M = Q * D * Q^T, returns Q * sqrt(D) *
 * Q^T, where sqrt(D) is the diagonal matrix whose entries are the square-root
 * of the eigenvalues and Q are the eigenvectors.
 *
 * \param[in] M - the symmetric d x d matrix
 *
 * \return exp(M) = Q * sqrt(D) * Q^T
 */
template <int N, typename T> syms<N, T> sqrtm(const syms<N, T>& M);

/**
 * \brief Computes the determinant of a symmetric matrix.
 *
 * \param[in] M - the symmetric d x d matrix
 *
 * \return det(M)
 */
template <int N, typename T> T det(const syms<N, T>& M);

/**
 * \brief Computes the inverse of a symmetric matrix.
 *
 * \param[in] M - the symmetric d x d matrix
 *
 * \return M^(-1), the d x d inverse of M
 */
template <int N, typename T> syms<N, T> inverse(const syms<N, T>& M);

/**
 * \brief Computes the eigendecomposition of a symmetric matrix,
 *        i.e. M = Q * D * Q^T, returns Q * exp(D) * Q^T,
 *        where D is a diagonal matrix with the eigenvalues
 *        and Q are the eigenvectors.
 *
 * \param[in] M - the symmetric d x d matrix
 * \param[out] D - d-vector of diagonal entries
 * \param[out] Q - d x d matrix with eigenvectors
 */
template <int N, typename T>
void eig(const syms<N, T>& M, vecs<N, T>& D, mats<N, N, T>& Q);

/**
 * \brief Computes the eigendecomposition of a symmetric matrix,
 *        i.e. M = Q * D * Q^T, returns Q * exp(D) * Q^T,
 *        where D is a diagonal matrix with the eigenvalues
 *        and Q are the eigenvectors.
 *
 * \param[in] M - the symmetric d x d matrix

 * \return a pair {D,Q} with first (D) being the d-vector of diagonal entries
           and second (Q) being the d x d matrix of eigenvectors.
 */
template <int N, typename T>
std::pair<vecs<N, T>, mats<N, N, T>> eig(const syms<N, T>& M);

/**
 * \brief Computes the Log-Euclidean weighted average of a set of symmetric
 * matrices.
 *
 * \param[in] alpha - vector of weights
 * \param[in] tensors - the set of tensors
 *
 * \return exp( \sum_i log(tensors_i) * alpha_i )
 */
template <int N, typename T>
syms<N, T> interp(const std::vector<double>& alpha,
                  const std::vector<syms<N, T>>& tensors);

}  // namespace vortex
