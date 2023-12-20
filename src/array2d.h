#pragma once

#include <fmt/format.h>

#include <iostream>
#include <vector>

#include "defs.h"
#include "log.h"

namespace vortex {

enum LayoutCategory {
  Layout_Jagged,      // size along second dimension is not constant
  Layout_Rectangular  // size along second dimension is constant
};

/**
 * \brief Represents a two-dimensional array, but is stored using
 *        one-dimensional arrays. Can represent rectangular arrays,
 *        in which the size along the second dimension is constant,
 *        or jagged arrays, in which the size along the second dimension
 *        is different for each element.
 */
template <typename T> class array2d {
 public:
  /**
   * \brief Constructs a 2d array, setting the layout and stride
   *
   * \param[in] stride - number of elements in along the second dimension
   *                     stride > 0 for rectangular arrays
   *                     stride < 0 for jagged arrays
   */
  array2d(int stride)
      : layout_((stride) > 0 ? Layout_Rectangular : Layout_Jagged),
        stride_(stride) {}

  /**
   * \brief Constructs a 2d array, setting the layout and stride and
   *        allocating the array, while setting all elements to a value.
   *
   * \param[in] stride - number of elements in along the second dimension
   *                     stride > 0 for rectangular arrays
   *                     stride < 0 for jagged arrays
   * \param[in] N - number of elements along the first dimension
   * \param[in] value - the value all (N x stride) components will be set to.
   */
  array2d(int stride, int N, T value)
      : layout_(Layout_Rectangular), stride_(stride) {
    ASSERT(stride > 0);
    std::vector<T> x(stride, value);
    for (int k = 0; k < N; k++) add(x.data());
  }

  /**
   * \brief Prevent copy constructor, i.e. array2d<int> a2 = a1;
   */
  array2d<T>& operator=(const array2d<T>&) = delete;

  /**
   * \brief Prevent copy assignment operator, i.e. array2d<int> a2 = a1;
   */
  array2d(const array2d<T>&) = delete;

  /**
   * \brief Sets the array layout
   *
   * \param[in] layout - either Layout_Jagged or Layout_Rectangular
   */
  void set_layout(LayoutCategory layout) { layout_ = layout; }

  /**
   * \brief Sets the array stride
   *
   * \param[in] stride - number of elements along second dimension
   *                     stride > 0 for rectangular arrays
   *                     stride < 0 for jagged arrays
   */
  void set_stride(int stride) { stride_ = stride; }

  /**
   * \brief Returns the number of elements along the first dimension.
   *
   * \return number of elements (such as number of vertices, or triangles).
   */
  size_t n() const {
    if (layout_ == Layout_Rectangular) {
      ASSERT(stride_ > 0);
      return data_.size() / stride_;
    }
    return length_.size();
  }

  /**
   * \brief Read access of a value at component (k,j)
   *
   * \param[in] k: 0 <= k < n()
   * \param[in] j: 0 <= j < length(j) or 0 <= j < stride
   *
   * \return const reference to component (k,j)
   */
  const T& operator()(size_t k, uint32_t j) const {
    ASSERT(k < n()) << fmt::format(
        "attempt to access element ({}, {}), but n = {}", k, j, n());
    if (layout_ == Layout_Rectangular) return data_[k * stride_ + j];
    ASSERT(layout_ == Layout_Jagged);
    ASSERT(j < length_[k]) << fmt::format(
        "attempt to access item {} but length({}) = {}", j, k, length_[k]);
    return data_[first_[k] + j];
  }

  /**
   * \brief Read/write access of a value at component (k,j)
   *
   * \param[in] k: 0 <= k < n()
   * \param[in] j: 0 <= j < length(j) or 0 <= j < stride
   *
   * \return reference to component (k,j)
   */
  T& operator()(size_t k, uint32_t j) {
    ASSERT(k < n());
    if (layout_ == Layout_Rectangular) return data_[k * stride_ + j];
    ASSERT(layout_ == Layout_Jagged);
    ASSERT(j < length_[k]);
    return data_[first_[k] + j];
  }

  /**
   * \brief Read access of pointer to element k
   *
   * \param[in] k: 0 <= k < n()
   *
   * \return const pointer to data at element k
   */
  const T* operator[](size_t k) const {
    ASSERT(k < n()) << fmt::format("attempt to access element {}, but n = {}",
                                   k, n());
    if (layout_ == Layout_Rectangular) return &data_[k * stride_];
    ASSERT(layout_ == Layout_Jagged);
    return &data_[first_[k]];
  }

  /**
   * \brief Read/write access of pointer to element k
   *
   * \param[in] k: 0 <= k < n()
   *
   * \return pointer to data at element k
   */
  T* operator[](size_t k) {
    ASSERT(k < n());
    if (layout_ == Layout_Rectangular) return &data_[k * stride_];
    ASSERT(layout_ == Layout_Jagged);
    return &data_[first_[k]];
  }

  /**
   * \brief Adds an element to a rectangular array.
   *
   * \param[in] x - values to add (there must be stride elements)
   */
  template <typename R> void add(const R* x) {
    ASSERT(layout_ == Layout_Rectangular);
    for (int j = 0; j < stride_; j++) data_.push_back(x[j]);
  }

  /**
   * \brief Adds an element to a jagged array.
   *
   * \param[in] x - values to add (there must be stride elements)
   * \param[in] n - number of values to add
   */
  template <typename R> void add(const R* x, int m) {
    if (layout_ == Layout_Rectangular) {
      add(x);
      return;
    }
    ASSERT(layout_ == Layout_Jagged);
    first_.push_back(data_.size());
    for (int j = 0; j < m; j++) data_.push_back(x[j]);
    length_.push_back(m);
  }

  /**
   * \brief Returns the number of items at element k
   *
   * \param[in] k: 0 <= k < n()
   *
   * \return number of items at k
   */
  int length(size_t k) const {
    ASSERT(k < n());
    if (layout_ == Layout_Rectangular) return stride_;
    return length_[k];
  }

  /**
   * \brief Removes all entries stored for element k
   *
   * \param[in]: k0, the element to remove
   *
   */
  void remove(size_t k0) {
    ASSERT(k0 < n()) << fmt::format(
        "attempt to remove element {}, but n() = {}", k0, n());
    int start, end;
    if (layout_ == Layout_Rectangular) {
      start = k0 * stride_;
      end = (k0 + 1) * stride_;
    } else {
      start = first_[k0];
      end = first_[k0] + length_[k0];
      for (size_t k = k0 + 1; k < n(); k++) {
        first_[k] -= length_[k0];  // every element after k0 is now shifted
      }
      first_.erase(first_.begin() + k0);
      length_.erase(length_.begin() + k0);
    }
    data_.erase(data_.begin() + start, data_.begin() + end);
  }

  /**
   * \brief Returns all the data (read-only)
   *
   * \return data vector
   */
  const std::vector<T>& data() const { return data_; }
  const std::vector<index_t>& first() const { return first_; }
  const std::vector<uint32_t>& length() const { return length_; }
  void set_data(const std::vector<T>& d) { data_ = d; }
  void set_first(const std::vector<index_t>& f) { first_ = f; }
  void set_length(const std::vector<uint32_t>& l) { length_ = l; }

  void reserve(int64_t n, int max_item = 10) {
    if (layout_ == Layout_Rectangular) {
      data_.reserve(n * stride_);
    } else {
      first_.reserve(n);
      length_.reserve(n);
      data_.reserve(10 * n);
    }
  }

  /**
   * \brief Clears all the data in this array.
   *
   * \param[in] reset_stride - whether stride should be set to 0.
   *            e.g. with 3d vertices, we still want the stride to be 3.
   */
  void clear(bool reset_stride = false) {
    // do not reset layout because this should not change
    if (reset_stride) stride_ = 0;
    data_.clear();
    first_.clear();
    length_.clear();
  }

  /**
   * \brief Copies the contents of this array to another one.
   *
   * \param[in] dst - destination of the copy.
   */
  void copy(array2d& dst) const {
    ASSERT(dst.stride() == stride_);
    dst.set_data(data_);
    dst.set_first(first_);
    dst.set_length(length_);
  }

  /**
   * \brief Prints out the array as a table.
   *
   * \param[in] label (optional) - prefix for array entries.
   */
  void print(const std::string& label = std::string()) const {
    for (size_t k = 0; k < n(); k++) {
      std::cout << (label.empty() ? "entry" : label) << "[";
      std::cout << k << "] = (";
      int m = length(k);
      for (int j = 0; j < length(k); j++) {
        std::cout << (*this)(k, j);
        if (j + 1 == m)
          std::cout << ")" << std::endl;
        else
          std::cout << ", ";
      }
    }
  }

  /**
   * \brief Returns the stride - only derived classes have access.
   *
   * \return stride (either -1 or > 0).
   */
  int stride() const { return stride_; }

  LayoutCategory layout() const { return layout_; }

 private:
  LayoutCategory layout_;         // either Layout_Jagged or Layout_Rectangular
  int stride_;                    // number of entries along second dimension
  std::vector<T> data_;           // the 1d data array
  std::vector<index_t> first_;    // index of first entry in jagged array
  std::vector<uint32_t> length_;  // number of entries for each element in
                                  // jagged array
};

}  // namespace vortex
