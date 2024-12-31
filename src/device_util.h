#pragma once
#include <cstdlib>

#include "log.h"

#define USE_CUDA 0

#if USE_CUDA
#define SHARED __shared__
#else
#define SHARED
#endif

#define MIN2(a, b) (a < b) ? a : b

namespace vortex {

template <typename T, size_t _NUM_THREADS, size_t _CAPACITY>
class thread_memory_pool {
 public:
  static constexpr size_t CAPACITY = _CAPACITY;
  static constexpr size_t NUM_THREADS = _NUM_THREADS;
  thread_memory_pool() {}

  T* block(size_t tid) { return data_ + CAPACITY * tid; }

 private:
  SHARED T data_[NUM_THREADS * CAPACITY];
};

template <typename T>
class device_vector {
 public:
  device_vector(size_t capacity) : resizable_(true) {
#if USE_CUDA
    cudaMalloc(&data_, sizeof(T) * capacity);
#else
    data_ = new T[capacity];
#endif
    size_ = 0;
    capacity_ = capacity;
  }

  device_vector(T* buffer, size_t capacity)
      : size_(0), capacity_(capacity), data_(buffer), resizable_(false) {}

  const T& operator[](size_t k) const { return data_[k]; }
  T& operator[](size_t k) { return data_[k]; }

  size_t size() const { return size_; }

  size_t capacity() const { return capacity_; }

  void push_back(const T& value) {
    if (size_ == capacity_) {
      ASSERT(resizable_) << capacity_;
      allocate(2 * capacity_);
    }
    data_[size_++] = value;
  }

  T& emplace_back() {
    if (size_ == capacity_) {
      ASSERT(resizable_);
      allocate(2 * capacity_);
    }
    ++size_;
    return data_[size_ - 1];
  }

  T*& data() { return data_; }

  T& back() { return data_[size_ - 1]; }
  const T& back() const { return data_[size_ - 1]; }
  const T& pop_back() { return data_[size_--]; }
  void clear() { size_ = 0; }
  void clear(const T& value) {
    for (size_t k = 0; k < capacity_; k++) data_[k] = value;
    size_ = 0;
  }
  bool empty() const { return size_ == 0; }

  void resize(size_t size) {
    if (size >= capacity_) {
      ASSERT(resizable_);
      allocate(size);
    }
    size_ = size;
  }

  void set_size(size_t size) {
    ASSERT(size_ <= capacity_);
    size_ = size;
  }

 private:
  void allocate(size_t capacity) {
    ASSERT(resizable_);
    T* new_data = nullptr;
#if USE_CUDA
    cudaMalloc(&new_data, sizeof(T) * capacity);
#else
    new_data = (T*)malloc(sizeof(T) * capacity);
#endif
    size_t limit = MIN2(capacity, capacity_);
    for (size_t k = 0; k < limit; k++) {
      new_data[k] = data_[k];
    }
#if USE_CUDA
    cudaFree(data_);
#else
    free(data_);
#endif
    data_ = new_data;
    capacity_ = capacity;
  }

  size_t size_{0};
  size_t capacity_{0};
  T* data_{nullptr};
  bool resizable_{false};
};

template <typename T>
size_t hash(const T& key) {
  return key;
}

template <typename T>
class device_hash_set {
 public:
  device_hash_set() {
    capacity_ = 16;
    size_ = 0;
#if USE_CUDA
    cudaMalloc(&data_, sizeof(T) * capacity_);
    cudaMalloc(&taken_, sizeof(bool) * capacity_);
#else
    data_ = (T*)malloc(sizeof(T) * capacity_);
    taken_ = (bool*)malloc(sizeof(bool) * capacity_);
#endif
    for (size_t k = 0; k < capacity_; k++) taken_[k] = false;
  }

  void insert(const T& key) {
    insert(key, data_, taken_, capacity_);
    if (2 * size_ > capacity_) {
      allocate(2 * capacity_);
    }
  }

  void insert(const T& key, T* dst, bool* tkn, const size_t m) {
    int64_t index = get_index(key, dst, tkn, m);
    tkn[index] = true;
    if (dst[index] == key) return;  // key already added
    dst[index] = key;
    ++size_;
  }

  bool contains(const T& key) {
    const size_t capacity_minus_one = capacity_ - 1;
    size_t index = hash(key) & capacity_minus_one;
    if (!taken_[index]) return false;
    while (taken_[index]) {
      if (data_[index] == key) return true;
      index = (index + 1) & capacity_minus_one;
    }
    return false;
  }

  void clear() {
    size_ = 0;
    for (size_t k = 0; k < capacity_; k++) {
      taken_[k] = false;
      data_[k] = std::numeric_limits<T>::max();
    }
  }

  void display(bool full = false) {
#if USE_CUDA == 0
    std::cout << "[";
    for (size_t k = 0; k < capacity_; k++) {
      if (!taken_[k] && !full) continue;
      if (taken_[k])
        std::cout << data_[k];
      else
        std::cout << "*";
      if (k + 1 < capacity_) std::cout << ", ";
    }
    std::cout << "]";
#endif
  }

 private:
  int64_t get_index(const T& key, T* dst, bool* tkn, const size_t m) {
    const size_t capacity_minus_one = m - 1;
    size_t index = hash(key) & capacity_minus_one;
    while (tkn[index] && dst[index] != key) {
      index = (index + 1) & capacity_minus_one;
    }
    return index;
  }

  void allocate(size_t capacity) {
    T* new_data = nullptr;
    bool* new_taken = nullptr;
#if USE_CUDA
    cudaMalloc(&new_data, sizeof(T) * capacity);
    cudaMalloc(&new_taken, sizeof(bool) * capacity);
#else
    new_data = (T*)malloc(sizeof(T) * capacity);
    new_taken = (bool*)malloc(sizeof(bool) * capacity);
#endif
    for (int k = 0; k < capacity; k++) {
      new_taken[k] = false;
      new_data[k] = std::numeric_limits<T>::max();
    }
    size_t limit = MIN2(capacity, capacity_);
    for (size_t k = 0; k < limit; k++) {
      if (!taken_[k]) continue;
      insert(data_[k], new_data, new_taken, capacity);
    }

#if USE_CUDA
    cudaFree(data_);
    cudaFree(taken_);
#else
    free(data_);
    free(taken_);
#endif
    data_ = new_data;
    taken_ = new_taken;
    capacity_ = capacity;
  }

  size_t size_;
  size_t capacity_;
  T* data_;
  bool* taken_;
};

}  // namespace vortex

#undef MIN2
