#pragma once

namespace vortex {

/*
 * \brief Represents the type of the result of an operation, by promoting the
 * operand types. For example, for x (float) and y (double), x * y would be a
 * double.
 */
template <typename S, typename T>
class result_of;

template <>
class result_of<float, int> {
 public:
  typedef float type;
};
template <>
class result_of<int, float> {
 public:
  typedef float type;
};
template <>
class result_of<int, double> {
 public:
  typedef double type;
};
template <>
class result_of<double, int> {
 public:
  typedef double type;
};
template <>
class result_of<double, double> {
 public:
  typedef double type;
};
template <>
class result_of<float, float> {
 public:
  typedef float type;
};
template <>
class result_of<float, double> {
 public:
  typedef double type;
};
template <>
class result_of<double, float> {
 public:
  typedef double type;
};

}  // namespace vortex
