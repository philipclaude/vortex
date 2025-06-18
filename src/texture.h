//
//  vortex: Voronoi mesher and fluid simulator for the Earth's oceans and
//  atmosphere.
//
//  Copyright 2023 - 2025 Philip Claude Caplan
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace vortex {

enum class TextureFormat : uint8_t {
  kRGB,
  kRGBA,
  kGrayscale,
};

static const std::map<TextureFormat, int> kFormat2Channels = {
    {TextureFormat::kRGB, 3},
    {TextureFormat::kRGBA, 4},
    {TextureFormat::kGrayscale, 1}};

struct TextureOptions {
  TextureFormat format{TextureFormat::kRGB};
  bool flipy{true};  // should the y-component be flipped?
};

/// @brief Represents an image that can be sampled to determine properties of a
/// point on a surface. For example, determining the color or height at a point
/// on a surface mesh.
class Texture {
 public:
  /// @brief Reads an image file
  /// @param filename Path to the image
  /// @param options (see above)
  Texture(const std::string& filename, TextureOptions options);

  /// @brief Sample the texture at a point in [0, 1] x [0, 1]
  /// @param s horizontal parameter coordinate (in [0, 1])
  /// @param t vertical parameter coordinate (in [0, 1])
  /// @param value sampled value
  void sample(double s, double t, double* value) const;

  /// @brief Returns the number of pixels in the horizontal direction.
  int width() const { return width_; }

  /// @brief Returns the number of pixels in the vertical direction.
  int height() const { return height_; }

  /// @brief Returns the number of channels (components) for each pixel: 1 for
  /// grayscale and 3 for RGB.
  int channels() const { return channels_; }

  /// @brief Caps the grayscale pixel values to be either min or max.
  /// Any value below threshold will be set to min, and above will be max.
  void make_binary(double threshold, double min, double max);

  /// @brief Force grayscale values to be in the [min, max] range.
  /// @param reverse option to flip the pixel values (255 - value).
  void limit(double min, double max, bool reverse = false);

  /// @brief Forces the pixel values at u = 0 and u = 1 to be the same.
  void make_periodic();

  /// @brief Applies n_iter iterations of Laplacian smoothing to the image.
  void smooth(int n_iter);

  /// @brief Writes the image to filename.
  void write(const std::string& filename) const;

  /// @brief Returns a pointer to the pixel data.
  const auto* data() const { return data_.data(); }

 private:
  void read(const std::string& filename, bool flipy);
  uint8_t channels_;
  int width_;
  int height_;
  std::vector<unsigned char> data_;
};

}  // namespace vortex