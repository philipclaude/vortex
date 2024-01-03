//
//  vortex: Voronoi mesher and fluid simulator for the Earth's oceans and
//  atmosphere.
//
//  Copyright 2023 - 2024 Philip Claude Caplan
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
  bool flipy{false};
};

class Texture {
 public:
  Texture(const std::string& filename, TextureOptions options);

  void sample(double s, double t, double* value) const;
  int width() const { return width_; }
  int height() const { return height_; }
  int channels() const { return channels_; }
  void make_binary(double threshold, double min, double max);
  void limit(double min, double max, bool reverse = false);
  void make_periodic();
  void smooth(int n_iter);
  void write(const std::string& filename) const;
  const auto* data() const { return data_.data(); }

 private:
  void read(const std::string& filename);
  uint8_t channels_;
  bool flipy_{false};
  int width_;
  int height_;
  std::vector<unsigned char> data_;
};

}  // namespace vortex