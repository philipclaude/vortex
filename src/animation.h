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

#include <string>
#include <vector>

#include "graphics.h"
#include "io.h"
#include "mesh.h"
#include "voronoi.h"

namespace vortex {

class ShaderLibrary;

class ParticleAnimation {
 public:
  ParticleAnimation(const ParticleAnimationParameters& params);

  // void setup();
  // void load();
  // void load_frame(size_t frame, size_t buffer_offset);
  // void reload();
  // void render(const ShaderLibrary& shaders, int time_step);

  void render(const ShaderLibrary& shaders, int time_step);
  void setup();
  void load();
  void reload();
  void load_frame(size_t frame, size_t buffer_offset);

  const auto& shader_names() const { return shader_names_; }
  const auto& params() const { return params_; }

 private:
  // size_t current_time_;
  int current_time_ = 0;
  int buffer_idx_ = 0;
  int current_buffer_ = 0;
  int loading_buffer_ = 1;
  size_t buffer_size_ = std::thread::hardware_concurrency();
  size_t loaded_frames_ = 0;
  ParticleAnimationParameters params_;
  std::vector<GLuint> point_buffer_;
  std::vector<GLuint> velocity_buffer_;
  std::vector<GLuint> density_buffer_;
  std::vector<GLuint> pressure_buffer_;
  std::vector<std::string> shader_names_;
  std::mutex load_mutex_;

  void multiThreadLoad(size_t start, size_t end, size_t buffer_offset);
};

}  // namespace vortex