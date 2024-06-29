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

namespace vortex {

class ShaderLibrary;

class ParticleAnimation {
 public:
  ParticleAnimation(const ParticleAnimationParameters& params);

  void setup();
  void load();
  void reload();
  void render(const ShaderLibrary& shaders, int time_step);

  const auto& shader_names() const { return shader_names_; }

 private:
  ParticleAnimationParameters params_;
  std::vector<double> time_;
  std::vector<GLuint> point_buffer_;
  std::vector<GLuint> velocity_buffer_;
  std::vector<GLuint> density_buffer_;
  std::vector<GLuint> pressure_buffer_;
  std::vector<std::string> shader_names_;
  int current_time_{-1};
};

}  // namespace vortex