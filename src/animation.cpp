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
#include "animation.h"

#include <fmt/format.h>

#include "log.h"
#include "wings/util/glm.h"
#include "wings/util/shader.h"

namespace vortex {

ParticleAnimation::ParticleAnimation(const ParticleAnimationParameters& params)
    : params_(params) {
  shader_names_.push_back("particles");
  // TODO add other shader names if necessary
  // for example, if there are other shaders for pathlines or other effects
}

void ParticleAnimation::render(const ShaderLibrary& shaders, int time_step) {
  const auto& particle_shader = shaders["particles"];
  particle_shader.use();

  if (time_step >= 0) {
    LOG << fmt::format("render particles @ time {}!", time_step);
    current_time_ = time_step;
  } else {
    // a time step < 0 indicates we are not animating, but we should
    // render whatever time step was rendered last, which might be
    // because the animation ended/paused and the user is modifying the view
    LOG << fmt::format("render particles @ time {}!", current_time_);
  }
}

}  // namespace vortex