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
#include "particles.h"
#include "wings/util/glm.h"
#include "wings/util/shader.h"

namespace vortex {

ParticleAnimation::ParticleAnimation(const ParticleAnimationParameters& params)
    : params_(params) {
  shader_names_.push_back("particles");
  // TODO add other shader names if necessary
  // for example, if there are other shaders for pathlines or other effects
  std::cout << "ParticleAnimation constructed with " << params_.num_particles
            << " particles.\n";

  //   // Setup the buffers
  //   std::cout << "Calling setup...\n";
  //   setup();

  //   std::cout << "Calling load...\n";
  //   // Load the initial data
  //   load();
}

void ParticleAnimation::render(const ShaderLibrary& shaders, int time_step) {
  const auto& particle_shader = shaders["particles"];
  particle_shader.use();
  particle_shader.set_uniform("u_hasDensity", params_.hasDensity_);

  if (time_step >= 0) {
    // LOG << fmt::format("render particles @ time {}!", time_step);
    current_time_ = time_step;
  } else {
    // a time step < 0 indicates we are not animating, but we should
    // render whatever time step was rendered last, which might be
    // because the animation ended/paused and the user is modifying the view
    current_time_ = 0;
    // LOG << fmt::format("render particles @ time {}!", current_time_);
  }

#if 0
  std::vector<GLfloat> points(params_.num_particles * 3);
  for (auto& x : points) x = double(rand()) / double(RAND_MAX);
  GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, point_buffer_[current_time_]));
  GL_CALL(glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(GLfloat),
                       points.data(), GL_STATIC_DRAW));
#endif

  // Bind buffers and render particles for the current time step
  GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, point_buffer_[current_time_]));
  GL_CALL(glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr));
  GL_CALL(glEnableVertexAttribArray(0));

  // Assuming other attributes like velocity, density, pressure need to be bound
  // Bind density buffer
  GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, density_buffer_[current_time_]));
  GL_CALL(glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, nullptr));
  GL_CALL(glEnableVertexAttribArray(1));

  // Draw particles
  GL_CALL(glDrawArrays(GL_POINTS, 0, params_.num_particles));

  GL_CALL(glDisableVertexAttribArray(0));
  GL_CALL(glDisableVertexAttribArray(1));
  GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, 0));
}

void ParticleAnimation::setup() {
  std::cout << "Entering ParticleAnimation setup..." << std::endl;
  // Ensure num_particles is set properly
  if (params_.num_particles == 0) {
    // Handle error: No particles to setup
    std::cerr << "Error: num_particles is zero in ParticleAnimation setup.\n";
    return;
  }
  std::cout << "Setting up ParticleAnimation with " << params_.num_particles
            << " particles.\n";

  // Resize buffers to accommodate multiple frames
  point_buffer_.resize(params_.time.size());
  velocity_buffer_.resize(params_.time.size());
  density_buffer_.resize(params_.time.size());
  pressure_buffer_.resize(params_.time.size());

  // Generate buffers for each frame
  GL_CALL(glGenBuffers(params_.time.size(), point_buffer_.data()));
  GL_CALL(glGenBuffers(params_.time.size(), velocity_buffer_.data()));
  GL_CALL(glGenBuffers(params_.time.size(), density_buffer_.data()));
  GL_CALL(glGenBuffers(params_.time.size(), pressure_buffer_.data()));

  // Allocate memory for each buffer
  for (size_t i = 0; i < params_.time.size(); ++i) {
    std::vector<GLfloat> tmp(params_.num_particles * sizeof(float) * 3, 0);
    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, point_buffer_[i]));
    GL_CALL(glBufferData(GL_ARRAY_BUFFER,
                         params_.num_particles * sizeof(float) * 3, tmp.data(),
                         GL_STATIC_DRAW));

    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, velocity_buffer_[i]));
    GL_CALL(glBufferData(GL_ARRAY_BUFFER,
                         params_.num_particles * sizeof(float) * 3, nullptr,
                         GL_STATIC_DRAW));

    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, density_buffer_[i]));
    GL_CALL(glBufferData(GL_ARRAY_BUFFER, params_.num_particles * sizeof(float),
                         tmp.data(), GL_STATIC_DRAW));

    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, pressure_buffer_[i]));
    GL_CALL(glBufferData(GL_ARRAY_BUFFER, params_.num_particles * sizeof(float),
                         nullptr, GL_STATIC_DRAW));
  }

  GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, 0));
}

void ParticleAnimation::load() {
  std::cout << "Entering ParticleAnimation load..." << std::endl;
  std::cout << params_.file_type;
  // Load data into buffers for each frame
  for (size_t t = 0; t < params_.time.size(); ++t) {
    std::string points_file = fmt::format(
        "{}/particles{}.{}", params_.points_prefix, t, params_.file_type);
    std::cout << points_file;
    Mesh mesh(3);
    std::vector<float> densities;  // To store the densities
    if (params_.file_type == "obj") {
      obj::read(points_file, mesh);
    } else if (params_.file_type == "mesh" || params_.file_type == "meshb") {
      meshb::read(points_file, mesh);
    } else if (params_.file_type == "sol" || params_.file_type == "solb") {
      meshb::read(points_file, mesh,
                  densities);  // Use the read method to read the .meshb file
    }

    const auto& vertices = mesh.vertices();

    std::vector<float> points(vertices.n() * 3);
    for (size_t i = 0; i < vertices.n(); ++i) {
      points[i * 3] = static_cast<float>(vertices[i][0]);
      points[i * 3 + 1] = static_cast<float>(vertices[i][1]);
      points[i * 3 + 2] = static_cast<float>(vertices[i][2]);
    }

    // Placeholder for velocities, densities, pressures
    std::vector<float> velocities(params_.num_particles * 3, 0.0f);
    std::vector<float> densities_out(params_.num_particles, 0.0f);
    params_.hasDensity_ = !densities.empty() ? 1 : 0;
    if (params_.hasDensity_) {
      // Store the densities read from the file
      for (size_t i = 0; i < vertices.n(); ++i) {
        densities_out[i] = densities[i];
      }
    }

    std::vector<float> pressures(params_.num_particles, 0.0f);

    // Upload data to GPU
    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, point_buffer_[t]));
    // GL_CALL(glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(float),
    //                      points.data(), GL_STATIC_DRAW));
    GL_CALL(glBufferSubData(GL_ARRAY_BUFFER, 0, points.size() * sizeof(float),
                            points.data()));

    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, velocity_buffer_[t]));
    GL_CALL(glBufferSubData(GL_ARRAY_BUFFER, 0,
                            velocities.size() * sizeof(float),
                            velocities.data()));

    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, density_buffer_[t]));
    GL_CALL(glBufferSubData(GL_ARRAY_BUFFER, 0,
                            densities_out.size() * sizeof(float),
                            densities_out.data()));

    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, pressure_buffer_[t]));
    GL_CALL(glBufferSubData(GL_ARRAY_BUFFER, 0,
                            pressures.size() * sizeof(float),
                            pressures.data()));
  }

  GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, 0));
}

}  // namespace vortex

/**

#include <fmt/format.h>

#include "animation.h"
#include "log.h"
#include "particles.h"
#include "wings/util/glm.h"
#include "wings/util/shader.h"

namespace vortex {

ParticleAnimation::ParticleAnimation(const ParticleAnimationParameters& params)
    : params_(params), current_buffer_(0), loading_buffer_(1) {
  shader_names_.push_back("particles");
  std::cout << "ParticleAnimation constructed with " << params_.num_particles
            << " particles.\n";
}

void ParticleAnimation::setup() {
  std::cout << "Entering ParticleAnimation setup..." << std::endl;
  if (params_.num_particles == 0) {
    std::cerr << "Error: num_particles is zero in ParticleAnimation setup.\n";
    return;
  }
  std::cout << "Setting up ParticleAnimation with " << params_.num_particles
            << " particles.\n";

  // Resize buffers to accommodate double buffering
  point_buffer_.resize(2 * buffer_size_);
  velocity_buffer_.resize(2 * buffer_size_);
  density_buffer_.resize(2 * buffer_size_);
  pressure_buffer_.resize(2 * buffer_size_);

  // Generate buffers for both sets
  GL_CALL(glGenBuffers(2 * buffer_size_, point_buffer_.data()));
  GL_CALL(glGenBuffers(2 * buffer_size_, velocity_buffer_.data()));
  GL_CALL(glGenBuffers(2 * buffer_size_, density_buffer_.data()));
  GL_CALL(glGenBuffers(2 * buffer_size_, pressure_buffer_.data()));

  // Allocate memory for each buffer
  for (size_t i = 0; i < 2 * buffer_size_; ++i) {
    std::vector<GLfloat> tmp(params_.num_particles * sizeof(float) * 3, 0);
    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, point_buffer_[i]));
    GL_CALL(glBufferData(GL_ARRAY_BUFFER,
                         params_.num_particles * sizeof(float) * 3, tmp.data(),
                         GL_STATIC_DRAW));

    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, velocity_buffer_[i]));
    GL_CALL(glBufferData(GL_ARRAY_BUFFER,
                         params_.num_particles * sizeof(float) * 3, nullptr,
                         GL_STATIC_DRAW));

    std::vector<GLfloat> tmp1(params_.num_particles * sizeof(float), 0);
    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, density_buffer_[i]));
    GL_CALL(glBufferData(GL_ARRAY_BUFFER, params_.num_particles * sizeof(float),
                         tmp1.data(), GL_STATIC_DRAW));

    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, pressure_buffer_[i]));
    GL_CALL(glBufferData(GL_ARRAY_BUFFER, params_.num_particles * sizeof(float),
                         nullptr, GL_STATIC_DRAW));
  }

  GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, 0));
}

void ParticleAnimation::load_frame(size_t frame, size_t buffer_offset) {
  if (frame >= params_.time.size()) {
    std::cerr << "Error: Frame index out of bounds: " << frame << std::endl;
    return;
  }

  std::string points_file = fmt::format(
      "{}/particles{}.{}", params_.points_prefix, frame, params_.file_type);
  Mesh mesh(3);
  std::vector<float> densities;

  if (params_.file_type == "obj") {
    obj::read(points_file, mesh);
  } else if (params_.file_type == "mesh" || params_.file_type == "meshb") {
    meshb::read(points_file, mesh);
  } else if (params_.file_type == "sol" || params_.file_type == "solb") {
    meshb::read(points_file, mesh, densities);
  }

  const auto& vertices = mesh.vertices();
  std::vector<float> points(vertices.n() * 3);
  for (size_t i = 0; i < vertices.n(); ++i) {
    points[i * 3] = static_cast<float>(vertices[i][0]);
    points[i * 3 + 1] = static_cast<float>(vertices[i][1]);
    points[i * 3 + 2] = static_cast<float>(vertices[i][2]);
  }

  std::vector<float> velocities(params_.num_particles * 3, 0.0f);
  std::vector<float> densities_out(params_.num_particles, 0.0f);
  params_.hasDensity_ = !densities.empty() ? 1 : 0;
  if (params_.hasDensity_) {
    for (size_t i = 0; i < vertices.n(); ++i) {
      densities_out[i] = densities[i];
    }
  }
  std::vector<float> pressures(params_.num_particles, 0.0f);

  size_t index = buffer_offset + (frame % buffer_size_);

  GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, point_buffer_[index]));
  GL_CALL(glBufferSubData(GL_ARRAY_BUFFER, 0, points.size() * sizeof(float),
                          points.data()));

  GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, velocity_buffer_[index]));
  GL_CALL(glBufferSubData(GL_ARRAY_BUFFER, 0, velocities.size() * sizeof(float),
                          velocities.data()));

  GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, density_buffer_[index]));
  GL_CALL(glBufferSubData(GL_ARRAY_BUFFER, 0,
                          densities_out.size() * sizeof(float),
                          densities_out.data()));

  GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, pressure_buffer_[index]));
  GL_CALL(glBufferSubData(GL_ARRAY_BUFFER, 0, pressures.size() * sizeof(float),
                          pressures.data()));

  GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, 0));
}

void ParticleAnimation::load() {
  std::cout << "Entering ParticleAnimation load..." << std::endl;
  for (size_t t = 0; t < buffer_size_; ++t) {
    load_frame(t, 0);
  }
  loaded_frames_ = buffer_size_;
}

void ParticleAnimation::reload() {
  LOG << "Reloading particle animation data...";

  size_t new_size =
      std::min<size_t>(current_time_ + buffer_size_ + 1, params_.time.size());

  if (new_size == params_.time.size()) {
    LOG << "Reached the last frame, starting from the beginning...";
    current_time_ = -1;
    new_size = buffer_size_;
  }

  LOG << fmt::format("Next frame range: {} - {}", current_time_ + 1, new_size);

  size_t buffer_offset = loading_buffer_ * buffer_size_;

  for (size_t t = current_time_ + 1; t < new_size; ++t) {
    load_frame(t, buffer_offset);
  }

  // Swap the buffers
  std::swap(current_buffer_, loading_buffer_);

  loaded_frames_ = new_size;
  std::cout << "There are " << loaded_frames_ << " frames." << std::endl;
}

void ParticleAnimation::render(const ShaderLibrary& shaders, int time_step) {
  const auto& particle_shader = shaders["particles"];
  particle_shader.use();
  particle_shader.set_uniform("u_hasDensity", params_.hasDensity_);

  if (time_step >= 0) {
    if (time_step >= loaded_frames_) {
      std::cerr << "Error: Attempting to render an unloaded frame: "
                << time_step << std::endl;
      reload();  // Attempt to reload more frames if necessary
      if (time_step >= loaded_frames_) {
        current_time_ = loaded_frames_ - 1;
      } else {
        current_time_ = time_step;
      }
    } else {
      current_time_ = time_step;
    }
  } else {
    current_time_ = 0;
  }

  size_t buffer_offset = current_buffer_ * buffer_size_;
  size_t index = buffer_offset + (current_time_ % buffer_size_);

  GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, point_buffer_[index]));
  GL_CALL(glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr));
  GL_CALL(glEnableVertexAttribArray(0));

  GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, density_buffer_[index]));
  GL_CALL(glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, nullptr));
  GL_CALL(glEnableVertexAttribArray(1));

  GL_CALL(glDrawArrays(GL_POINTS, 0, params_.num_particles));

  GL_CALL(glDisableVertexAttribArray(0));
  GL_CALL(glDisableVertexAttribArray(1));
  GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, 0));
}

}  // namespace vortex
*/