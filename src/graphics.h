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
#include <memory>
#include <string>
#include <vector>

#define GL_GLEXT_PROTOTYPES
#ifdef __APPLE__
#define __gl_h_
#define GL_DO_NOT_WARN_IF_MULTI_GL_VERSION_HEADERS_INCLUDED
#define GL_SILENCE_DEPRECATION
#include <OpenGL/OpenGL.h>
#include <OpenGL/gl3.h>
#else
#include <GL/gl.h>
#endif

#define VORTEX_GL_VERSION_MAJOR 3
#define VORTEX_GL_VERSION_MINOR 3

#ifndef VORTEX_SOURCE_DIR
#define VORTEX_SOURCE_DIR "./"
#endif

namespace wings {
class RenderingServer;
class ShaderProgram;
}  // namespace wings

namespace vortex {

struct ParticleAnimationParameters {
  std::string points_prefix;
  std::string velocity_prefix;
  std::string density_prefix;
  std::string pressure_prefix;
  std::string file_type;
  std::vector<double> time;
  int num_particles{100};
  int total_frames;
  int hasDensity_{0};

  bool active() const { return !points_prefix.empty(); }
};

class Mesh;
class MeshScene;

class ShaderLibrary {
 public:
  void create();

  void add(const std::string& name, const std::string& prefix,
           const std::vector<std::string>& macros, bool with_gs = false);

  const wings::ShaderProgram& operator[](const std::string& name) const;

 private:
  std::map<std::string, wings::ShaderProgram> shaders_;
};
class Viewer {
 public:
  Viewer(const Mesh& mesh, const ParticleAnimationParameters& particle_params,
         int port);
  ~Viewer();

 private:
  std::unique_ptr<MeshScene> scene_;
  std::unique_ptr<wings::RenderingServer> renderer_;
  std::unique_ptr<ShaderLibrary> shaders_;
};

#define GL_CALL(X)                                                          \
  {                                                                         \
    (X);                                                                    \
    [[maybe_unused]] bool error = false;                                    \
    for (GLenum glerr = glGetError(); glerr != GL_NO_ERROR;                 \
         glerr = glGetError()) {                                            \
      const char* message = "";                                             \
      switch (glerr) {                                                      \
        case GL_INVALID_ENUM:                                               \
          message = "invalid enum";                                         \
          break;                                                            \
        case GL_INVALID_VALUE:                                              \
          message = "invalid value";                                        \
          break;                                                            \
        case GL_INVALID_OPERATION:                                          \
          message = "invalid operation";                                    \
          break;                                                            \
        case GL_INVALID_FRAMEBUFFER_OPERATION:                              \
          message = "invalid framebuffer operation";                        \
          break;                                                            \
        case GL_OUT_OF_MEMORY:                                              \
          message = "out of memory";                                        \
          break;                                                            \
        default:                                                            \
          message = "unknown error";                                        \
      }                                                                     \
      std::cout << "OpenGL error at " << __FILE__ << "(" << __LINE__ << ")" \
                << message << std::endl;                                    \
      error = true;                                                         \
    }                                                                       \
    ASSERT(!error);                                                         \
  }

}  // namespace vortex