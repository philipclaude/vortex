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

#include <memory>
#include <string>

namespace wings {
class RenderingServer;
}

namespace argparse {
class ArgumentParser;
}

namespace vortex {

class Mesh;
class MeshScene;
class ShaderLibrary;
class Viewer {
 public:
  Viewer(const Mesh& mesh, int port, const std::string view = "");
  void save(const std::string& filename,
            const argparse::ArgumentParser& argparse);
  ~Viewer();

 private:
#if VORTEX_WITH_VIZ
  std::unique_ptr<MeshScene> scene_;
  std::unique_ptr<wings::RenderingServer> renderer_;
  std::unique_ptr<ShaderLibrary> shaders_;
#endif
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