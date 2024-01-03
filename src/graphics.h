#pragma once

#include <memory>

namespace wings {
class RenderingServer;
}

namespace vortex {

class Mesh;
class MeshScene;
class ShaderLibrary;
class Viewer {
 public:
  Viewer(const Mesh& mesh, int port, bool earth);
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