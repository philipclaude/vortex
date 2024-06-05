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
#include <algorithm>
#include <string>

#include "mesh.h"

namespace vortex {

class Mesh;

namespace meshb {

/// @brief Reads a .meshb file and stores it in a Mesh object.
/// @param filename Path to the .meshb file.
/// @param mesh Destination of the mesh.
void read(const std::string& filename, Mesh& mesh);

/// @brief Writes a .meshb file.
/// @param mesh Mesh object to write.
/// @param filename Path to the output .meshb file.
/// @param twod whether we are writing vertices in 2d.
void write(const Mesh& mesh, const std::string& filename, bool twod = false);

}  // namespace meshb

namespace obj {

/// @brief Reads a .obj file and stores it in a Mesh object.
/// @param filename Path to the .obj file.
/// @param mesh Destination of the mesh.
void read(const std::string& filename, Mesh& mesh);

/// @brief Writes a .obj file.
/// @param mesh Mesh object to write.
/// @param filename Path to the output .obj file.
void write(const Mesh& mesh, const std::string& filename);

}  // namespace obj

/// @brief Retrieves the extension of the file.
/// @param filename
/// @return string containing the characeters after the last '.'
inline std::string get_file_ext(const std::string& filename) {
  std::string::size_type idx;
  idx = filename.rfind('.');  // find the '.' in reverse order
  if (idx != std::string::npos) return filename.substr(idx + 1);
  return "";
}

/// @brief Reads a .meshb or .obj file and stores it in a Mesh object.
/// @param filename Path to the mesh file.
/// @param mesh Destination of the mesh.
inline void read_mesh(const std::string& filename, Mesh& mesh) {
  std::string ext = get_file_ext(filename);
  std::transform(ext.begin(), ext.end(), ext.begin(),
                 [](unsigned char c) { return std::tolower(c); });

  if (ext == "mesh" || ext == "meshb")
    meshb::read(filename, mesh);
  else if (ext == "obj")
    obj::read(filename, mesh);
  else
    NOT_IMPLEMENTED;
}

namespace vtk {

/// @brief Writes a .vtk file.
/// @param vertices Vertices object to write.
/// @param filename Path to the output .vtk file.
void write(const Vertices& vertices, const std::string& filename);

}  // namespace vtk

namespace io {
template <typename T>
void swap_end(T& var) {
  // https://stackoverflow.com/questions/10913666/error-writing-binary-vtk-files
  char* varArray = reinterpret_cast<char*>(&var);
  for (long i = 0; i < static_cast<long>(sizeof(var) / 2); i++)
    std::swap(varArray[sizeof(var) - 1 - i], varArray[i]);
}
}  // namespace io

}  // namespace vortex