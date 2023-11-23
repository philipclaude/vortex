#include "graphics.h"

#include <fstream>
#include <set>

#include "mesh.h"
#include "shaders/colormaps.h"
#include "wings.h"
#include "wings/util/glm.h"
#include "wings/util/shader.h"

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

namespace vortex {

enum BufferOption {
  JAGGED_TEXTURE = 1,
  RECTANGULAR_TEXTURE = 2,
  NOT_TEXTURED = 0
};

enum TextureIndex {
  POINT_TEXTURE = 0,
  NORMAL_TEXTURE = 1,
  COLORMAP_TEXTURE = 3,
  INDEX_TEXTURE = 4,
  FIRST_TEXTURE = 5,
  LENGTH_TEXTURE = 6,
  FIELD_TEXTURE = 7,
  IMAGE_TEXTURE = 8
};

struct AABB {
  AABB() {}
  wings::vec3f min{1e20f, 1e20f, 1e20f};
  wings::vec3f max{-1e20f, -1e20f, -1e20f};
};

class GLPrimitive {
 public:
  template <typename T>
  GLPrimitive(const std::string& title, const Topology<T>& topology,
              const std::string& name, enum BufferOption option)
      : title_(title), name_(name), option_(option) {
    write(topology);
  }

  template <typename T>
  void write(const Topology<T>& topology) {
    if (option_ == JAGGED_TEXTURE) {
      std::vector<GLuint> indices(topology.data().begin(),
                                  topology.data().end());
      GL_CALL(glGenBuffers(1, &index_buffer_));
      GL_CALL(glBindBuffer(GL_TEXTURE_BUFFER, index_buffer_));
      GL_CALL(glBufferData(GL_TEXTURE_BUFFER, sizeof(GLuint) * indices.size(),
                           indices.data(), GL_STATIC_DRAW));
      GL_CALL(glBindBuffer(GL_TEXTURE_BUFFER, 0));

      std::vector<GLuint> first(topology.first().begin(),
                                topology.first().end());
      GL_CALL(glGenBuffers(1, &first_buffer_));
      GL_CALL(glBindBuffer(GL_TEXTURE_BUFFER, first_buffer_));
      GL_CALL(glBufferData(GL_TEXTURE_BUFFER, sizeof(GLuint) * first.size(),
                           first.data(), GL_STATIC_DRAW));
      GL_CALL(glBindBuffer(GL_TEXTURE_BUFFER, 0));

      std::vector<GLuint> length(topology.length().begin(),
                                 topology.length().end());
      GL_CALL(glGenBuffers(1, &length_buffer_));
      GL_CALL(glBindBuffer(GL_TEXTURE_BUFFER, length_buffer_));
      GL_CALL(glBufferData(GL_TEXTURE_BUFFER, sizeof(GLuint) * length.size(),
                           length.data(), GL_STATIC_DRAW));
      GL_CALL(glBindBuffer(GL_TEXTURE_BUFFER, 0));

      n_draw_ = first.size();
    } else if (option_ == RECTANGULAR_TEXTURE) {
      // get all the group indices of the elements in the topology
      std::set<int> groups;
      for (int k = 0; k < topology.n(); k++) groups.insert(topology.group(k));
      std::vector<GLuint> order(topology.n());
      int count = 0;
      int igroup = 0;
      n_draw_group_.resize(groups.size());
      for (int group : groups) {
        for (int k = 0; k < topology.n(); k++) {
          if (topology.group(k) != group) continue;
          order[count++] = k;
        }
        n_draw_group_[igroup++] = count;
      }
      std::vector<GLuint> indices(topology.data().begin(),
                                  topology.data().end());
      GL_CALL(glGenBuffers(1, &index_buffer_));
      GL_CALL(glBindBuffer(GL_TEXTURE_BUFFER, index_buffer_));
      GL_CALL(glBufferData(GL_TEXTURE_BUFFER, sizeof(GLuint) * indices.size(),
                           indices.data(), GL_STATIC_DRAW));
      GL_CALL(glBindBuffer(GL_TEXTURE_BUFFER, 0));
      n_draw_ = topology.n();
    } else {
      std::vector<GLuint> indices(topology.data().begin(),
                                  topology.data().end());
      GL_CALL(glGenBuffers(1, &index_buffer_));
      GL_CALL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_));
      GL_CALL(glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                           sizeof(GLuint) * indices.size(), indices.data(),
                           GL_STATIC_DRAW));
      GL_CALL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0));
      n_draw_ = indices.size();
    }
    LOG << fmt::format("wrote {} {}", n_draw_, title_);

    // generate a buffer for the field
    GL_CALL(glGenBuffers(1, &field_buffer_));
  }
  void write_field(const Field& field, int rank) {
    int n_basis = -1;
    const array2d<coord_t>* data = nullptr;
    if (title_ == "Lines") data = &(field.lines()), n_basis = field.lines().m();
    if (title_ == "Triangles")
      data = &(field.triangles()), n_basis = field.triangles().m();
    if (title_ == "Quads") data = &(field.quads()), n_basis = field.quads().m();
    if (title_ == "Polygons")
      data = &(field.polygons()), n_basis = field.polygons().m();
    ASSERT(data != nullptr);

    // extract the appropriate rank and buffer the data
    std::vector<GLfloat> u(data->n() * n_basis, 0.0f);
    for (int k = 0; k < data->n(); k++) {
      for (int j = 0; j < n_basis; j++) {
        u[k * n_basis + j] = (*data)[k][rank * n_basis + j];
      }
    }
    if (u.size() == 0) return;

    umin_ = *std::min_element(u.begin(), u.end());
    umax_ = *std::max_element(u.begin(), u.end());

    // bind the data to the buffer
    GL_CALL(glBindBuffer(GL_TEXTURE_BUFFER, field_buffer_));
    GL_CALL(glBufferData(GL_TEXTURE_BUFFER, sizeof(GLfloat) * u.size(),
                         u.data(), GL_STATIC_DRAW));
    GL_CALL(glBindBuffer(GL_TEXTURE_BUFFER, 0));
  }

  void draw(const wings::ShaderProgram& shader, GLuint point_buffer,
            GLuint field_texture, GLuint index_texture, GLuint first_texture,
            GLuint length_texture) const {
    shader.use();

    // bind the field buffer to the field texture
    GL_CALL(glActiveTexture(GL_TEXTURE0 + FIELD_TEXTURE));
    GL_CALL(glBindTexture(GL_TEXTURE_BUFFER, field_texture));
    GL_CALL(glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, field_buffer_));
    shader.set_uniform("field", int(FIELD_TEXTURE));

    if (option_ == JAGGED_TEXTURE) {
      GL_CALL(glActiveTexture(GL_TEXTURE0 + INDEX_TEXTURE));
      GL_CALL(glBindTexture(GL_TEXTURE_BUFFER, index_texture));
      GL_CALL(glTexBuffer(GL_TEXTURE_BUFFER, GL_R32UI, index_buffer_));
      shader.set_uniform("index", int(INDEX_TEXTURE));

      GL_CALL(glActiveTexture(GL_TEXTURE0 + FIRST_TEXTURE));
      GL_CALL(glBindTexture(GL_TEXTURE_BUFFER, first_texture));
      GL_CALL(glTexBuffer(GL_TEXTURE_BUFFER, GL_R32UI, first_buffer_));
      shader.set_uniform("first", int(FIRST_TEXTURE));

      GL_CALL(glActiveTexture(GL_TEXTURE0 + LENGTH_TEXTURE));
      GL_CALL(glBindTexture(GL_TEXTURE_BUFFER, length_texture));
      GL_CALL(glTexBuffer(GL_TEXTURE_BUFFER, GL_R32UI, length_buffer_));
      shader.set_uniform("count", int(LENGTH_TEXTURE));

      GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, point_buffer));
      GL_CALL(glDrawArrays(GL_POINTS, 0, n_draw_));
    } else if (option_ == RECTANGULAR_TEXTURE) {
      GL_CALL(glActiveTexture(GL_TEXTURE0 + INDEX_TEXTURE));
      GL_CALL(glBindTexture(GL_TEXTURE_BUFFER, index_texture));
      GL_CALL(glTexBuffer(GL_TEXTURE_BUFFER, GL_R32UI, index_buffer_));
      shader.set_uniform("index", int(INDEX_TEXTURE));

      GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, point_buffer));
      GL_CALL(glBindBuffer(GL_TEXTURE_BUFFER, index_buffer_));
      GL_CALL(glDrawArrays(GL_POINTS, 0, n_draw_));
    } else {
      GL_CALL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_));
      GL_CALL(glDrawElements(GL_TRIANGLES, n_draw_, GL_UNSIGNED_INT, 0));
    }
  }
  const std::string& name() const { return name_; }
  const std::string& title() const { return title_; }

  float umin() const { return umin_; }
  float umax() const { return umax_; }

 private:
  index_t n_draw_;
  std::vector<index_t> n_draw_group_;
  std::string title_;
  std::string name_;
  BufferOption option_;

  GLuint index_buffer_;
  GLuint first_buffer_;
  GLuint length_buffer_;
  GLuint field_buffer_;

  float umin_;
  float umax_;

  std::set<int> groups_;
};

class ShaderLibrary {
 public:
  void create();

  void add(const std::string& name, const std::string& prefix,
           bool with_geometry, bool with_tessellation,
           const std::vector<std::string>& macros = {}) {
    shaders_.insert({name, wings::ShaderProgram()});
    std::string base = std::string(VORTEX_SOURCE_DIR) + "/shaders/";
    shaders_[name].set_source(base, name, with_geometry, with_tessellation,
                              macros);
  }

  const wings::ShaderProgram& operator[](const std::string& name) const {
    ASSERT(shaders_.find(name) != shaders_.end())
        << "could not find shader " << name;
    return shaders_.at(name);
  }

 private:
  std::map<std::string, wings::ShaderProgram> shaders_;
};

void ShaderLibrary::create() {
  std::string version = "#version " + std::to_string(VORTEX_GL_VERSION_MAJOR) +
                        std::to_string(VORTEX_GL_VERSION_MINOR) + "0";
  add("points", "points", false, false, {version, "#define WITH_GS 0"});
  add("nodes", "points", true, false, {version, "#define WITH_GS 1"});
  add("edges-q1-p0", "edges", true, false, {version, "#define ORDER 0"});
  add("triangles-q1-p0", "triangles", true, false,
      {version, "#define ORDER 0"});
  add("triangles-q1-p1", "triangles", true, false,
      {version, "#define ORDER 1"});
  add("triangles-q1-pt", "triangles", true, false,
      {version, "#define ORDER -1"});
  add("quads-q1-p0", "quads", true, false, {version, "#define ORDER 0"});
  add("quads-q1-p1", "quads", true, false, {version, "#define ORDER 1"});
  add("quads-q1-pt", "quads", true, false, {version, "#define ORDER -1"});
  add("polygons-q1-p0", "polygons", true, false, {version, "#define ORDER 0"});
}

struct PickableObject {
  template <typename T>
  PickableObject(const Vertices& vertices, const Topology<T>& topology,
                 index_t k, const std::string& name);

  template <typename T>
  void save_points(const Vertices& vertices, const Topology<T>& topology,
                   index_t k);

  double intersection(const wings::vec3f& point, const wings::vec3f& ray,
                      const wings::mat4f& model_matrix) const;
  double intersection(int k, const wings::vec3f& point, const wings::vec3f& ray,
                      const wings::mat4f& model_matrix) const;

  int n_triangles() const { return triangles.size() / 3; }

  std::string name;
  std::vector<wings::vec4f> points;
  std::vector<index_t> triangles;
  std::vector<index_t> nodes;
  index_t index;
};

template <typename T>
void PickableObject::save_points(const Vertices& vertices,
                                 const Topology<T>& topology, index_t k) {
  index = k;
  int dim = (vertices.dim() >= 3) ? 3 : vertices.dim();
  points.resize(topology.length(k));
  nodes.resize(points.size());
  for (index_t j = 0; j < topology.length(k); j++) {
    nodes[j] = topology[k][j];
    for (int d = 0; d < dim; d++) points[j][d] = vertices[nodes[j]][d];
    for (int d = dim; d < 3; d++) points[j][d] = 0.0;  // in case the mesh is 2d
    points[j][3] = 1.0;
  }
}

template <>
PickableObject::PickableObject(const Vertices& vertices,
                               const Topology<Triangle>& topology, index_t k,
                               const std::string& n)
    : name(n) {
  save_points(vertices, topology, k);
  triangles = {0, 1, 2};
}

template <>
PickableObject::PickableObject(const Vertices& vertices,
                               const Topology<Quad>& topology, index_t k,
                               const std::string& n)
    : name(n) {
  save_points(vertices, topology, k);
  triangles = {0, 1, 2, 0, 2, 3};
}

template <>
PickableObject::PickableObject(const Vertices& vertices,
                               const Topology<Polygon>& topology, index_t k,
                               const std::string& n)
    : name(n) {
  save_points(vertices, topology, k);
  triangles.resize(3 * (topology.length(k) - 2));
  int m = 0;
  for (int j = 2; j < topology.length(k); j++) {
    triangles[m++] = 0;
    triangles[m++] = j - 1;
    triangles[m++] = j;
  }
}

double PickableObject::intersection(int k, const wings::vec3f& eye,
                                    const wings::vec3f& ray,
                                    const wings::mat4f& model_matrix) const {
  int i0 = triangles[3 * k];
  int i1 = triangles[3 * k + 1];
  int i2 = triangles[3 * k + 2];

  wings::vec4f q0 = model_matrix * points[i0];
  wings::vec4f q1 = model_matrix * points[i1];
  wings::vec4f q2 = model_matrix * points[i2];

  wings::mat3f A;
  wings::vec3f b;
  for (int d = 0; d < 3; d++) {
    A(d, 0) = q0[d] - q2[d];
    A(d, 1) = q1[d] - q2[d];
    A(d, 2) = -ray[d];
    b[d] = eye[d] - q2[d];
  }
  wings::vec3f c = wings::glm::inverse(A) * b;
  float alpha = c[0], beta = c[1], gamma = 1.0 - alpha - beta;
  float t = c[2];
  if (alpha < 0.0 || alpha > 1.0) return 1e20;
  if (beta < 0.0 || beta > 1.0) return 1e20;
  if (gamma < 0.0 || gamma > 1.0) return 1e20;
  if (t < 0.0) return 1e20;
  return t;
}

double PickableObject::intersection(const wings::vec3f& eye,
                                    const wings::vec3f& ray,
                                    const wings::mat4f& model_matrix) const {
  float tmin = 1e20;
  for (int k = 0; k < n_triangles(); k++) {
    float tk = intersection(k, eye, ray, model_matrix);
    if (tk < tmin) tmin = tk;
  }
  return tmin;
}

class MeshScene : public wings::Scene {
  struct ClientView {
    wings::mat4f model_matrix;
    wings::mat4f view_matrix;
    wings::mat4f projection_matrix;
    wings::mat4f center_translation, inverse_center_translation;
    wings::mat4f translation_matrix;
    wings::vec3f center, eye;
    int width{800};
    int height{600};
    float size{1.0};
    float fov{45};
    double x{0}, y{0};
    GLuint vertex_array;
    std::unordered_map<std::string, bool> active = {
        {"Points", false},     {"Nodes", false},  {"Lines", false},
        {"Triangles", true},   {"Quads", true},   {"Polygons", true},
        {"Tetrahedra", false}, {"Prisms", false}, {"Pyramids", false},
        {"Polyhedra", false}};
    int show_wireframe{1};
    float transparency{1.0};
    int lighting{1};
    bool culling{false};
    const PickableObject* picked{nullptr};
    int field_mode{0};
    int field_index{0};
  };

 public:
  MeshScene(const Mesh& mesh) : mesh_(mesh) {
    context_ =
        wings::RenderingContext::create(wings::RenderingContextType::kOpenGL);
    context_->print();
    shaders_.create();
    write(&mesh_.fields());
    build_pickables();
  }

  void build_pickables() {
    pickables_.clear();
    auto add_pickables = [&](const Vertices& vertices, const auto& topology,
                             const std::string& name) {
      for (int k = 0; k < topology.n(); k++) {
        pickables_.emplace_back(vertices, topology, k, name);
      }
    };

    add_pickables(mesh_.vertices(), mesh_.triangles(), "Triangles");
    add_pickables(mesh_.vertices(), mesh_.quads(), "Quads");
    add_pickables(mesh_.vertices(), mesh_.polygons(), "Polygons");
  }

  const PickableObject* pick(float x, float y, const ClientView& view) {
    // calculate the basis for the camera transformation which is needed
    // to calculate the ray direction
    wings::vec3f g = view.center - view.eye;
    wings::vec3f up = {0, 1, 0};
    wings::vec3f w = -1.0f * unit_vector(g);
    wings::vec3f u = unit_vector(cross(up, w));
    wings::vec3f v = cross(w, u);
    wings::mat3f basis;
    for (int d = 0; d < 3; d++) {
      basis(d, 0) = u[d];
      basis(d, 1) = v[d];
      basis(d, 2) = w[d];
    }

    // computes the 3d world coordinates of a pixel
    // can save some computation by not adding eye, but leaving it for now
    auto pixel2world = [&](double u, double v) {
      double d = length(view.center - view.eye);
      double a = double(view.width) / double(view.height);
      double h = 2.0 * d * tan(view.fov / 2.0);
      double w = a * h;

      float pu = -0.5 * w + w * u;
      float pv = -0.5 * h + h * v;
      float pw = -d;
      wings::vec3f q = {pu, pv, pw};

      return basis * q + view.eye;
    };
    wings::vec3f ray = unit_vector(
        pixel2world(x / view.width, /*1.0 - */ y / view.height) - view.eye);

    // find the closest element
    double tmin = 1e20;
    const PickableObject* picked = nullptr;
    for (size_t k = 0; k < pickables_.size(); k++) {
      const PickableObject& object = pickables_[k];
      if (!view.active.at(object.name)) continue;

      double t = object.intersection(view.eye, ray, view.model_matrix);
      if (t < tmin) {
        tmin = t;
        picked = &object;
      }
    }
    if (picked != nullptr) {
      // TODO more general element printing
      // LOG << fmt::format("picked element {}", picked->index);
    }
    return picked;
  }

  void write(const FieldLibrary* fields = nullptr) {
    context_->make_context_current();

    // generate textures
    GL_CALL(glGenTextures(1, &point_texture_));
    GL_CALL(glGenTextures(1, &index_texture_));
    GL_CALL(glGenTextures(1, &first_texture_));
    GL_CALL(glGenTextures(1, &length_texture_));
    GL_CALL(glGenTextures(1, &field_texture_));
    GL_CALL(glGenTextures(1, &image_texture_));

    // write the point data
    std::vector<GLfloat> coordinates(3 * mesh_.vertices().n());
    for (int i = 0; i < mesh_.vertices().n(); i++)
      for (int d = 0; d < 3; d++) {
        float x = mesh_.vertices()[i][d];
        coordinates[3 * i + d] = x;
        if (x < aabb_.min[d]) aabb_.min[d] = x;
        if (x > aabb_.max[d]) aabb_.max[d] = x;
      }
    GL_CALL(glGenBuffers(1, &point_buffer_));
    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, point_buffer_));
    GL_CALL(glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * coordinates.size(),
                         coordinates.data(), GL_STATIC_DRAW));
    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, 0));

    coordinates.clear();
    n_nodes_ = 0;
    for (int i = 0; i < mesh_.vertices().n(); i++) {
      // TODO get corners
      n_nodes_++;
      for (int d = 0; d < 3; d++) coordinates.push_back(mesh_.vertices()[i][d]);
    }
    if (n_nodes_ > 0) {
      GL_CALL(glGenBuffers(1, &node_buffer_));
      GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, node_buffer_));
      GL_CALL(glBufferData(GL_ARRAY_BUFFER,
                           sizeof(GLfloat) * coordinates.size(),
                           coordinates.data(), GL_STATIC_DRAW));
      GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, 0));
    }

    // generate a texture to hold the mesh coordinates
    GL_CALL(glActiveTexture(GL_TEXTURE0 + POINT_TEXTURE));
    GL_CALL(glBindTexture(GL_TEXTURE_BUFFER, point_texture_));
    GL_CALL(glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, point_buffer_));

    // write the primitives (edges, triangles, quads, polygons)
    primitives_.reserve(4);

    if (mesh_.triangles().n() > 0) {
      primitives_.emplace_back("Triangles", mesh_.triangles(), "triangles-q1",
                               RECTANGULAR_TEXTURE);
    }

    if (mesh_.quads().n() > 0) {
      primitives_.emplace_back("Quads", mesh_.quads(), "quads-q1",
                               RECTANGULAR_TEXTURE);
    }

    if (mesh_.polygons().n() > 0) {
      primitives_.emplace_back("Polygons", mesh_.polygons(), "polygons-q1",
                               JAGGED_TEXTURE);
    }

    if (mesh_.lines().n() > 0) {
      primitives_.emplace_back("Lines", mesh_.lines(), "edges-q1",
                               RECTANGULAR_TEXTURE);
    }

    if (fields != nullptr) {
      // pick one of the fields and activate it
      fields_ = fields;
      const auto& f = fields->fields().begin()->second;
      change_field(f, 0);  // activate rank 0

      field_names_.clear();
      for (const auto& [name, f] : fields->fields()) {
        for (int i = 0; i < f.ranks(); i++) {
          field_names_.push_back({name, i});
        }
      }
    }

    draw_order_.resize(primitives_.size());
    for (size_t k = 0; k < primitives_.size(); k++) draw_order_[k] = k;

    // generate initial buffer & texture for the colormap
    GL_CALL(glGenBuffers(1, &colormap_buffer_));
    GL_CALL(glGenTextures(1, &colormap_texture_));
    change_colormap("giraffe");
  }

  void change_colormap(const std::string& name) {
    index_t n_color = 256 * 3;
    const float* colormap = nullptr;
    if (name == "giraffe") colormap = colormaps::color_giraffe;
    if (name == "viridis") colormap = colormaps::color_viridis;
    if (name == "parula") colormap = colormaps::color_parula;
    if (name == "blue-white-red") colormap = colormaps::color_bwr;
    if (name == "blue-green-red") colormap = colormaps::color_bgr;
    if (name == "jet") colormap = colormaps::color_jet;
    if (name == "hot") colormap = colormaps::color_hot;
    if (name == "hsv") colormap = colormaps::color_hsv;
    ASSERT(colormap != nullptr);

    GL_CALL(glBindBuffer(GL_TEXTURE_BUFFER, colormap_buffer_));
    GL_CALL(glBufferData(GL_TEXTURE_BUFFER, sizeof(GLfloat) * n_color, colormap,
                         GL_STATIC_DRAW));
    GL_CALL(glActiveTexture(GL_TEXTURE0 + COLORMAP_TEXTURE));
    GL_CALL(glBindTexture(GL_TEXTURE_BUFFER, colormap_texture_));
    GL_CALL(glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, colormap_buffer_));
  }

  void change_field(const Field& field, int rank) {
    // change the field in all the primitives
    for (size_t k = 0; k < primitives_.size(); k++) {
      primitives_[k].write_field(field, rank);
    }
  }

  void center_view(ClientView& view, const wings::vec3f& point) {
    wings::vec4f p(point.data(), 3);
    p[3] = 1.0;
    wings::vec3f q = (view.model_matrix * p).xyz();
    float len = length(view.center - view.eye);
    wings::vec3f dir = unit_vector(view.center - view.eye);
    view.center = q;
    view.eye = view.center - len * dir;

    const wings::vec3f up = {0, 1, 0};
    view.view_matrix = wings::glm::lookat(view.eye, view.center, up);

    view.center_translation.eye();
    view.center_translation(0, 3) = view.center[0];
    view.center_translation(1, 3) = view.center[1];
    view.center_translation(2, 3) = view.center[2];
  }

  void center_view(ClientView& view) {
    if (!view.picked) return;
    wings::vec4f point;
    for (const auto& p : view.picked->points) point = point + p;
    point = (1.0f / view.picked->points.size()) * point;
    center_view(view, point.xyz());
  }

  bool render(const wings::ClientInput& input, int client_idx,
              std::string* msg) {
    ClientView& view = view_[client_idx];
    bool updated = false;
    switch (input.type) {
      case wings::InputType::MouseMotion: {
        if (input.dragging) {
          if (!input.modifier) {
            double dx = (view.x - input.x) / view.width;
            double dy = (view.y - input.y) / view.height;
            wings::mat4f R = view.center_translation * view.translation_matrix *
                             wings::glm::rotation(dx, dy) *
                             wings::glm::inverse(view.translation_matrix *
                                                 view.center_translation);
            view.model_matrix = R * view.model_matrix;
          } else {
            double dx = -(view.x - input.x) / view.width;
            double dy = -(view.y - input.y) / view.height;
            dx *= view.size;
            dy *= view.size;
            wings::mat4f T = wings::glm::translation(dx, dy);
            view.translation_matrix = T * view.translation_matrix;
            view.model_matrix = T * view.model_matrix;
          }
          updated = true;
        }
        view.x = input.x;
        view.y = input.y;
        break;
      }
      case wings::InputType::DoubleClick: {
        view.picked = pick(input.x, input.y, view);
        if (view.picked) {
          std::string info = fmt::format("*picked element {} ({}): (",
                                         view.picked->index, view.picked->name);
          size_t i = 0;
          for (auto p : view.picked->nodes) {
            info += std::to_string(p);
            if (i + 1 < view.picked->nodes.size())
              info += ", ";
            else
              info += ")";
            i++;
          }
          *msg = info;
        }
        updated = true;
        break;
      }
      case wings::InputType::KeyValueInt: {
        updated = true;
        if (input.key == 'Q')
          quality_ = input.ivalue;
        else if (input.key == 'n')
          view.active["Nodes"] = input.ivalue > 0;
        else if (input.key == 'v')
          view.active["Points"] = input.ivalue > 0;
        else if (input.key == 'e')
          view.active["Lines"] = input.ivalue > 0;
        else if (input.key == 't')
          view.active["Triangles"] = input.ivalue > 0;
        else if (input.key == 'T')
          view.active["Tetrahedra"] = input.ivalue > 0;
        else if (input.key == 'p')
          view.active["Polygons"] = input.ivalue > 0;
        else if (input.key == 'y')
          view.active["Prisms"] = input.ivalue > 0;
        else if (input.key == 'Y')
          view.active["Pyramids"] = input.ivalue > 0;
        else if (input.key == 'P')
          view.active["Polyhedra"] = input.ivalue > 0;
        else if (input.key == 'a')
          view.transparency = 0.01 * input.ivalue;
        else if (input.key == 'w')
          view.show_wireframe = input.ivalue > 0;
        else {
          updated = false;
        }
        break;
      }
      case wings::InputType::KeyValueStr: {
        if (input.key == 'C') {  // colormap change
          change_colormap(input.svalue);
          updated = true;
        } else if (input.key == 'p') {
          center_view(view);
          updated = true;
        } else if (input.key == 'f') {
          if (fields_) {
            if (view.field_mode == 0) {
              view.field_mode = 1;
              view.field_index = 0;
            } else if (size_t(view.field_index + 1) == field_names_.size()) {
              view.field_mode = 0;
              view.field_index = 0;
            } else {
              view.field_mode = 1;
              view.field_index++;  // = (view.field_index + 1) % view.field_;
            }
            auto field_info = field_names_[view.field_index];
            const auto& field = fields_->fields().at(field_info.first);
            int rank = field_info.second;
            change_field(field, rank);
            if (view.field_mode == 1)
              *msg =
                  fmt::format("*plotting \"{}\" ({})", field_info.first, rank);
            else
              *msg = "*no field";
            updated = true;
          }
        }
        break;
      }
      case wings::InputType::Scroll: {
        wings::vec3f direction = view.eye - view.center;
        view.eye = view.center + direction * 1.0f * input.fvalue;
        wings::vec3f up{0, 1, 0};
        view.view_matrix = wings::glm::lookat(view.eye, view.center, up);
        updated = true;
        break;
      }
      default:
        break;
    }
    if (!updated) return false;

    // write shader uniforms
    GL_CALL(glViewport(0, 0, view.width, view.height));
    GL_CALL(glClearColor(1.0f, 1.0f, 1.0f, 1.0f));
    GL_CALL(glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT));
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_POLYGON_SMOOTH);
    glDepthMask(GL_TRUE);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1, 1);
    glEnable(GL_BLEND);
    // glEnable(GL_STENCIL_TEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_DEPTH_CLAMP);
    float alpha = view.transparency;
    int u_lighting = view.lighting;
    // glDepthRange(0.1, 0.9);
    if (alpha < 1.0) {
      glDisable(GL_CULL_FACE);
      // glDepthMask(GL_FALSE);
      // glDepthFunc(GL_LEQUAL);
      u_lighting = 0;
    }

    // compute the matrices
    wings::mat4f model_view_matrix = view.view_matrix * view.model_matrix;
    wings::mat4f normal_matrix =
        wings::glm::inverse(wings::glm::transpose(model_view_matrix));
    wings::mat4f mvp_matrix = view.projection_matrix * model_view_matrix;

    if (view.active["Nodes"] && n_nodes_ > 0) {
      // bind which attributes we want to draw
      GL_CALL(glBindVertexArray(view.vertex_array));
      GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, node_buffer_));
      GL_CALL(glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0));
      GL_CALL(glEnableVertexAttribArray(0));
      GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, 0));

      const auto& shader = shaders_["points"];
      shader.use();
      shader.set_uniform("u_ModelViewProjectionMatrix", mvp_matrix);
      shader.set_uniform("u_length", view.size);
      GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, node_buffer_));
      GL_CALL(glDrawArrays(GL_POINTS, 0, n_nodes_));
      GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, 0));
    }

    // bind which attributes we want to draw
    GL_CALL(glBindVertexArray(view.vertex_array));
    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, point_buffer_));
    GL_CALL(glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0));
    GL_CALL(glEnableVertexAttribArray(0));
    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, 0));

    if (view.active["Points"]) {
      const auto& shader = shaders_["points"];
      shader.use();
      shader.set_uniform("u_ModelViewProjectionMatrix", mvp_matrix);
      GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, point_buffer_));
      GL_CALL(glDrawArrays(GL_POINTS, 0, mesh_.vertices().n()));
      GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, 0));
    }

    // generate a texture to hold the mesh coordinates
    GL_CALL(glActiveTexture(GL_TEXTURE0 + POINT_TEXTURE));
    GL_CALL(glBindTexture(GL_TEXTURE_BUFFER, point_texture_));
    GL_CALL(glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, point_buffer_));

    for (size_t k = 0; k < primitives_.size(); k++) {
      const GLPrimitive& primitive = primitives_[draw_order_[k]];
      if (!view.active[primitive.title()]) continue;

      int order = 0;  // parameters["field order"];
      const auto& shader = select_shader(primitive, order);
      shader.use();

      // bind the image texture (if any)
      glActiveTexture(GL_TEXTURE0 + IMAGE_TEXTURE);
      GL_CALL(glBindTexture(GL_TEXTURE_2D, image_texture_));
      shader.set_uniform("image", int(IMAGE_TEXTURE));

      // bind the desired colormap
      glActiveTexture(GL_TEXTURE0 + COLORMAP_TEXTURE);
      GL_CALL(glBindTexture(GL_TEXTURE_BUFFER, colormap_texture_));
      shader.set_uniform("colormap", int(COLORMAP_TEXTURE));

      // set the uniforms for the shader
      shader.set_uniform("u_edges", view.show_wireframe);
      shader.set_uniform("u_lighting", u_lighting);
      shader.set_uniform("u_alpha", alpha);

      shader.set_uniform("u_ModelViewProjectionMatrix", mvp_matrix);
      shader.set_uniform("u_ModelViewMatrix", model_view_matrix);
      shader.set_uniform("u_NormalMatrix", normal_matrix);
      // shader.set_uniform("u_ViewportSize", screen_size_);

      shader.set_uniform("u_umin", primitive.umin());
      shader.set_uniform("u_umax", primitive.umax());
      shader.set_uniform("u_field_mode", view.field_mode);
      shader.set_uniform("u_picking", int(-1));

      if (view.picked) {
        if (view.picked->name == primitive.title()) {
          shader.set_uniform("u_picking", int(view.picked->index));
        }
      }

      primitive.draw(shader, point_buffer_, field_texture_, index_texture_,
                     first_texture_, length_texture_);
    }

    // save the pixels in the wings::Scene
    GLsizei channels = 3;
    GLsizei stride = channels * view.width;
    stride += (stride % 4) ? (4 - stride % 4) : 0;
    pixels_.resize(stride * view.height);
    GL_CALL(glPixelStorei(GL_PACK_ALIGNMENT, 4));
    GL_CALL(glReadPixels(0, 0, view.width, view.height, GL_RGB,
                         GL_UNSIGNED_BYTE, pixels_.data()));
    glFinish();

    return true;
  }

  void onconnect() {
    // set up the view
    view_.emplace_back();
    ClientView& view = view_.back();
    view.eye = {0, 0, 0};
    view.center = {0, 0, 0};
    wings::vec3f xmin{1e20f, 1e20f, 1e20f}, xmax{-1e20f, -1e20f, -1e20f};
    for (int i = 0; i < mesh_.vertices().n(); i++) {
      for (int d = 0; d < 3; d++) {
        float x = mesh_.vertices()(i, d);
        view.center[d] += x;
        if (x > xmax[d]) xmax[d] = x;
        if (x < xmin[d]) xmin[d] = x;
      }
    }
    view.center = view.center / (1.0f * mesh_.vertices().n());
    view.center_translation.eye();
    view.inverse_center_translation.eye();
    for (int d = 0; d < 3; d++) {
      view.center_translation(d, 3) = view.center[d];
      view.inverse_center_translation(d, 3) = -view.center[d];
    }

    wings::vec3f scale = aabb_.max - aabb_.min;
    view.size = std::max(std::max(scale[0], scale[1]), scale[2]);
    float d = scale[2] / 2.0 + scale[0] / (2.0 * tan(view.fov / 2.0));

    wings::vec3f dir{0, 0, 1};
    view.eye = view.center + 1.05f * d * dir;
    view.center = view.eye - 2.1f * d * dir;

    view.model_matrix.eye();
    wings::vec3f up{0, 1, 0};
    view.view_matrix = wings::glm::lookat(view.eye, view.center, up);
    view.projection_matrix = wings::glm::perspective(
        view.fov, float(view.width) / float(view.height), 0.1f, 10000.0f);
    view.translation_matrix.eye();

    // vertex arrays are not shared between OpenGL contexts in different threads
    // (buffers & textures are though). Each client runs in a separate thread
    // so we need to create a vertex array upon each client connection.
    GL_CALL(glGenVertexArrays(1, &view.vertex_array));
    view.field_mode = 0;
  }

  const wings::ShaderProgram& select_shader(const GLPrimitive& primitive,
                                            int order) const {
    std::string suffix = (order < 0) ? "t" : std::to_string(order);
    return shaders_[primitive.name() + "-p" + suffix];
  }

 private:
  const Mesh& mesh_;
  const FieldLibrary* fields_{nullptr};
  std::vector<std::pair<std::string, int>> field_names_;

  std::vector<ClientView> view_;
  AABB aabb_;

  GLuint point_buffer_;
  GLuint point_texture_;
  GLuint node_buffer_;
  GLuint n_nodes_{0};

  GLuint index_texture_;
  GLuint first_texture_;
  GLuint length_texture_;

  GLuint image_texture_;
  GLuint field_texture_;
  GLuint colormap_texture_;
  GLuint colormap_buffer_;

  std::vector<GLPrimitive> primitives_;
  std::vector<int> draw_order_;
  ShaderLibrary shaders_;
  std::vector<PickableObject> pickables_;
};

Viewer::Viewer(const Mesh& mesh, int port) {
  scene_ = std::make_unique<MeshScene>(mesh);
  renderer_ = std::make_unique<wings::RenderingServer>(*scene_, port);
}

Viewer::~Viewer() {}

}  // namespace vortex