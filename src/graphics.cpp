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

enum TextureIndex {
  kPoint = 0,
  kColormap = 1,
  kIndex = 2,
  kVisibility = 3,
  kPrimitive2Cell = 4,
  kField = 5,
  kImage = 6
};

enum RenderingPipeline {
  kElements,
  kArrays,
};

struct AABB {
  AABB() {}
  wings::vec3f min{1e20f, 1e20f, 1e20f};
  wings::vec3f max{-1e20f, -1e20f, -1e20f};
};

class GLPrimitive {
 public:
  template <typename T>
  GLPrimitive(const std::string& title, const Vertices& vertices,
              const Topology<T>& topology, const std::string& name,
              RenderingPipeline pipeline, GLenum type)
      : title_(title), name_(name), pipeline_(pipeline), type_(type) {
    // write the topology and save the number to draw
    write(vertices, topology);

    // generate a buffer for the field
    GL_CALL(glGenBuffers(1, &field_buffer_));
  }

  void write(const std::vector<GLuint>& indices,
             const std::vector<char>& visibility,
             const std::vector<GLuint>& primitive2cell) {
    GL_CALL(glGenBuffers(1, &index_buffer_));
    GL_CALL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_));
    GL_CALL(glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                         sizeof(GLuint) * indices.size(), indices.data(),
                         GL_STATIC_DRAW));
    GL_CALL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0));

    GL_CALL(glGenBuffers(1, &visibility_buffer_));
    GL_CALL(glBindBuffer(GL_TEXTURE_BUFFER, visibility_buffer_));
    GL_CALL(glBufferData(GL_TEXTURE_BUFFER, sizeof(char) * visibility.size(),
                         visibility.data(), GL_STATIC_DRAW));
    GL_CALL(glBindBuffer(GL_TEXTURE_BUFFER, 0));

    if (!primitive2cell.empty()) {
      GL_CALL(glGenBuffers(1, &primitive2cell_buffer_));
      GL_CALL(glBindBuffer(GL_TEXTURE_BUFFER, primitive2cell_buffer_));
      GL_CALL(glBufferData(GL_TEXTURE_BUFFER,
                           sizeof(GLuint) * primitive2cell.size(),
                           primitive2cell.data(), GL_STATIC_DRAW));
      GL_CALL(glBindBuffer(GL_TEXTURE_BUFFER, 0));
    }
  }

  template <typename T>
  void write(const Vertices& vertices, const Topology<T>& topology) {
    // write the elements
    std::vector<GLuint> indices(topology.data().begin(), topology.data().end());
    std::vector<char> visibility(topology.data().size(), char(0));
    write(indices, visibility, {});
    n_draw_ = topology.data().size();
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
    for (size_t k = 0; k < data->n(); k++) {
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
            GLuint field_texture, GLuint index_texture,
            GLuint visibility_texture, GLuint primitive2cell_texture) const {
    shader.use();

    // bind the field buffer to the field texture
    GL_CALL(glActiveTexture(GL_TEXTURE0 + kField));
    GL_CALL(glBindTexture(GL_TEXTURE_BUFFER, field_texture));
    GL_CALL(glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, field_buffer_));
    shader.set_uniform("field", int(kField));

    GL_CALL(glActiveTexture(GL_TEXTURE0 + kIndex));
    GL_CALL(glBindTexture(GL_TEXTURE_BUFFER, index_texture));
    GL_CALL(glTexBuffer(GL_TEXTURE_BUFFER, GL_R32UI, index_buffer_));
    shader.set_uniform("index", int(kIndex));

    if (type_ == GL_TRIANGLES) {
      GL_CALL(glActiveTexture(GL_TEXTURE0 + kVisibility));
      GL_CALL(glBindTexture(GL_TEXTURE_BUFFER, visibility_texture));
      GL_CALL(glTexBuffer(GL_TEXTURE_BUFFER, GL_R8, visibility_buffer_));
      shader.set_uniform("visibility", int(kVisibility));
    }

    if (split_primitives_) {
      GL_CALL(glActiveTexture(GL_TEXTURE0 + kPrimitive2Cell));
      GL_CALL(glBindTexture(GL_TEXTURE_BUFFER, primitive2cell_texture));
      GL_CALL(glTexBuffer(GL_TEXTURE_BUFFER, GL_R32UI, primitive2cell_buffer_));
      shader.set_uniform("primitive2cell", int(kPrimitive2Cell));
    }

    if (pipeline_ == kArrays) {
      GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, point_buffer));
      GL_CALL(glDrawArrays(type_, 0, n_draw_));
    } else if (pipeline_ == kElements) {
      GL_CALL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_));
      GL_CALL(glDrawElements(type_, n_draw_, GL_UNSIGNED_INT, 0));
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

  GLuint index_buffer_;
  GLuint field_buffer_;
  GLuint visibility_buffer_;
  GLuint primitive2cell_buffer_;
  bool split_primitives_{false};

  float umin_;
  float umax_;

  std::set<int> groups_;
  RenderingPipeline pipeline_;
  GLenum type_;
};

template <>
void GLPrimitive::write(const Vertices& vertices,
                        const Topology<Polygon>& polygons) {
  // reserve enough space for the triangles
  std::vector<GLuint> indices, primitive2cell;
  std::vector<char> visibility;
  indices.reserve(polygons.data().size() * 3);
  visibility.reserve(indices.size());
  primitive2cell.reserve(polygons.data().size());

  // TODO(philip): for each polygon, first check if basic
  // tessellation can be used, if not:
  // project all points to the tangent plane of average point
  // and perform ear clipping in the tangent plane
  for (size_t k = 0; k < polygons.n(); k++) {
    auto* pk = polygons[k];
    auto nk = polygons.length(k);

    for (int j = 2; j < nk; j++) {
      indices.push_back(pk[0]);
      indices.push_back(pk[j - 1]);
      indices.push_back(pk[j]);

      int e0 = 0;  // always visible
      int e1 = (j + 1 == nk) ? 0 : 1;
      int e2 = (j == 2) ? 0 : 1;
      visibility.push_back(e0);
      visibility.push_back(e1);
      visibility.push_back(e2);

      primitive2cell.push_back(k);
    }
  }
  write(indices, visibility, primitive2cell);

  n_draw_ = indices.size();
  split_primitives_ = true;
}

template <>
void GLPrimitive::write(const Vertices& vertices, const Topology<Quad>& quads) {
  // reserve enough space for the triangles
  std::vector<GLuint> indices(quads.n() * 6), primitive2cell(2 * quads.n());
  std::vector<char> visibility(indices.size());

  size_t t = 0;
  for (size_t k = 0; k < quads.n(); k++) {
    auto* qk = quads[k];

    indices[3 * t + 0] = qk[0];
    indices[3 * t + 1] = qk[1];
    indices[3 * t + 2] = qk[2];
    visibility[3 * t + 0] = 0;
    visibility[3 * t + 1] = 1;
    visibility[3 * t + 2] = 0;
    primitive2cell[t] = k;
    t++;

    indices[3 * t + 0] = qk[0];
    indices[3 * t + 1] = qk[2];
    indices[3 * t + 2] = qk[3];
    visibility[3 * t + 0] = 0;
    visibility[3 * t + 1] = 0;
    visibility[3 * t + 2] = 1;
    primitive2cell[t] = k;
    t++;
  }
  write(indices, visibility, primitive2cell);

  n_draw_ = indices.size();
  split_primitives_ = true;
}

class ShaderLibrary {
 public:
  void create();

  void add(const std::string& name, const std::string& prefix,
           const std::vector<std::string>& macros, bool with_gs = false) {
    shaders_.insert({name, wings::ShaderProgram()});
    std::string base = std::string(VORTEX_SOURCE_DIR) + "/shaders/";
    shaders_[name].set_source(base, prefix, with_gs, false, macros);
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
  add("points", "points", {version});
  add("edges-q1-p0", "edges", {version, "#define ORDER 0"}, true);
  add("triangles-q1-p0", "triangles", {version, "#define ORDER 0"});
  add("triangles-q1-p1", "triangles", {version, "#define ORDER 1"});
  add("quads-q1-p0", "triangles",
      {version, "#define ORDER 0", "#define SPLIT_PRIMITIVES"});
  add("quads-q1-p1", "triangles",
      {version, "#define ORDER 1", "#define SPLIT_PRIMITIVES"});
  add("polygons-q1-p0", "triangles",
      {version, "#define ORDER 0", "#define SPLIT_PRIMITIVES"});
  add("polygons-q1-p1", "triangles",
      {version, "#define ORDER 1", "#define SPLIT_PRIMITIVES"});
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
  for (int j = 0; j < topology.length(k); j++) {
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
    float size{1.0};
    float fov{45};
    double x{0}, y{0};
    GLuint vertex_array;
    std::unordered_map<std::string, bool> active = {
        {"Points", false},   {"Nodes", false}, {"Lines", false},
        {"Triangles", true}, {"Quads", true},  {"Polygons", true}};
    int show_wireframe{1};
    float transparency{1.0};
    int lighting{1};
    bool culling{false};
    const PickableObject* picked{nullptr};
    int field_mode{0};
    int field_index{0};
    wings::glCanvas canvas{800, 600};
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
      for (size_t k = 0; k < topology.n(); k++) {
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
      double a = double(view.canvas.width) / double(view.canvas.height);
      double h = 2.0 * d * tan(view.fov / 2.0);
      double w = a * h;

      float pu = -0.5 * w + w * u;
      float pv = -0.5 * h + h * v;
      float pw = -d;
      wings::vec3f q = {pu, pv, pw};

      return basis * q + view.eye;
    };
    wings::vec3f ray =
        unit_vector(pixel2world(x / view.canvas.width,
                                /*1.0 - */ y / view.canvas.height) -
                    view.eye);

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
    GL_CALL(glGenTextures(1, &visibility_texture_));
    GL_CALL(glGenTextures(1, &primitive2cell_texture_));
    GL_CALL(glGenTextures(1, &field_texture_));
    GL_CALL(glGenTextures(1, &image_texture_));

    // write the point data
    std::vector<GLfloat> coordinates(3 * mesh_.vertices().n());
    for (size_t i = 0; i < mesh_.vertices().n(); i++)
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
    for (size_t i = 0; i < mesh_.vertices().n(); i++) {
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
    GL_CALL(glActiveTexture(GL_TEXTURE0 + kPoint));
    GL_CALL(glBindTexture(GL_TEXTURE_BUFFER, point_texture_));
    GL_CALL(glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, point_buffer_));

    // write the primitives (edges, triangles, quads, polygons)
    primitives_.reserve(4);

    if (mesh_.triangles().n() > 0) {
      primitives_.emplace_back("Triangles", mesh_.vertices(), mesh_.triangles(),
                               "triangles-q1", kElements, GL_TRIANGLES);
    }

    if (mesh_.quads().n() > 0) {
      primitives_.emplace_back("Quads", mesh_.vertices(), mesh_.quads(),
                               "quads-q1", kElements, GL_TRIANGLES);
    }

    if (mesh_.polygons().n() > 0) {
      primitives_.emplace_back("Polygons", mesh_.vertices(), mesh_.polygons(),
                               "polygons-q1", kElements, GL_TRIANGLES);
    }

    if (mesh_.lines().n() > 0) {
      primitives_.emplace_back("Lines", mesh_.vertices(), mesh_.lines(),
                               "edges-q1", kArrays, GL_POINTS);
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
    if (name == "blue-white-red") colormap = colormaps::color_bwr;
    if (name == "blue-green-red") colormap = colormaps::color_bgr;
    if (name == "jet") colormap = colormaps::color_jet;
    if (name == "hot") colormap = colormaps::color_hot;
    if (name == "hsv") colormap = colormaps::color_hsv;
    ASSERT(colormap != nullptr);

    GL_CALL(glBindBuffer(GL_TEXTURE_BUFFER, colormap_buffer_));
    GL_CALL(glBufferData(GL_TEXTURE_BUFFER, sizeof(GLfloat) * n_color, colormap,
                         GL_STATIC_DRAW));
    GL_CALL(glActiveTexture(GL_TEXTURE0 + kColormap));
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
            double dx = (view.x - input.x) / view.canvas.width;
            double dy = (view.y - input.y) / view.canvas.height;
            wings::mat4f R = view.center_translation * view.translation_matrix *
                             wings::glm::rotation(dx, dy) *
                             wings::glm::inverse(view.translation_matrix *
                                                 view.center_translation);
            view.model_matrix = R * view.model_matrix;
          } else {
            double dx = -(view.x - input.x) / view.canvas.width;
            double dy = -(view.y - input.y) / view.canvas.height;
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
        else if (input.key == 'p')
          view.active["Polygons"] = input.ivalue > 0;
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
    GL_CALL(glViewport(0, 0, view.canvas.width, view.canvas.height));
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
    GL_CALL(glActiveTexture(GL_TEXTURE0 + kPoint));
    GL_CALL(glBindTexture(GL_TEXTURE_BUFFER, point_texture_));
    GL_CALL(glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, point_buffer_));

    for (size_t k = 0; k < primitives_.size(); k++) {
      const GLPrimitive& primitive = primitives_[draw_order_[k]];
      if (!view.active[primitive.title()]) continue;

      int order = 0;
      const auto& shader = select_shader(primitive, order);
      shader.use();

      // bind the image texture (if any)
      glActiveTexture(GL_TEXTURE0 + kImage);
      GL_CALL(glBindTexture(GL_TEXTURE_2D, image_texture_));
      shader.set_uniform("image", int(kImage));

      // bind the desired colormap
      glActiveTexture(GL_TEXTURE0 + kColormap);
      GL_CALL(glBindTexture(GL_TEXTURE_BUFFER, colormap_texture_));
      shader.set_uniform("colormap", int(kColormap));

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
                     visibility_texture_, primitive2cell_texture_);
    }

    // save the pixels in the wings::Scene
    GLsizei channels = 3;
    GLsizei stride = channels * view.canvas.width;
    stride += (stride % 4) ? (4 - stride % 4) : 0;
    pixels_.resize(stride * view.canvas.height);
    GL_CALL(glPixelStorei(GL_PACK_ALIGNMENT, 4));
    GL_CALL(glReadPixels(0, 0, view.canvas.width, view.canvas.height, GL_RGB,
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
    for (size_t i = 0; i < mesh_.vertices().n(); i++) {
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
        view.fov, float(view.canvas.width) / float(view.canvas.height), 0.1f,
        10000.0f);
    view.translation_matrix.eye();

    // vertex arrays are not shared between OpenGL contexts in different
    // threads (buffers & textures are though). Each client runs in a separate
    // thread so we need to create a vertex array upon each client connection.
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
  GLuint visibility_texture_;
  GLuint primitive2cell_texture_;

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