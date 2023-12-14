
#include "io.h"

#include <fmt/format.h>
#include <libmeshb7.h>
#include <shapelib/shapefil.h>
#include <tinyobjloader/tiny_obj_loader.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <unordered_set>

#include "mesh.h"
#include "numerics.h"

namespace vortex {

namespace meshb {
void read_vertices(int64_t fid, int version, Vertices& vertices) {
  int status = 1;
  int dim = vertices.dim();
  double dvalues[4] = {0, 0, 0, 0};
  float fvalues[4] = {0, 0, 0, 0};

  ASSERT(GmfGotoKwd(fid, GmfVertices) > 0);
  int n = GmfStatKwd(fid, GmfVertices);
  vertices.reserve(n);
  LOG << fmt::format("reading {} vertices", n);

  for (int k = 0; k < n; k++) {
    int domain = -1;
    if (dim == 2) {
      if (version == 1)
        status = GmfGetLin(fid, GmfVertices, &fvalues[0], &fvalues[1], &domain);
      else
        status = GmfGetLin(fid, GmfVertices, &dvalues[0], &dvalues[1], &domain);
    } else if (dim == 3) {
      if (version == 1)
        GmfGetLin(fid, GmfVertices, &fvalues[0], &fvalues[1], &fvalues[2],
                  &domain);
      else
        GmfGetLin(fid, GmfVertices, &dvalues[0], &dvalues[1], &dvalues[2],
                  &domain);
    }
    ASSERT(status == 1);

    if (version == 1) {
      for (int d = 0; d < 3; d++) dvalues[d] = double(fvalues[d]);
    }
    vertices.add(dvalues, domain - 1);
  }
}

void read_edges(int64_t fid, Mesh& mesh) {
  int status = 1;
  int n;
  int data[3] = {-1, -1, -1};
  index_t edge[3] = {0, 0, 0};

  if (GmfGotoKwd(fid, GmfEdges) <= 0) return;

  n = GmfStatKwd(fid, GmfEdges);
  mesh.lines().reserve(n);
  LOG << fmt::format("reading {} edges", n);

  for (int k = 0; k < n; k++) {
    status = GmfGetLin(fid, GmfEdges, &data[0], &data[1], &data[2]);
    ASSERT(status == 1);

    for (int j = 0; j < 2; j++) edge[j] = data[j] - 1;

    mesh.lines().add(edge);
    mesh.lines().set_group(k, data[2]);
  }
}

void read_triangles(int64_t fid, Mesh& mesh) {
  int status = 1;
  int n;
  int data[4] = {-1, -1, -1, -1};
  index_t triangle[4] = {0, 0, 0, 0};

  if (GmfGotoKwd(fid, GmfTriangles) <= 0) return;

  n = GmfStatKwd(fid, GmfTriangles);
  mesh.triangles().reserve(n);
  LOG << fmt::format("reading {} triangles", n);

  for (int k = 0; k < n; k++) {
    status =
        GmfGetLin(fid, GmfTriangles, &data[0], &data[1], &data[2], &data[3]);
    ASSERT(status == 1);

    for (int j = 0; j < 3; j++) triangle[j] = data[j] - 1;
    if (data[3] == 0) {
      // triangle[0] = triangle[1] = triangle[2] = 0;
    }
    mesh.triangles().add(triangle);
    mesh.triangles().set_group(k, data[3]);
  }
}

void read_quads(int64_t fid, Mesh& mesh) {
  int status = 1;
  int n;
  int data[5] = {-1, -1, -1, -1, -1};
  index_t quad[5] = {0, 0, 0, 0, 0};

  if (GmfGotoKwd(fid, GmfQuadrilaterals) <= 0) return;

  n = GmfStatKwd(fid, GmfQuadrilaterals);
  mesh.quads().reserve(n);
  LOG << fmt::format("reading {} quads", n);

  for (int k = 0; k < n; k++) {
    status = GmfGetLin(fid, GmfQuadrilaterals, &data[0], &data[1], &data[2],
                       &data[3], &data[4]);
    ASSERT(status == 1);

    for (int j = 0; j < 4; j++) quad[j] = data[j] - 1;

    mesh.quads().add(quad);
    mesh.quads().set_group(k, data[4]);
  }
}

void read_polydata(int64_t fid, GmfKwdCod header_type, GmfKwdCod entity_type,
                   std::function<void(std::vector<long>&, int)> add_entity_fn) {
  // read the polygon/polyhedra headers: (a, b)
  // a: first index of the vertex in the polygon/polyhedron
  // b: material reference (which will be used as the group number)
  // GmfGetLin doesn't seem to support GmfInnerPolygonHeaders, so
  // a block-based method is used similar to the utility in
  // extern/libmeshb/utilities/libmeshb7_helpers.c
  int64_t n_poly = GmfStatKwd(fid, header_type);
  std::vector<std::array<long, 2>> header(n_poly);
  ASSERT(GmfGetBlock(fid, header_type, 1, n_poly, 0, nullptr, nullptr,
                     GmfLongVec, 2, &header.front(), &header.back()));

  // read the polygon vertices or polyhedron faces
  int64_t n_entities = GmfStatKwd(fid, entity_type);
  std::vector<long> entities(n_entities);
  ASSERT(GmfGetBlock(fid, entity_type, 1, n_entities, 0, nullptr, nullptr,
                     GmfLong, &entities.front(), &entities.back()));

  // add the entities to the topology using the add_entity callback
  LOG << "creating entities";
  std::vector<long> entity(32);
  for (int64_t k = 0; k < n_poly; k++) {
    index_t u = header[k][0] - 1;
    index_t v = (k + 1 == n_poly) ? n_entities : (header[k + 1][0] - 1);
    entity.resize(v - u);
    for (index_t m = u; m < v; m++) entity[m - u] = entities[m] - 1;
    add_entity_fn(entity, header[k][1]);
  }
}

void read_polygons(int64_t fid, Mesh& mesh) {
  bool have_polygons = GmfGotoKwd(fid, GmfPolygons) > 0;
  bool have_boundary_polygons = GmfGotoKwd(fid, GmfBoundaryPolygonHeaders) > 0;
  if (!have_polygons && !have_boundary_polygons) return;

  if (have_boundary_polygons) {
    ASSERT(!have_polygons) << "should only have one set of polygons";
  }
  if (have_polygons) NOT_IMPLEMENTED;

  // read the boundary polygon data
  read_polydata(fid, GmfBoundaryPolygonHeaders, GmfBoundaryPolygonVertices,
                [&mesh](std::vector<long>& polygon, int group) {
                  std::reverse(polygon.begin(), polygon.end());
                  int64_t k = mesh.polygons().n();
                  mesh.polygons().add(polygon.data(), polygon.size());
                  if (group >= 0) mesh.polygons().set_group(k, group);
                });
}

void write_polygons(int64_t fid, const Mesh& mesh) {
  if (mesh.polygons().n() == 0) return;

  index_t n_boundary = mesh.polygons().n();
  GmfSetKwd(fid, GmfBoundaryPolygonHeaders, n_boundary);
  index_t m = 1;
  for (int k = 0; k < mesh.polygons().n(); k++) {
    GmfSetLin(fid, GmfBoundaryPolygonHeaders, m, mesh.polygons().group(k));
    m += mesh.polygons().length(k);
  }

  GmfSetKwd(fid, GmfBoundaryPolygonVertices, m - 1);
  for (int k = 0; k < mesh.polygons().n(); k++) {
    for (int j = 0; j < mesh.polygons().length(k); j++)
      GmfSetLin(fid, GmfBoundaryPolygonVertices, mesh.polygons()(k, j) + 1);
  }
}

template <typename T>
void write_simplex(int64_t fid, const Topology<T>& topology) {
  const int nv = T::n_vertices;
  int type;
  if (nv == 1)
    type = GmfCorners;
  else if (nv == 2)
    type = GmfEdges;
  else if (nv == 3)
    type = GmfTriangles;
  else if (nv == 4)
    type = GmfTetrahedra;
  else
    NOT_POSSIBLE;

  int n = topology.n();
  GmfSetKwd(fid, type, n);

  int indices[5];
  for (int k = 0; k < n; k++) {
    for (int j = 0; j < nv; j++) indices[j] = topology(k, j) + 1;
    indices[nv] = topology.group(k) + 1;

    if (nv == 1)
      GmfSetLin(fid, type, indices[0]);
    else if (nv == 2)
      GmfSetLin(fid, type, indices[0], indices[1], indices[2]);
    else if (nv == 3)
      GmfSetLin(fid, type, indices[0], indices[1], indices[2], indices[3]);
    else if (nv == 4)
      GmfSetLin(fid, type, indices[0], indices[1], indices[2], indices[3],
                indices[4]);
    else
      NOT_POSSIBLE;
  }
}

void read(const std::string& filename, Mesh& mesh) {
  // open the file
  int version;
  int dim;
  int64_t fid = GmfOpenMesh(filename.c_str(), GmfRead, &version, &dim);
  ASSERT(fid) << "could not open mesh file " << filename;
  mesh.vertices().set_dim(dim);

  read_edges(fid, mesh);
  read_triangles(fid, mesh);
  read_quads(fid, mesh);
  read_polygons(fid, mesh);
  read_vertices(fid, version, mesh.vertices());

  GmfCloseMesh(fid);
}

void write(const Mesh& mesh, const std::string& filename, bool twod) {
  int dim = mesh.vertices().dim();
  if (twod) dim = 2;

  int64_t fid = GmfOpenMesh(filename.c_str(), GmfWrite, GmfDouble, dim);
  ASSERT(fid);

  GmfSetKwd(fid, GmfVertices, mesh.vertices().n());
  for (int k = 0; k < mesh.vertices().n(); k++) {
    int ref = mesh.vertices().group(k) + 1;
    if (dim == 2)
      GmfSetLin(fid, GmfVertices, mesh.vertices()[k][0], mesh.vertices()[k][1],
                ref);
    else if (dim == 3)
      GmfSetLin(fid, GmfVertices, mesh.vertices()[k][0], mesh.vertices()[k][1],
                mesh.vertices()[k][2], ref);
  }

  // write the elements
  write_simplex(fid, mesh.lines());
  write_simplex(fid, mesh.triangles());
  write_polygons(fid, mesh);
  GmfCloseMesh(fid);
}
}  // namespace meshb

namespace obj {

std::string get_base_dir(const std::string& filepath) {
  if (filepath.find_last_of("/\\") != std::string::npos)
    return filepath.substr(0, filepath.find_last_of("/\\"));
  return "";
}

void read(const std::string& filename, Mesh& mesh) {
  tinyobj::attrib_t attrib;
  std::vector<tinyobj::shape_t> shapes;
  std::vector<tinyobj::material_t> materials;

  std::string base_dir = get_base_dir(filename.c_str());
  if (base_dir.empty()) {
    base_dir = ".";
  }
  base_dir += "/";

  std::string warn;
  std::string err;
  bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err,
                              filename.c_str(), base_dir.c_str(), false);
  if (!warn.empty()) {
    std::cout << "WARN: " << warn << std::endl;
  }
  if (!err.empty()) {
    std::cerr << err << std::endl;
  }
  ASSERT(ret) << "failed to load " << filename.c_str();

  LOG << fmt::format("# of vertices  = {}", attrib.vertices.size() / 3);
  LOG << fmt::format("# of normals   = {}", attrib.normals.size() / 3);
  LOG << fmt::format("# of texcoords = {}", attrib.texcoords.size() / 2);
  LOG << fmt::format("# of materials = {}", materials.size());
  LOG << fmt::format("# of shapes    = {}", shapes.size());

  // read the vertices
  mesh.vertices().set_dim(3);
  mesh.vertices().reserve(attrib.vertices.size() / 3);
  for (size_t v = 0; v < attrib.vertices.size() / 3; v++) {
    mesh.vertices().add(&attrib.vertices[3 * v]);
  }

  // read the texture coordinates
  array2d<coord_t> texcoord(2);
  texcoord.reserve(attrib.texcoords.size() / 2);
  for (size_t t = 0; t < attrib.texcoords.size() / 2; t++) {
    texcoord.add(&attrib.texcoords[2 * t]);
  }

  // read the shapes
  std::vector<index_t> indices;
  for (size_t i = 0; i < shapes.size(); i++) {
    size_t index_offset = 0;

    ASSERT(shapes[i].mesh.num_face_vertices.size() ==
           shapes[i].mesh.material_ids.size());
    ASSERT(shapes[i].mesh.num_face_vertices.size() ==
           shapes[i].mesh.smoothing_group_ids.size());

    // for each face
    for (size_t f = 0; f < shapes[i].mesh.num_face_vertices.size(); f++) {
      size_t fnum = shapes[i].mesh.num_face_vertices[f];
      indices.resize(fnum);
      for (size_t j = 0; j < fnum; j++)
        indices[j] = shapes[i].mesh.indices[index_offset + j].vertex_index;

      if (fnum == 3) {
        mesh.triangles().add(indices.data());
      } else if (fnum == 4) {
        mesh.quads().add(indices.data());
      } else {
        mesh.polygons().add(indices.data(), indices.size());
      }

      index_offset += fnum;
    }
  }  // loop over shape

// initialize the field of texture coordinates (after building the mesh)
#if 0
  mesh.fields().add("texcoord", 2, 1);

  Field& field = mesh.fields().fields()["texcoord"];
  field.triangles().reserve(mesh.triangles().n());
  // field.quads().reserve(mesh.quads().n());

  std::vector<coord_t> data;
  for (size_t i = 0; i < shapes.size(); i++) {
    size_t index_offset = 0;

    // for each face
    for (size_t f = 0; f < shapes[i].mesh.num_face_vertices.size(); f++) {
      size_t fnum = shapes[i].mesh.num_face_vertices[f];
      data.resize(2 * fnum);
      for (size_t j = 0; j < fnum; j++) {
        data[0 * fnum + j] =
            texcoord[shapes[i].mesh.indices[index_offset + j].texcoord_index]
                    [0];
        data[1 * fnum + j] =
            texcoord[shapes[i].mesh.indices[index_offset + j].texcoord_index]
                    [1];
      }

      if (fnum == 3) {
        field.triangles().add(data.data());
      } else if (fnum == 4) {
        // field.quads().add(data.data());
      } else {
        field.polygons().add(data.data(), data.size());
      }

      index_offset += fnum;
    }
  }  // loop over shape

  ASSERT(field.triangles().n() == mesh.triangles().n());
#endif
}

void write(const Mesh& mesh, const std::string& filename) {
  FILE* out = fopen(filename.c_str(), "w");
  ASSERT(out) << "failed to open file: " << filename;

  // header
  fprintf(out, "# obj mesh from terra\n");

  // write vertices
  for (int k = 0; k < mesh.vertices().n(); k++) {
    const double* p = mesh.vertices()[k];
    double z = (mesh.vertices().dim() == 2) ? 0.0 : p[2];
    fprintf(out, "v %.10f %.10f %.10f\n", p[0], p[1], z);
  }

  // write triangles
  for (int k = 0; k < mesh.triangles().n(); k++) {
    fprintf(out, "f ");
    for (int j = 0; j < 3; j++) {
      int idx = mesh.triangles()[k][j];
      fprintf(out, "%d ", idx + 1);
    }
    fprintf(out, "\n");
  }

  // write quads
  for (int k = 0; k < mesh.quads().n(); k++) {
    fprintf(out, "f ");
    for (int j = 0; j < 4; j++) {
      int idx = mesh.quads()[k][j];
      fprintf(out, "%d ", idx + 1);
    }
    fprintf(out, "\n");
  }

  if (mesh.polygons().n() > 0)
    std::cout << "[warning] skipping polygons for now" << std::endl;

  fclose(out);
}
}  // namespace obj

namespace shp {

static const std::map<int, std::string> shpType2Name = {
    {1, "point"}, {3, "arcs"}, {5, "polygon"}, {8, "multipoint"}};

void read(const std::string& filename, Mesh& mesh) {
  SHPHandle handle = SHPOpen(filename.c_str(), "r");

  int n_entities, shape_type;
  double adf_min[4], adf_max[4];
  SHPGetInfo(handle, &n_entities, &shape_type, adf_min, adf_max);
  int shx = (handle->fpSHX) ? *handle->fpSHX : 0;

  // all natural-earth-vector data is represented in WGS84
  // https://www.naturalearthdata.com/features/
  LOG << fmt::format("# entities = {} shape type = {}, SHX = {}", n_entities,
                     shape_type, shx);

  std::map<int, std::array<double, 2>> border;
  int group = 0;
  for (int i = 0; i < n_entities; i++) {
    SHPObject* object = SHPReadObject(handle, i);
    ASSERT(object);
    int n_vertices = object->nVertices;
    int n_parts = object->nParts;
    int type = object->nSHPType;
    // LOG << fmt::format(
    //     "reading object {} with {} parts with {} vertices, type = {}", i,
    //     n_parts, n_vertices, shpType2Name.at(type));

    const double* x = object->padfX;
    const double* y = object->padfY;
    const double* z = object->padfZ;
    // const double* m = object->padfM;

    // add vertices
    double tol = 1e-6;
    int offset = mesh.vertices().n();
    mesh.vertices().reserve(mesh.vertices().n() + n_vertices);
    for (int j = 0; j < n_vertices; j++) {
      double lon = M_PI - M_PI * x[j] / 180;
      double lat = M_PI * y[j] / 180;
      double X = cos(lat) * cos(lon);
      double Y = cos(lat) * sin(lon);
      double Z = sin(lat);
      // const double point[3] = {x[j], y[j], z[j]};
      const double point[3] = {X, Y, Z};
      if (std::fabs(x[j] - 180) < tol || std::fabs(180 - x[j]) < tol)
        border.insert({j + offset, {lon, lat}});
      mesh.vertices().add(point);
    }

    // add parts
    const int* part_start = object->panPartStart;
    for (int k = 0; k < n_parts; k++) {
      group++;
      int j0 = part_start[k];
      int j1 = (k + 1 == n_parts) ? n_vertices : part_start[k + 1];
      int e0 = j0 + offset;
      int first = mesh.lines().n();
      for (int j = j0; j < j1; j++) {
        int e1 = (j + 1 == j1) ? e0 : j + 1 + offset;
        int edge[2] = {e0, e1};
        int n = mesh.lines().n();
        mesh.lines().add(edge);
        mesh.lines().set_group(n, group);
        e0 = e1;
      }

      int last = mesh.lines().n() - 1;
      vec3d p(mesh.vertices()[mesh.lines()(first, 0)]);
      vec3d q(mesh.vertices()[mesh.lines()(last, 1)]);
      double d = length(q - p);
      if (d < 1e-12) {
        mesh.lines()(last, 1) = mesh.lines()(first, 0);
      } else
        LOG << fmt::format("entity {}, part {}, distance = {}", i, k, d);
      // ASSERT(mesh.lines()(last, 1) == mesh.lines()(first, 0));
    }
  }
  LOG << fmt::format("# border vertices = {}", border.size());
  for (auto& p : border) {
    LOG << fmt::format("vtx = {}, lon = {}, lat = {}", p.first, p.second[0],
                       p.second[1]);
  }

  SHPClose(handle);
}
}  // namespace shp

}  // namespace vortex
