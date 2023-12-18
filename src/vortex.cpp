
#include <fmt/format.h>

#include <argparse/argparse.hpp>
#include <queue>
#include <unordered_set>

#include "halfedges.h"
#include "io.h"
#include "log.h"
#include "mesher.h"
#include "numerics.h"
#include "texture.h"

namespace vortex {

namespace {

void apply_mask(const std::string& input, double tmin, double tmax,
                Mesh& mesh) {
  TextureOptions tex_opts{.format = TextureFormat::kGrayscale};
  Texture mask(input, tex_opts);
  mask.make_binary(tmin, tmin, tmax);

  // tag vertices and triangles
  for (size_t k = 0; k < mesh.vertices().n(); k++) {
    vec3d p(mesh.vertices()[k]);
    vec3d u;
    sphere_params(p, u);
    double t;
    mask.sample(u[0], u[1], &t);
    if (t <= tmin)
      mesh.vertices().set_group(k, 1);
    else
      mesh.vertices().set_group(k, 2);
  }

  for (size_t k = 0; k < mesh.triangles().n(); k++) {
    mesh.triangles().set_group(k, 1);
    auto* t = mesh.triangles()[k];
    for (int j = 0; j < 3; j++) {
      if (mesh.vertices().group(t[j]) == 2) mesh.triangles().set_group(k, 2);
    }
  }
}

void run_mesher(argparse::ArgumentParser& program) {
  TextureOptions tex_opts{.format = TextureFormat::kGrayscale};
  std::string filename = program.get<std::string>("metric");
  bool no_smooth = program.get<bool>("--no_smooth");
  Texture texture(filename, tex_opts);
  double tmin = 10, tmax = 255;
  if (!no_smooth) {
    LOG << "smoothing image ...";
    texture.make_binary(tmin, tmin, tmax);
    texture.smooth(10);
  }
  texture.make_periodic();
  std::string output_metric = program.get<std::string>("--output_metric");
  if (!output_metric.empty()) texture.write(output_metric);

  double hmin = program.get<double>("--hmin");
  double hmax = program.get<double>("--hmax");
  if (hmax < 0) hmax = 2.0 * hmin;
  int n_iter = program.get<int>("--num_iter");
  std::string output_filename = program.get<std::string>("--output");
  ASSERT(!output_filename.empty());

  LOG << fmt::format("generating mesh (# iter = {}) with hmin = {}, hmax = {}",
                     n_iter, hmin, hmax);

  MeshingParameters msh_opts{.max_iter = n_iter, .h_min = hmin, .h_max = hmax};
  EarthMesher mesher(texture);
  mesher.generate(msh_opts);

  LOG << "writing " << output_filename;
  Mesh output_mesh(3);
  mesher.extract(output_mesh);

  // assign groups to cells
  apply_mask(filename, tmin, tmax, output_mesh);
  meshb::write(output_mesh, output_filename);
}

void run_extract(argparse::ArgumentParser& program) {
  // read input mesh
  Mesh mesh(3);
  const auto& input = program.get<std::string>("input");
  read_mesh(input, mesh);

  //
  const auto& input_mask = program.get<std::string>("--mask");
  if (input_mask != "groups") {
    // read the image and apply the mask
    double tmin = 10, tmax = 255;  // TODO make user inputs
    apply_mask(input_mask, tmin, tmax, mesh);
  }

  // create a halfmesh (there should be no boundary)
  HalfMesh hmesh(mesh);

  // build a list of boundary edges to process
  std::unordered_set<half_t> edges;
  edges.reserve(hmesh.edges().size());
  for (auto& e : hmesh.edges()) {
    int gl = e.get_face().group();
    int gr = e.get_twin().get_face().group();
    if (gl < gr) edges.insert(e.index());
  }
  LOG << fmt::format("detected {} edges on boundary", edges.size());
  auto bnd = edges;  // make a copy for separating land + water later

  // utility to add an edge with a group
  Mesh water(3), land(3);
  auto add_line = [](Mesh& mesh, int* e, int g) {
    size_t id = mesh.lines().n();
    mesh.lines().add(e);
    mesh.lines().set_group(id, g);
  };

  // form loops
  int n_group = 0;
  std::vector<half_t> loops;
  while (!edges.empty()) {
    n_group++;
    // get the next available edge to start a loop
    auto eit = edges.begin();
    half_t first = *eit;
    loops.push_back(first);

    // remove the edge + twin from the set of edges to process
    edges.erase(eit);
    eit = edges.find(hmesh.edges()[first].twin());
    if (eit != edges.end()) edges.erase(eit);

    // continue to next edge until we return to the first edge
    half_t current = first;
    int group = hmesh.edges()[current].get_face().group();
    do {
      // add this edge to the lines topology
      int p = hmesh.edges()[current].node();
      int q = hmesh.edges()[current].get_twin().node();
      int line[2] = {p, q};
      add_line(mesh, line, n_group);
      add_line(water, line, n_group);
      add_line(land, line, n_group);

      // iterate through the edges-onering of q
      current = hmesh.edges()[current].twin();
      while (true) {
        half_t next = hmesh.edges()[current].get_twin().next();
        int g_next = hmesh.edges()[next].get_face().group();
        if (g_next != group) break;
        current = next;
      };

      // remove the edge and twin from the list
      eit = edges.find(current);
      if (eit != edges.end()) edges.erase(eit);
      eit = edges.find(hmesh.edges()[current].twin());
      if (eit != edges.end()) edges.erase(eit);

    } while (current != first);
  }
  LOG << fmt::format("formed {} groups", n_group);

  std::vector<bool> visited(hmesh.faces().size(), false);
  std::unordered_set<half_t> land_faces;
  land_faces.reserve(hmesh.faces().size());
  for (size_t i = 0; i < loops.size(); i++) {
    visited.resize(hmesh.faces().size(), false);
    half_t e = loops[i];
    std::queue<half_t> queue;
    queue.push(hmesh.edges()[e].face());
    while (!queue.empty()) {
      half_t f = queue.front();
      queue.pop();
      HalfFace& face = hmesh.faces()[f];
      face.set_group(i);
      mesh.triangles().set_group(f, i);
      visited[f] = true;
      land_faces.insert(f);

      // loop through the neighbors
      e = face.edge();
      for (int j = 0; j < 3; j++) {
        // check if this is an edge of the boundary
        half_t t = hmesh.edges()[e].twin();
        if (bnd.find(e) != bnd.end() || bnd.find(t) != bnd.end()) {
          e = hmesh.edges()[e].next();
          continue;
        }
        // check if the opposite face has been visited
        half_t fj = hmesh.edges()[e].get_twin().face();
        if (!visited[fj]) {
          visited[fj] = true;
          queue.push(fj);
        }
        e = hmesh.edges()[e].next();
      }
    }
  }

  hmesh.activate_faces_by([&land_faces](half_t k) {
    return land_faces.find(k) == land_faces.end();
  });
  hmesh.extract(water);
  meshb::write(water, "water.meshb");

  hmesh.activate_faces_by([&land_faces](half_t k) {
    return land_faces.find(k) != land_faces.end();
  });
  hmesh.extract(land);
  meshb::write(land, "land.meshb");
}
}  // namespace
}  // namespace vortex

int main(int argc, char** argv) {
  argparse::ArgumentParser program("vortex", "1.0");

  argparse::ArgumentParser cmd_mesh("mesh");
  cmd_mesh.add_description("generate mesh on the unit sphere");
  cmd_mesh.add_argument("metric").help("input metric image (.png, .jpg)");
  cmd_mesh.add_argument("-o", "--output").help("output mesh file");
  cmd_mesh.add_argument("--output_metric")
      .help("output metric file")
      .default_value("");
  cmd_mesh.add_argument("--no_smooth")
      .help("disable smoothing of input metric")
      .default_value(false)
      .implicit_value(true);
  cmd_mesh.add_argument("--hmin")
      .help("minimum mesh size")
      .default_value(0.005)
      .scan<'g', double>();
  cmd_mesh.add_argument("--hmax")
      .help("maximum mesh size")
      .default_value(-1.0)
      .scan<'g', double>();
  cmd_mesh.add_argument("--num_iter")
      .help("# adaptation iterations")
      .default_value(5)
      .scan<'d', int>();
  program.add_subparser(cmd_mesh);

  argparse::ArgumentParser cmd_extract("extract");
  cmd_extract.add_description("extract ocean and continents meshes");
  cmd_extract.add_argument("input").help("input mesh file (.obj, .meshb)");
  cmd_extract.add_argument("--mask")
      .help(
          "how to separate continents and oceans: can be an image (e.g. "
          "oceans_2048.png) or 'groups' to use existing grouping in input")
      .default_value("groups");
  cmd_extract.add_argument("--oceans")
      .help("output mesh file for oceans ([prefer] .meshb or .obj)")
      .default_value("oceans.meshb");
  cmd_extract.add_argument("--continents")
      .help("output mesh file for continents ([prefer] .meshb or .obj)")
      .default_value("continents.meshb");
  program.add_subparser(cmd_extract);

  try {
    program.parse_args(argc, argv);
  } catch (const std::runtime_error& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    std::exit(1);
  }

  if (program.is_subcommand_used("mesh")) {
    vortex::run_mesher(program.at<argparse::ArgumentParser>("mesh"));
  } else if (program.is_subcommand_used("extract")) {
    vortex::run_extract(program.at<argparse::ArgumentParser>("extract"));
  } else {
    std::cout << program.help().str();
  }

  return 0;
}
