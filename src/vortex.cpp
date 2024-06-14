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
#include <fmt/format.h>

#include <argparse/argparse.hpp>
#include <queue>
#include <unordered_set>

#include "graphics.h"
#include "halfedges.h"
#include "io.h"
#include "library.h"
#include "log.h"
#include "mesher.h"
#include "numerics.h"
#include "particles.h"
#include "texture.h"
#include "util.h"
#include "voronoi.h"

#ifndef VORTEX_SOURCE_DIR
#define VORTEX_SOURCE_DIR "./"
#endif

namespace vortex {

namespace {

void run_visualizer(argparse::ArgumentParser &program) {
  std::string filename = program.get<std::string>("input");

  Mesh mesh(3);
  read_mesh(filename, mesh);

  // randomize polygon colors a bit, otherwise neighboring cells
  // will have similar colors and won't visually stand out
  size_t n_colors = 20;
  int n_groups = 1;
  for (int k = 0; k < mesh.polygons().n(); k++)
    n_groups = std::max(n_groups, mesh.polygons().group(k));
  std::vector<int> site2color(n_groups);
  for (size_t k = 0; k < n_groups; k++)
    site2color[k] = int(n_colors * double(rand()) / double(RAND_MAX));
  for (size_t k = 0; k < mesh.polygons().n(); k++) {
    int group = mesh.polygons().group(k); // the group is the site
    mesh.polygons().set_group(k, site2color[group]);
  }

  mesh.fields().set_defaults(mesh);
  Viewer viewer(mesh, 7681);
}

void apply_mask(const std::string &input, double tmin, double tmax,
                Mesh &mesh) {
  TextureOptions tex_opts;
  tex_opts.format = TextureFormat::kGrayscale;
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
    auto *t = mesh.triangles()[k];
    for (int j = 0; j < 3; j++) {
      if (mesh.vertices().group(t[j]) == 2)
        mesh.triangles().set_group(k, 2);
    }
  }
}

void run_mesher(argparse::ArgumentParser &program) {
  TextureOptions tex_opts;
  tex_opts.format = TextureFormat::kGrayscale;
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
  if (!output_metric.empty())
    texture.write(output_metric);

  double hmin = program.get<double>("--hmin");
  double hmax = program.get<double>("--hmax");
  if (hmax < 0)
    hmax = 2.0 * hmin;
  int n_iter = program.get<int>("--num_iter");
  std::string output_filename = program.get<std::string>("--output");
  ASSERT(!output_filename.empty());

  LOG << fmt::format("generating mesh (# iter = {}) with hmin = {}, hmax = {}",
                     n_iter, hmin, hmax);

  MeshingParameters msh_opts;
  msh_opts.max_iter = n_iter;
  msh_opts.h_min = hmin;
  msh_opts.h_max = hmax;
  EarthMesher mesher(texture);
  mesher.generate(msh_opts);

  LOG << "writing " << output_filename;
  Mesh output_mesh(3);
  mesher.extract(output_mesh);

  // assign groups to cells
  apply_mask(filename, tmin, tmax, output_mesh);
  meshb::write(output_mesh, output_filename);
}

void run_extract(argparse::ArgumentParser &program) {
  // read input mesh
  Mesh mesh(3);
  const auto &input = program.get<std::string>("input");
  read_mesh(input, mesh);

  //
  const auto &input_mask = program.get<std::string>("--mask");
  if (input_mask != "groups") {
    // read the image and apply the mask
    double tmin = 10, tmax = 255; // TODO make user inputs
    apply_mask(input_mask, tmin, tmax, mesh);
  }

  // create a halfmesh (there should be no boundary)
  HalfMesh hmesh(mesh);

  // build a list of boundary edges to process
  std::unordered_set<half_t> edges;
  edges.reserve(hmesh.edges().size());
  for (auto &e : hmesh.edges()) {
    int gl = e.get_face().group();
    int gr = e.get_twin().get_face().group();
    if (gl < gr)
      edges.insert(e.index());
  }
  LOG << fmt::format("detected {} edges on boundary", edges.size());
  auto bnd = edges; // make a copy for separating land + water later

  // utility to add an edge with a group
  Mesh water(3), land(3);
  auto add_line = [](Mesh &mesh, int *e, int g) {
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
    if (eit != edges.end())
      edges.erase(eit);

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
        if (g_next != group)
          break;
        current = next;
      };

      // remove the edge and twin from the list
      eit = edges.find(current);
      if (eit != edges.end())
        edges.erase(eit);
      eit = edges.find(hmesh.edges()[current].twin());
      if (eit != edges.end())
        edges.erase(eit);

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
      HalfFace &face = hmesh.faces()[f];
      face.set_group(i);
      // mesh.polygons().set_group(f, i);
      //  mesh.triangles().set_group(f, i);
      visited[f] = true;
      land_faces.insert(f);

      // loop through the neighbors
      e = face.edge();
      for (int j = 0; j < face.n(); j++) {
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
  auto arg_water = program.get<std::string>("--oceans");
  meshb::write(water, arg_water);

  hmesh.activate_faces_by([&land_faces](half_t k) {
    return land_faces.find(k) != land_faces.end();
  });
  hmesh.extract(land);
  auto arg_land = program.get<std::string>("--continents");
  meshb::write(land, arg_land);
}

void run_voronoi(argparse::ArgumentParser &program) {
  auto arg_domain = program.get<std::string>("--domain");
  auto arg_points = program.get<std::string>("--points");
  auto resolution = program.get<double>("--resolution");
  size_t n_points;
  double earth_area = 4 * M_PI * std::pow(6378, 2);
  if (arg_domain == "sphere" && resolution > 0) {
    n_points = (int)(earth_area / std::pow(resolution, 2));
  } else {
    n_points = program.get<int>("--n_points");
  }
  auto n_smooth = program.get<int>("--n_smooth");

  // set up the mesh if using a triangle mesh
  bool use_mesh = true;
  Mesh background_mesh(3);
  if (arg_domain == "sphere" || arg_domain == "square") {
    // nothing to prepare
    use_mesh = false;
  } else if (arg_domain == "icosahedron") {
    SubdividedIcosahedron sphere(program.get<int>("--n_subdiv"));
    sphere.vertices().copy(background_mesh.vertices());
    sphere.triangles().copy(background_mesh.triangles());
  } else {
    read_mesh(arg_domain, background_mesh);
  }
  if (use_mesh)
    LOG << "# triangles in mesh = " << background_mesh.triangles().n();

  auto irand = [](int min, int max) {
    return min + double(rand()) / (double(RAND_MAX) + 1.0) * (max - min);
  };

  // TODO, should this program also accept weights for SDOT?
  int dim = 3;
  Vertices sample(3);
  if (arg_points == "vertices") {
    background_mesh.vertices().copy(sample);
  } else if (arg_points == "random") {
    sample.reserve(n_points);
    if (arg_domain == "sphere") {
      double x[3];
      for (size_t k = 0; k < n_points; k++) {
        coord_t theta = 2.0 * M_PI * irand(0, 1);
        coord_t phi = acos(2.0 * irand(0, 1) - 1.0);
        x[0] = cos(theta) * sin(phi);
        x[1] = sin(theta) * sin(phi);
        x[2] = cos(phi);
        sample.add(x);
      }
    } else if (arg_domain == "square") {
      double x[3] = {0, 0, 0};
      for (size_t k = 0; k < n_points; k++) {
        for (int d = 0; d < 2; d++)
          x[d] = double(rand()) / double(RAND_MAX);
        sample.add(x);
      }
    } else
      sample_surface(background_mesh, sample, n_points);
  } else if (arg_points == "random_oceans") {
    sample.reserve(n_points);
    std::string tex_file =
        std::string(VORTEX_SOURCE_DIR) + "/../data/oceans_2048.png";
    TextureOptions tex_opts;
    tex_opts.format = TextureFormat::kGrayscale;
    Texture texture(tex_file, tex_opts);
    texture.make_binary(10, 10, 255);

    vec3d x, uv;
    while (sample.n() < n_points) {
      coord_t theta = 2.0 * M_PI * irand(0, 1);
      coord_t phi = acos(2.0 * irand(0, 1) - 1.0);
      x[0] = cos(theta) * sin(phi);
      x[1] = sin(theta) * sin(phi);
      x[2] = cos(phi);

      // calculate (u, v) consistent with other algorithms
      sphere_params(x, uv);
      double t;
      texture.sample(uv[0], uv[1], &t);
      if (t < 50)
        continue;

      sample.add(&x[0]);
    }
  } else {
    Mesh tmp(dim);
    read_mesh(arg_points, tmp);
    tmp.vertices().copy(sample);
  }
  n_points = sample.n();
  LOG << fmt::format("initialized {} points", n_points);

  std::vector<index_t> order(n_points);
  auto ordering = program.get<std::string>("--reorder");
  if (ordering == "morton")
    sort_points_on_zcurve(sample[0], n_points, dim, order);
  else {
    for (size_t k = 0; k < n_points; k++)
      order[k] = k;
  }

  Vertices points(dim);
  points.reserve(n_points);
  coord_t x[dim];
  for (size_t i = 0; i < n_points; i++) {
    for (int d = 0; d < dim; d++)
      x[d] = sample[order[i]][d];
    points.add(x);
  }

  VoronoiDiagram voronoi(dim, points[0], n_points);
  VoronoiDiagramOptions options;
  options.n_neighbors = program.get<int>("--n_neighbors");
  options.allow_reattempt = false;
  options.parallel = true;

  auto quiet = program.get<bool>("--quiet");
  auto verbose = program.get<bool>("--verbose");
  auto save = program.present<std::string>("--output");
  auto on_sphere = program.get<bool>("--on_sphere");
  auto calculate_voronoi_diagram = [&voronoi, &options, &points, n_smooth, save,
                                    quiet, verbose, on_sphere](auto &domain) {
    int n_iter = n_smooth;
    for (int iter = 1; iter <= n_iter; ++iter) {
      options.store_mesh = (iter == n_iter) && save;
      options.verbose = (verbose || iter == 1 || iter == n_iter) && !quiet;
      voronoi.vertices().clear();
      voronoi.vertices().set_dim(3);
      voronoi.polygons().clear();
      voronoi.triangles().clear();
      voronoi.compute(domain, options);

      // move each site to the centroid of the corresponding cell
      voronoi.smooth(points, on_sphere);
      auto props = voronoi.analyze();
      if (!quiet)
        LOG << fmt::format("iter = {}, area = {}", iter, props.area);
    }
  };

  // calculate!
  if (arg_domain == "sphere") {
    SphereDomain domain;
    if (n_points <= 10000)
      domain.set_initialization_fraction(0.7);
    calculate_voronoi_diagram(domain);
  } else if (arg_domain == "square") {
    SquareDomain domain;
    calculate_voronoi_diagram(domain);
  } else {
    const auto *p = background_mesh.vertices()[0];
    size_t np = background_mesh.vertices().n();
    const auto *t = background_mesh.triangles()[0];
    size_t nt = background_mesh.triangles().n();
    TriangulationDomain domain(p, np, t, nt);
    calculate_voronoi_diagram(domain);
  }
  if (save)
    voronoi.merge();

  if (program.present<std::string>("--output")) {
    LOG << fmt::format("writing {} polygons", voronoi.polygons().n());
    auto arg_output = program.get<std::string>("--output");
    if (voronoi.polygons().n() > 0)
      meshb::write(voronoi, arg_output);
  }
  if (program.present<std::string>("--output_points")) {
    Mesh tmp(dim);
    points.copy(tmp.vertices());
    meshb::write(tmp, program.get<std::string>("--output_points"));
  }
}

void run_merge(argparse::ArgumentParser &program) {
  auto arg_input = program.get<std::string>("input");
  auto arg_output = program.get<std::string>("--output");
  auto combine = program.get<bool>("--combine");

  double tol = 1e-10;
  Mesh input_mesh(3);
  read_mesh(arg_input, input_mesh);
  input_mesh.merge(tol);

  if (!combine) {
    meshb::write(input_mesh, arg_output);
    return;
  }

  Mesh output_mesh(input_mesh.vertices().dim());
  input_mesh.vertices().copy(output_mesh.vertices());
  input_mesh.separate_polygons_into_connected_components(
      output_mesh.polygons());
  meshb::write(output_mesh, arg_output);
}

void run_simulation(argparse::ArgumentParser &program) {
  const int dim = 3;
  size_t n_points = program.get<int>("--n_particles");
  auto corners = program.get<std::vector<double>>("--corners");

  Vertices sample(3);
  auto arg_points = program.get<std::string>("--particles");
  auto arg_domain = program.get<std::string>("--domain");
  auto irand = [](int min, int max) {
    return min + double(rand()) / (double(RAND_MAX) + 1.0) * (max - min);
  };
  if (arg_points == "random") {
    sample.reserve(n_points);
    if (arg_domain == "sphere") {
      double x[3];
      for (size_t k = 0; k < n_points; k++) {
        coord_t theta = 2.0 * M_PI * irand(0, 1);
        coord_t phi = acos(2.0 * irand(0, 1) - 1.0);
        x[0] = cos(theta) * sin(phi);
        x[1] = sin(theta) * sin(phi);
        x[2] = cos(phi);
        sample.add(x);
      }
    } else if (arg_domain == "rectangle") {
      SquareDomain domain({corners[0], corners[1], 0},
                          {corners[2], corners[3], 0});
      for (size_t k = 0; k < n_points; k++) {
        auto x = domain.random_point();
        // for (int d = 0; d < 2; d++) x[d] = double(rand()) / double(RAND_MAX);
        sample.add(&x[0]);
      }
    } else
      NOT_IMPLEMENTED;
  } else if (arg_points == "random_oceans") {
    sample.reserve(n_points);
    std::string tex_file =
        std::string(VORTEX_SOURCE_DIR) + "/../data/oceans_2048.png";
    TextureOptions tex_opts;
    tex_opts.format = TextureFormat::kGrayscale;
    Texture texture(tex_file, tex_opts);
    texture.make_binary(10, 10, 255);

    vec3d x, uv;
    while (sample.n() < n_points) {
      coord_t theta = 2.0 * M_PI * irand(0, 1);
      coord_t phi = acos(2.0 * irand(0, 1) - 1.0);
      x[0] = cos(theta) * sin(phi);
      x[1] = sin(theta) * sin(phi);
      x[2] = cos(phi);

      // calculate (u, v) consistent with other algorithms
      sphere_params(x, uv);
      double t;
      texture.sample(uv[0], uv[1], &t);
      if (t < 50)
        continue;

      sample.add(&x[0]);
    }
  } else {
    Mesh tmp(dim);
    read_mesh(arg_points, tmp);
    tmp.vertices().copy(sample);
  }
  n_points = sample.n();
  LOG << fmt::format("initialized {} points", n_points);

  // construct a better ordering of the points
  std::vector<index_t> order(n_points);
  sort_points_on_zcurve(sample[0], n_points, dim, order);
  Vertices vertices(dim);
  vertices.reserve(n_points);
  coord_t x[dim];
  for (size_t i = 0; i < n_points; i++) {
    for (int d = 0; d < dim; d++)
      x[d] = sample[order[i]][d];
    vertices.add(x);
  }

  auto run_case = [&](auto &domain, const auto &velocity, const auto &density,
                      const auto &force) {
    using Domain_t = typename std::remove_reference<decltype(domain)>::type;

    // smooth the initial point distribution with Lloyd relaxation
    VoronoiDiagram smoother(dim, vertices[0], n_points);
    VoronoiDiagramOptions options;
    options.n_neighbors = 100;
    options.allow_reattempt = false;
    options.parallel = true;
    options.store_facet_data = true;
    int n_iter = program.get<int>("--n_smooth");

    for (int iter = 1; iter <= n_iter; ++iter) {
      options.store_mesh = false;
      options.verbose = false;
      smoother.compute(domain, options); // calculate voronoi diagram
      smoother.smooth(vertices, false);  // move sites to centroids
      for (size_t k = 0; k < vertices.n(); k++)
        project_point<Domain_t>(vertices[k]);
    }
    LOG << "done smoothing points";

    // set up the fluid simulator
    SpringParticles<Domain_t> solver(domain, n_points, vertices[0],
                                     vertices.dim());

    // assign initial velocities
    solver.particles().set_velocity(velocity);
    solver.particles().set_density(density);

    // set up fluid and solver properties
    FluidProperties props;
    SimulationOptions solver_opts;
    solver.initialize(domain, solver_opts);
    LOG << "initialized simulation";

    int nt = program.get<int>("--total_time_steps");
    double hn = solver.voronoi().max_radius();
    solver_opts.epsilon = program.get<double>("--epsilon_scale") * hn;
    solver_opts.time_step = program.get<double>("--time_step_scale") *
                            std::pow(solver_opts.epsilon, 2);
    LOG << fmt::format("hn = {:1.3e}, eps = {:1.3e}, dt = {:1.3e}", hn,
                       solver_opts.epsilon, solver_opts.time_step);
    double force_eps = program.get<double>("--force_epsilon");
    if (force_eps > 0)
      solver_opts.epsilon = force_eps;
    double force_dt = program.get<double>("--force_time_step");
    if (force_dt > 0)
      solver_opts.time_step = force_dt;
    solver_opts.verbose = false;
    solver_opts.backtrack = false;
    solver_opts.reflection_boundary_condition =
        program.get<bool>("--boundary_reflection");
    ASSERT(!solver_opts.reflection_boundary_condition);
    solver_opts.advect_from_centroid = !program.get<bool>("--advect_from_site");
    const auto directory = program.get<std::string>("--output_directory");
    int s = 0;
    Timer timer;
    timer.start();
    for (int t = 0; t < nt; t++) {
      if (t % program.get<int>("--save_every") == 0)
        solver.particles().save(
            fmt::format("{}/particles{}.vtk", directory, s++));
      solver.step(force, props, solver_opts);
      solver_opts.iteration++;
    }
    timer.stop();
    solver.print_footer();
    LOG << fmt::format("done! total time = {} seconds.", timer.seconds());
  };

  const auto density_ratio = program.get<double>("--density_ratio");
  if (arg_domain == "rectangle") {
    typedef SquareDomain Domain_t;
    auto velocity = [](const double *x) -> vec3 { return {0, 0, 0}; };
    auto force = [](const Particle &p) -> vec3d {
      vec3d f;
      f[1] -= 9.81 * p.mass;
      return f;
    };
    auto density = [&](const double *x) -> double {
      double f = -0.2 * cos(M_PI * x[0]);
      return x[1] > f ? density_ratio : 1;
    };
    Domain_t domain({corners[0], corners[1], 0}, {corners[2], corners[3], 0});
    run_case(domain, velocity, density, force);
  } else if (arg_domain == "sphere") {
    typedef SphereDomain Domain_t;
    vec3d omega{0., 0., 0.};
    omega[1] = program.get<double>("--omega");
    auto velocity = [omega](const double *x) -> vec3d {
      vec3d p(x);
      vec3d v = cross(omega, p);
      return v;
    };
    auto force = [omega](const Particle &particle) -> vec3d {
      vec3d f;
      double m = particle.mass;
      const auto &v = particle.velocity;
      const auto &p = particle.position;
      f = -2 * m * omega[1] * p[1] * cross(p, v);
      return f;
    };
    auto density = [](const double *x) -> double {
      return std::fabs(x[1]) > 0.25 ? 1 : 10;
    };
    Domain_t domain;
    run_case(domain, velocity, density, force);
  }
}

} // namespace

} // namespace vortex

int main(int argc, char **argv) {
  argparse::ArgumentParser program("vortex", "1.0");

  argparse::ArgumentParser cmd_viz("viz");
  cmd_viz.add_description("visualize a mesh");
  cmd_viz.add_argument("input").help("input file");
  program.add_subparser(cmd_viz);

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
      .help("how to separate continents and oceans: can be an image (e.g. "
            "oceans_2048.png) or 'groups' to use existing grouping in input")
      .default_value("groups");
  cmd_extract.add_argument("--oceans")
      .help("output mesh file for oceans ([prefer] .meshb or .obj)")
      .default_value("oceans.meshb");
  cmd_extract.add_argument("--continents")
      .help("output mesh file for continents ([prefer] .meshb or .obj)")
      .default_value("continents.meshb");
  program.add_subparser(cmd_extract);

  argparse::ArgumentParser cmd_voronoi("voronoi");
  cmd_voronoi.add_description("calculate Voronoi diagram on surface");
  cmd_voronoi.add_argument("--domain")
      .help("input surface: (.obj or .meshb), sphere or icosahedron");
  cmd_voronoi.add_argument("--n_subdiv")
      .help("number of subdivisions of the sphere mesh (only applicable to "
            "icosahedron domain)")
      .default_value(4)
      .scan<'i', int>();
  cmd_voronoi.add_argument("--points")
      .help("sampling technique to use (random, vertices, or specify mesh file "
            "to use "
            "vertices)")
      .default_value("random");
  cmd_voronoi.add_argument("--n_points")
      .help("# Voronoi sites, only applicable for random --points option")
      .default_value(10000)
      .scan<'i', int>();
  cmd_voronoi.add_argument("--resolution")
      .help("an approximate size of the cell in kilometers, only applicable "
            "for sphere --domain option")
      .default_value(0.0)
      .scan<'g', double>();
  cmd_voronoi.add_argument("--n_smooth")
      .help("# iterations of Lloyd relaxation")
      .default_value(1)
      .scan<'i', int>();
  cmd_voronoi.add_argument("--output")
      .help("output mesh file ([prefer] .meshb or .obj)");
  cmd_voronoi.add_argument("--output_points")
      .help("output points filename ([prefer] .meshb or .obj)");
  cmd_voronoi.add_argument("--n_neighbors")
      .help("number of nearest neighbors to use when calculating the Voronoi "
            "diagram")
      .default_value(50)
      .scan<'i', int>();
  cmd_voronoi.add_argument("--verbose")
      .help("print timing information at each smoothing iteration")
      .flag();
  cmd_voronoi.add_argument("--quiet").flag();
  cmd_voronoi.add_argument("--on_sphere")
      .help("project sites to the sphere during smoothing")
      .flag();
  cmd_voronoi.add_argument("--reorder")
      .help("ordering method (morton, none)")
      .default_value("morton");
  program.add_subparser(cmd_voronoi);

  argparse::ArgumentParser cmd_merge("merge");
  cmd_merge.add_description(
      "merge nearby vertices in the mesh, renumbering mesh elements");
  cmd_merge.add_argument("input").help("path to input mesh file (.meshb)");
  cmd_merge.add_argument("--combine")
      .help("combine polygons with the same group and separate them into "
            "connected components (useful for Voronoi diagrams restricted to a "
            "triangulation)")
      .flag();
  cmd_merge.add_argument("--output").help("path to output mesh file (.meshb)");
  program.add_subparser(cmd_merge);

  argparse::ArgumentParser cmd_sim("simulate");
  cmd_sim.add_description("simulate fluid flow");
  cmd_sim.add_argument("--domain").help("simulation domain: rectangle, sphere");
  cmd_sim.add_argument("--particles")
      .help("sampling technique to use for initial particle positions (random, "
            "random_oceans)")
      .default_value("random");
  cmd_sim.add_argument("--n_particles")
      .help("# particles to use in the simulation")
      .default_value(10000)
      .scan<'i', int>();
  cmd_sim.add_argument("--resolution")
      .help("an approximate size of the cell in kilometers, only applicable "
            "for sphere --domain option")
      .default_value(100.0)
      .scan<'g', double>();
  cmd_sim.add_argument("--omega")
      .help("y-component of rotational velocity")
      .default_value(0)
      .scan<'g', double>();
  cmd_sim.add_argument("--boundary_reflection")
      .help("option to add reflection boundary conditions (only for the "
            "rectangular domain)")
      .flag();
  cmd_sim.add_argument("--density_ratio")
      .default_value(10.0)
      .scan<'g', double>();
  cmd_sim.add_argument("--n_smooth")
      .help("# iterations of Lloyd relaxation for initial points")
      .default_value(100)
      .scan<'i', int>();
  cmd_sim.add_argument("--epsilon_scale")
      .help(
          "scaling factor used to compute the inverse spring coefficient eps = "
          "epsilon_scale * hn")
      .default_value(10.0)
      .scan<'g', double>();
  cmd_sim.add_argument("--time_step_scale")
      .help(
          "scaling factor used to compute the time step: dt = time_step_scale "
          "* eps * eps")
      .default_value(0.15)
      .scan<'g', double>();
  cmd_sim.add_argument("--output_directory")
      .help("output location for .vtk files")
      .default_value("");
  cmd_sim.add_argument("--save_every")
      .help("# time steps after which the solution is saved")
      .default_value(50)
      .scan<'i', int>();
  cmd_sim.add_argument("--total_time_steps")
      .help("total # time steps in simulation")
      .default_value(100)
      .scan<'i', int>();
  cmd_sim.add_argument("--advect_from_site").flag();
  cmd_sim.add_argument("--force_time_step")
      .default_value(-1.0)
      .scan<'g', double>();
  cmd_sim.add_argument("--force_epsilon")
      .default_value(-1.0)
      .scan<'g', double>();
  cmd_sim.add_argument("--corners")
      .help("xmin ymin xmax ymax")
      .nargs(4)
      .default_value(std::vector<double>{0.0, 0.0, 1.0, 1.0})
      .scan<'g', double>();

  program.add_subparser(cmd_sim);

  try {
    program.parse_args(argc, argv);
  } catch (const std::runtime_error &err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    std::exit(1);
  }

  if (program.is_subcommand_used("viz")) {
    vortex::run_visualizer(program.at<argparse::ArgumentParser>("viz"));
  } else if (program.is_subcommand_used("mesh")) {
    vortex::run_mesher(program.at<argparse::ArgumentParser>("mesh"));
  } else if (program.is_subcommand_used("extract")) {
    vortex::run_extract(program.at<argparse::ArgumentParser>("extract"));
  } else if (program.is_subcommand_used("voronoi")) {
    vortex::run_voronoi(program.at<argparse::ArgumentParser>("voronoi"));
  } else if (program.is_subcommand_used("merge")) {
    vortex::run_merge(program.at<argparse::ArgumentParser>("merge"));
  } else if (program.is_subcommand_used("simulate")) {
    vortex::run_simulation(program.at<argparse::ArgumentParser>("simulate"));
  } else {
    std::cout << program.help().str();
  }

  return 0;
}
