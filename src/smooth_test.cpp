/*
Author: Col McDermott
Scrpit: Test file for voronoi diagram energy caclulation and optimized site
smoothing
*/

#include <Predicates_psm.h>
#include <limits.h>

#include <cmath>
#include <nlopt.hpp>

#include "graphics.h"
#include "io.h"
#include "library.h"
#include "log.h"
#include "math/vec.hpp"
#include "quadrature.h"
#include "quadrature.hpp"
#include "tester.h"
#include "util.h"
#include "voronoi.h"

using namespace vortex;

double energy(VoronoiDiagram &voronoi, int t_type, int order, bool grad,
              int dim) {
  double energy = 0.0;
  Topology<Polygon> const &polygons = voronoi.polygons();

  // iterate through each voronoi cell
  for (int i = 0; i < polygons.n(); i++) {
    coord_t cell_energy = 0.0;

    // get site
    auto site_idx = polygons.group(i);
    vec3d site(voronoi.vertices()[site_idx]);

    // triangulate cell
    const auto *p1 = voronoi.vertices()[polygons(i, 0)];
    for (auto j = 1; j < polygons.length(i) - 1; j++) {
      const auto *p2 = voronoi.vertices()[polygons(i, j)];
      const auto *p3 = voronoi.vertices()[polygons(i, j + 1)];
      auto f = [&site](const vec3d &x) -> coord_t {
        return pow(length(x - site), 2);
      };
      if (t_type == 0) {  // clean this up later
        TriangleQuadrature<Triangle> quad(order);
        cell_energy += quad.integrate(f, p1, p2, p3);
      } else {
        TriangleQuadrature<SphericalTriangle> quad(order);
        cell_energy += quad.integrate(f, p1, p2, p3);
      }
    }
    energy += cell_energy;
  }

  if (grad) {
    double norm = 0;
    VoronoiCellProperties props;
    for (int i = 0; i < polygons.n(); i++) {
      props = voronoi.properties()[i];
      vec3 site(voronoi.vertices()[props.site][0],
                voronoi.vertices()[props.site][1],
                voronoi.vertices()[props.site][2]);
      vec3 grads(0., 0., 0.);
      grads = (2 * ((props.volume * site) - props.moment));
      norm += length(grads);
    }

    LOG << fmt::format("Grad. Norm = {:1.10e}", norm);
  }

  return energy;
}

// helper method for energy_objective
void new_sites(Vertices &vertices, const double *prev_sites, int dim,
               unsigned n_sites) {
  for (auto i = 0; i < n_sites; i++) {
    for (int j = 0; j < dim; j++) {
      vertices[i][j] = prev_sites[(dim * i) + j];
    }
  }
}

// Data for site optimization
template <typename D>
struct energy_data {
  double first_e;
  double prev_e;
  VoronoiDiagram &vor;
  D &domain;
  VoronoiDiagramOptions &options;
  Vertices &verts;
  std::vector<double> &norms;
  double norm;
  int dim;
  int it{0};
};

// Objective function to minimize
double energy_objective(unsigned n, const double *sites, double *grad,
                        void *ed) {
  // get data
  energy_data<SphereDomain> &data_e =
      *static_cast<energy_data<SphereDomain> *>(ed);

  // increase iteration
  data_e.it++;

  // reset sites and recompute diagram
  data_e.vor.vertices().clear();
  data_e.vor.polygons().clear();
  data_e.vor.triangles().clear();
  new_sites(data_e.verts, sites, data_e.dim, n / data_e.dim);
  data_e.vor.compute(data_e.domain, data_e.options);

  // calculate energy and record current guess
  double curr_e = energy(data_e.vor, 0, 4, false, data_e.dim);
  if ((data_e.prev_e > curr_e) &&
      (((data_e.prev_e - curr_e) / data_e.prev_e) * 100) > 0.006) {
    LOG << fmt::format(
        "IT. {} -> Current Energy = {:1.10e}: Energy Decreased by {:2.3f}%",
        data_e.it, curr_e, ((data_e.prev_e - curr_e) / data_e.prev_e) * 100);
  }
  if ((curr_e > data_e.prev_e) &&
      (((curr_e - data_e.prev_e) / curr_e) * 100) > 0.006) {
    LOG << fmt::format(
        "IT. {} -> Current Energy = {:1.10e}: Energy Increased by {:2.3f}%",
        data_e.it, curr_e, ((curr_e - data_e.prev_e) / curr_e) * 100);
  }
  if (((data_e.prev_e > curr_e) &&
       (((data_e.prev_e - curr_e) / data_e.prev_e) * 100) < 0.006) ||
      (((curr_e > data_e.prev_e) &&
        (((curr_e - data_e.prev_e) / curr_e) * 100) < 0.006))) {
    LOG << fmt::format("IT. {} -> Current Energy = {:1.10e}", data_e.it,
                       curr_e);
  }

  // calculate gradient values
  VoronoiCellProperties props;
  if (grad) {
    data_e.norm = 0;
    for (int i = 0; i < (n / data_e.dim); i++) {
      props = data_e.vor.properties()[i];
      vec3 site(data_e.verts[props.site][0], data_e.verts[props.site][1],
                data_e.verts[props.site][2]);
      vec3 grads(0., 0., 0.);
      grads = (2 * ((props.volume * site) - props.moment));
      for (int j = 0; j < data_e.dim; j++) {
        grad[(3 * i) + j] = grads[j];
      }
      data_e.norm += pow(length(grads), 2);
    }
    data_e.norm = sqrt(data_e.norm);
    data_e.norms.push_back(data_e.norm);
    LOG << fmt::format("Grad. Norm = {:1.10e}",
                       data_e.norms[data_e.norms.size() - 1]);
  }

  data_e.prev_e = curr_e;
  return curr_e;
}

UT_TEST_SUITE(smooth_test_suite)

UT_TEST_CASE(optimized_points_square) {
  // set up optimization
  SquareDomain domain;
  VoronoiDiagramOptions options;
  options.n_neighbors = 75;
  options.parallel = true;
  options.verbose = false;
  options.store_mesh = true;
  size_t n_sites = 1e3;
  int dim = 3;
  std::vector<coord_t> sites(n_sites * dim, 0.0);
  for (size_t k = 0; k < n_sites; k++) {
    sites[k * dim + 0] = double(rand()) / double(RAND_MAX);
    sites[k * dim + 1] = double(rand()) / double(RAND_MAX);
    sites[k * dim + 2] = 0.0;
  }
  std::vector<index_t> order(n_sites);
  sort_points_on_zcurve(sites.data(), n_sites, dim, order);
  Vertices vertices(dim);
  vertices.reserve(n_sites);
  coord_t x[dim];
  for (size_t i = 0; i < n_sites; i++) {
    for (int d = 0; d < dim; d++) x[d] = sites[dim * order[i] + d];
    vertices.add(x);
  }
  VoronoiDiagram voronoi(dim, vertices[0], n_sites);
  voronoi.compute(domain, options);

  double first_e = energy(voronoi, 0, 4, false, dim);

  std::vector<double> norms(n_sites, 0.);
  // fill optimization data
  energy_data<SquareDomain> data_e = {
      first_e, first_e, voronoi, domain, options, vertices, norms, 0, dim, 0};
  int n = sites.size();
  nlopt::opt opt_sites(nlopt::LD_LBFGS, n);
  opt_sites.set_min_objective(&energy_objective, static_cast<void *>(&data_e));

  // optimizaition params
  opt_sites.set_xtol_rel(1e-15);  // adjust
  opt_sites.set_ftol_rel(1e-6);

  // set the lower and upper bounds
  std::vector<double> lower_bound(n, -3.);
  std::vector<double> upper_bound(n, 3.);
  opt_sites.set_lower_bounds(lower_bound);
  opt_sites.set_upper_bounds(upper_bound);

  LOG << fmt::format("Initial Energy = {:1.10e}", first_e);
  double e_opt;
  std::vector<double> X = vertices.data();
  nlopt::result result = opt_sites.optimize(X, e_opt);
  LOG << fmt::format(
      "CVT Diagram Energy = {:1.10e} -> {} Iterations to Converge", e_opt,
      data_e.it);
  LOG << fmt::format(
      "Percentage Decrease Relative to Initial Energy = {:2.3f}%",
      ((first_e - e_opt) / first_e) * 100);

  ASSERT(result == nlopt::SUCCESS);
}
UT_TEST_CASE_END(optimized_points_square)

/*
UT_TEST_CASE(optimized_points_sphere) {
  // creating a voronoi diagram on a sphere domain (copied from
  // voronoi_test.cpp)
  auto irand = [](int min, int max) {
    return min + double(rand()) / (double(RAND_MAX) + 1.0) * (max - min);
  };
  static const int dim = 3;
  size_t n_sites = 1e3;
  std::vector<coord_t> sites(n_sites * dim, 0.0);
  for (size_t k = 0; k < n_sites; k++) {
    coord_t theta = 2.0 * M_PI * irand(0, 1);
    coord_t phi = acos(2.0 * irand(0, 1) - 1.0);
    sites[k * dim + 0] = cos(theta) * sin(phi);
    sites[k * dim + 1] = sin(theta) * sin(phi);
    sites[k * dim + 2] = cos(phi);
  }
  std::vector<index_t> order(n_sites);
  sort_points_on_zcurve(sites.data(), n_sites, dim, order);
  Vertices vertices(dim);
  vertices.reserve(n_sites);
  coord_t x[dim];
  for (size_t i = 0; i < n_sites; i++) {
    for (int d = 0; d < dim; d++) x[d] = sites[dim * order[i] + d];
    vertices.add(x);
  }
  SphereDomain domain;
  VoronoiDiagram voronoi(dim, vertices[0], n_sites);
  VoronoiDiagramOptions options;
  options.n_neighbors = 75;
  options.parallel = true;
  options.verbose = false;
  options.store_mesh = true;
  voronoi.compute(domain, options);
  double first_e = energy(voronoi, 1, 4, false, dim);
  LOG << fmt::format("Initial Energy = {:1.10e}", first_e);

  std::vector<double> norms(n_sites, 0.);
  // fill optimization data
  energy_data<SphereDomain> data_e = {
      first_e, first_e, voronoi, domain, options, vertices, norms, 0, dim, 0};
  int n = sites.size();
  nlopt::opt opt_sites(nlopt::LD_LBFGS, n);
  opt_sites.set_min_objective(&energy_objective, static_cast<void *>(&data_e));

  // optimizaition params
  opt_sites.set_xtol_rel(1e-15);  // adjust
  opt_sites.set_ftol_rel(1e-15);
  // opt_sites.set_maxeval(20);

  // set the lower and upper bounds
  std::vector<double> lower_bound(n, -3.);
  std::vector<double> upper_bound(n, 3.);
  opt_sites.set_lower_bounds(lower_bound);
  opt_sites.set_upper_bounds(upper_bound);

  double e_opt;
  std::vector<double> X = vertices.data();
  nlopt::result result = opt_sites.optimize(X, e_opt);
  LOG << fmt::format(
      "CVT Diagram Energy = {:1.10e} -> {} Iterations to Converge", e_opt,
      data_e.it);
  LOG << fmt::format(
      "Percentage Decrease Relative to Initial Energy = {:2.3f}%",
      ((first_e - e_opt) / first_e) * 100);

  ASSERT(result == nlopt::SUCCESS);
}
UT_TEST_CASE_END(optimized_points_sphere)
*/

UT_TEST_CASE(energy_test_square) {
  double first_e;
  double curr_e;
  double prev_e;
  // creating a voronoi diagram on a square domain (copied from
  // voronoi_test.cpp)
  static const int dim = 3;
  size_t n_sites = 1e3;
  std::vector<coord_t> sites(n_sites * dim, 0.0);
  for (size_t k = 0; k < n_sites; k++) {
    sites[k * dim + 0] = double(rand()) / double(RAND_MAX);
    sites[k * dim + 1] = double(rand()) / double(RAND_MAX);
    sites[k * dim + 2] = 0.0;
  }
  std::vector<index_t> order(n_sites);
  sort_points_on_zcurve(sites.data(), n_sites, dim, order);
  Vertices vertices(dim);
  vertices.reserve(n_sites);
  coord_t x[dim];
  for (size_t i = 0; i < n_sites; i++) {
    for (int d = 0; d < dim; d++) x[d] = sites[dim * order[i] + d];
    vertices.add(x);
  }
  SquareDomain domain;
  VoronoiDiagram voronoi(dim, vertices[0], n_sites);
  VoronoiDiagramOptions options;
  options.n_neighbors = 75;
  options.parallel = true;
  options.store_mesh = true;
  options.verbose = false;
  voronoi.compute(domain, options);
  curr_e = energy(voronoi, 0, 4, true, dim);
  int n_iter = 10;
  for (int iter = 1; iter <= n_iter; ++iter) {
    voronoi.smooth(vertices, false);

    // Examining energy decrase as diagram gets smoothed
    if (iter == 1) {
      LOG << fmt::format("IT. {} -> Initial Energy = {:1.10e}", iter, curr_e);
      first_e = curr_e;
    } else {
      if (prev_e > curr_e) {
        LOG << fmt::format(
            "IT. {} -> Energy = {:1.10e}: Energy Decreased by {:2.3f}%", iter,
            curr_e, ((prev_e - curr_e) / prev_e) * 100);
      }
      if (prev_e <= curr_e && prev_e != 0.) {
        LOG << fmt::format(
            "IT. {} -> Energy = {:1.10e} Energy Increased by {:2.3f}%", iter,
            curr_e, ((curr_e - prev_e) / curr_e) * 100);
      }
      UT_ASSERT(prev_e > curr_e);
    }
    voronoi.vertices().clear();
    voronoi.polygons().clear();
    voronoi.triangles().clear();
    voronoi.compute(domain, options);
    prev_e = curr_e;
    curr_e = energy(voronoi, 0, 4, true, dim);
  }
  LOG << fmt::format(
      "Percentage Decrease Relative to Initial Energy = {:2.3f}%",
      ((first_e - curr_e) / first_e) * 100);
}
UT_TEST_CASE_END(energy_test_square)

/*
UT_TEST_CASE(energy_test_sphere) {
  double first_e;
  double curr_e;
  double prev_e;
  // creating a voronoi diagram on a sphere domain (copied from
  // voronoi_test.cpp)
  auto irand = [](int min, int max) {
    return min + double(rand()) / (double(RAND_MAX) + 1.0) * (max - min);
  };
  static const int dim = 3;
  size_t n_sites = 1e3;
  std::vector<coord_t> sites(n_sites * dim, 0.0);
  for (size_t k = 0; k < n_sites; k++) {
    coord_t theta = 2.0 * M_PI * irand(0, 1);
    coord_t phi = acos(2.0 * irand(0, 1) - 1.0);
    sites[k * dim + 0] = cos(theta) * sin(phi);
    sites[k * dim + 1] = sin(theta) * sin(phi);
    sites[k * dim + 2] = cos(phi);
  }
  std::vector<index_t> order(n_sites);
  sort_points_on_zcurve(sites.data(), n_sites, dim, order);
  Vertices vertices(dim);
  vertices.reserve(n_sites);
  coord_t x[dim];
  for (size_t i = 0; i < n_sites; i++) {
    for (int d = 0; d < dim; d++) x[d] = sites[dim * order[i] + d];
    vertices.add(x);
  }
  SphereDomain domain;
  VoronoiDiagram voronoi(dim, vertices[0], n_sites);
  VoronoiDiagramOptions options;
  options.n_neighbors = 75;
  options.parallel = true;
  options.verbose = false;
  options.store_mesh = true;
  voronoi.compute(domain, options);
  curr_e = energy(voronoi, 1, 4, true, dim);
  int iter = 1;
  while (iter == 1 || ((prev_e - curr_e) / prev_e) * 100 > 0.005) {
    voronoi.smooth(vertices, true);
    // Examining energy decrase as diagram gets smoothed
    if (iter == 1) {
      LOG << fmt::format("IT. {} -> Initial Energy = {:1.10e}", iter, curr_e);
      first_e = curr_e;
    } else {
      if (prev_e > curr_e) {
        LOG << fmt::format(
            "IT. {} -> Energy = {:1.10e}: Energy Decreased by {:2.3f}%", iter,
            curr_e, ((prev_e - curr_e) / prev_e) * 100);
      }
      if (prev_e <= curr_e && prev_e != 0.) {
        LOG << fmt::format(
            "IT. {} -> Energy = {:1.10e} Energy Increased by {:2.3f}%", iter,
            curr_e, ((curr_e - prev_e) / curr_e) * 100);
      }
      UT_ASSERT(prev_e > curr_e);
    }
    voronoi.vertices().clear();
    voronoi.polygons().clear();
    voronoi.triangles().clear();
    voronoi.compute(domain, options);
    prev_e = curr_e;
    curr_e = energy(voronoi, 1, 4, true, dim);

    iter++;
  }
  LOG << fmt::format(
      "Percentage Decrease Relative to Initial Energy = {:2.3f}%",
      ((first_e - curr_e) / first_e) * 100);
}
UT_TEST_CASE_END(energy_test_sphere)
*/

UT_TEST_SUITE_END(smooth_test_suite)