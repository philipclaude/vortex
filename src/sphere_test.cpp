/*
Author: Col McDermott

Test file for voronoi diagram energy caclulation and optimized site
smoothing using spherical coordinates
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

// Data for site optimization
template <typename D>
struct energy_data {
  VoronoiDiagram &vor;
  D &domain;
  VoronoiDiagramOptions &options;
  Vertices &sites;
  std::vector<double> &norms;
  double first_e;
  double percent_change;
  double norm;
  int dim;
  int it{0};
};

// helper method for energy_objective (mapping sites from spherical coords to
// Cartesian coords)
void new_sites(Vertices &sites, const double *guess, int dim,
               unsigned n_sites) {
  double x, y, z;
  for (auto i = 0; i < n_sites; i++) {
    x = cos(guess[(dim * i)]) * sin(guess[(dim * i) + 1]);
    y = sin(guess[(dim * i)]) * sin(guess[(dim * i) + 1]);
    z = cos(guess[(dim * i) + 1]);

    sites[i][0] = x;
    sites[i][1] = y;
    sites[i][2] = z;
  }
}

double energy(VoronoiDiagram &voronoi, int t_type, int order, bool grad) {
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
  LOG << fmt::format("Energy = {:1.10e}", energy);

  if (grad) {
    double norm = 0;
    VoronoiCellProperties props;
    for (int i = 0; i < polygons.n(); i++) {  // might need to hard-coded but
                                              // should be polygons.n();
      props = voronoi.properties()[i];
      vec3 site(voronoi.vertices()[props.site][0],
                voronoi.vertices()[props.site][1],
                voronoi.vertices()[props.site][2]);
      vec3 grads(0., 0., 0.);
      grads = (2 * ((props.volume * site) - props.moment));
      norm += pow(length(grads), 2);
    }

    norm = sqrt(norm);
    LOG << fmt::format("Grad. Norm = {:1.10e}", norm);
  }

  return energy;
}

// Objective function to minimize
double energy_objective(unsigned n, const double *guess, double *grad,
                        void *ed) {
  // get data
  energy_data<SphereDomain> &data_e =
      *static_cast<energy_data<SphereDomain> *>(ed);

  // increase iteration
  data_e.it++;
  LOG << fmt::format("IT. {}", data_e.it);

  // reset sites and recompute diagram
  data_e.vor.vertices().clear();
  data_e.vor.polygons().clear();
  data_e.vor.triangles().clear();
  new_sites(data_e.sites, guess, data_e.dim, n / data_e.dim);
  data_e.vor.compute(data_e.domain, data_e.options);

  // calculate energy (record current guess)
  double curr_e = energy(data_e.vor, 1, 10, false);

  // calculate gradient values
  VoronoiCellProperties props;
  if (grad) {  // change
    data_e.norm = 0;
    for (int i = 0; i < (n / data_e.dim); i++) {
      props = data_e.vor.properties()[i];
      vec3 site(data_e.sites[props.site][0], data_e.sites[props.site][1],
                data_e.sites[props.site][2]);
      vec3 grad_i(0., 0., 0.);
      grad_i = (2 * ((props.volume * site) - props.moment));

      grad[(data_e.dim * i)] = (grad_i[1] * site[0]) - (grad_i[0] * site[1]);
      grad[(data_e.dim * i) + 1] =
          (grad_i[0] *
           (cos(guess[data_e.dim * i]) * cos(guess[data_e.dim * i + 1]))) +
          (grad_i[1] *
           (sin(guess[data_e.dim * i]) * cos(guess[data_e.dim * i + 1]))) -
          (grad_i[2] * sin(guess[data_e.dim * i + 1]));

      data_e.norm +=
          pow(grad[data_e.dim * i], 2) + pow(grad[data_e.dim * i + 1], 2);
    }
    data_e.norm = sqrt(data_e.norm);
    data_e.norms.push_back(data_e.norm);
    LOG << fmt::format("Grad. Norm = {:1.10e}",
                       data_e.norms[data_e.norms.size() - 1]);
  }

  data_e.percent_change = ((data_e.first_e - curr_e) / data_e.first_e) * 100;

  LOG << fmt::format("Percent Decrease in Energy = {:2.3f}%",
                     ((data_e.first_e - curr_e) / data_e.first_e) * 100);

  if (data_e.percent_change > 30) {
    return -1 * INFINITY;
  } else {
    return curr_e;
  }
}

UT_TEST_SUITE(smooth_test_suite)

UT_TEST_CASE(mapping_test) {
  /*** adjustments ***/
  size_t n_sites = 1e4;
  const int dim = 2;
  bool pre = false;
  std::string str = "\n";
  if (pre) str = "w/ Pre-Processing (50 iter. of Lloyd Relaxation)\n";
  double xtol = 1e-12;
  double ftol = 1e-12;
  int eval = 100;
  /********************/

  SphereDomain domain;
  VoronoiDiagramOptions options;
  options.n_neighbors = 75;
  options.parallel = true;
  options.verbose = false;
  options.store_mesh = true;

  auto irand = [](int min, int max) {
    return min + double(rand()) / (double(RAND_MAX) + 1.0) * (max - min);
  };
  std::vector<coord_t> sites(n_sites * dim, 0.);
  std::vector<coord_t> sites_cart(n_sites * 3, 0.);

  for (size_t k = 0; k < n_sites; k++) {
    coord_t theta = 2.0 * M_PI * irand(0, 1);
    coord_t phi = acos(2.0 * irand(0, 1) - 1.0);

    sites[k * dim] = theta;
    sites[k * dim + 1] = phi;

    sites_cart[k * 3] = cos(theta) * sin(phi);
    sites_cart[(k * 3) + 1] = sin(theta) * sin(phi);
    sites_cart[(k * 3) + 2] = cos(phi);
  }
  std::vector<index_t> order(n_sites);
  sort_points_on_zcurve(sites_cart.data(), n_sites, 3, order);

  // for debugging
  for (int i = 0; i < n_sites; i++) {
    order[i] = i;
  }

  Vertices vertices_o(3);
  Vertices vertices_l(3);
  vertices_o.reserve(n_sites);
  vertices_l.reserve(n_sites);
  coord_t x[3];
  coord_t y[3];
  for (size_t i = 0; i < n_sites; i++) {
    for (int d = 0; d < 3; d++) {
      x[d] = sites_cart[3 * order[i] + d];
      if (d < 2) {
        y[d] = sites[dim * order[i] + d];
      } else {
        y[d] = 0.;
      }
    }
    vertices_l.add(x);
    vertices_o.add(y);
  }

  std::cout << "\n\nLloyd Relaxation\n" << std::endl;

  VoronoiDiagram voronoi_l(3, vertices_l[0], n_sites);
  voronoi_l.compute(domain, options);
  double e_lloyd = energy(voronoi_l, 1, 10, true);
  for (int i = 1; i < eval + 1; i++) {
    voronoi_l.smooth(vertices_l, true);
    voronoi_l.vertices().clear();
    voronoi_l.polygons().clear();
    voronoi_l.triangles().clear();
    voronoi_l.compute(domain, options);
    LOG << fmt::format("IT. {}", i);
    e_lloyd = energy(voronoi_l, 1, 10, true);
  }

  std::cout << "\n\nOptimization Method " << str << std::endl;

  auto map = [sites, &vertices_o, n_sites, dim]() {
    double x, y, z;
    for (int i = 0; i < n_sites; i++) {
      x = cos(sites[(dim * i)]) * sin(sites[(dim * i) + 1]);
      y = sin(sites[(dim * i)]) * sin(sites[(dim * i) + 1]);
      z = cos(sites[(dim * i) + 1]);

      vertices_o[i][0] = x;
      vertices_o[i][1] = y;
      vertices_o[i][2] = z;
    }
  };
  map();

  VoronoiDiagram voronoi_o(3, vertices_o[0], n_sites);
  voronoi_o.compute(domain, options);
  // pre-processing with Lloyd Relaxation
  if (pre) {
    for (int i = 0; i < 50; i++) {
      voronoi_o.smooth(vertices_o, true);
      voronoi_o.vertices().clear();
      voronoi_o.polygons().clear();
      voronoi_o.triangles().clear();
      voronoi_o.compute(domain, options);
    }

    for (int i = 0; i < n_sites; i++) {
      sites[dim * i] = atan2(vertices_o[i][1], vertices_o[i][0]);
      sites[dim * i + 1] = acos(vertices_o[i][2] / 1);
    }
  }
  double e_opt = energy(voronoi_o, 1, 10, true);

  std::vector<double> norms(n_sites, 0.);

  // fill optimization data
  energy_data<SphereDomain> data_e = {
      voronoi_o, domain, options, vertices_o, norms, e_opt, 0, 0, dim, 0};
  int n = n_sites * dim;

  nlopt::opt opt_sites(nlopt::LD_LBFGS, n);
  opt_sites.set_min_objective(&energy_objective, static_cast<void *>(&data_e));

  // optimizaition params
  opt_sites.set_xtol_rel(xtol);
  opt_sites.set_ftol_rel(ftol);
  opt_sites.set_maxeval(eval);
  opt_sites.set_stopval(-1 * INFINITY);

  nlopt::result result = opt_sites.optimize(sites, e_opt);
  UT_ASSERT_EQUALS(result, nlopt::SUCCESS);
}
UT_TEST_CASE_END(mapping_test)

UT_TEST_SUITE_END(smooth_test_suite)