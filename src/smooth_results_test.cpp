/*
Author: Col McDermott

Test file for voronoi diagram energy caclulation and optimized site
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
  LOG << fmt::format("Energy = {:1.10e}", energy);

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
      norm += pow(length(grads), 2);
    }
    norm = sqrt(norm);
    LOG << fmt::format("Grad. Norm = {:1.10e}", norm);
  }

  return energy;
}

// helper method for energy_objective (proecting new site guesses onto the unit
// sphere domain)
void new_sites(Vertices &vertices, const double *guess, int dim,
               unsigned n_sites, bool proj) {
  if (proj) {
    for (auto i = 0; i < n_sites; i++) {
      vec<double, 3> proj(guess[dim * i], guess[dim * i + 1],
                          guess[dim * i + 2]);
      proj = unit_vector(proj);
      for (int j = 0; j < dim; j++) {
        vertices[i][j] = proj[j];
      }
    }
  } else {
    for (auto i = 0; i < n_sites; i++) {
      for (int j = 0; j < dim; j++) {
        vertices[i][j] = guess[(3 * i) + j];
      }
    }
  }
}

// Data for site optimization
template <typename D>
struct energy_data {
  VoronoiDiagram &vor;
  D &domain;
  VoronoiDiagramOptions &options;
  Vertices &verts;  // sites
  std::vector<double> &norms;
  double norm;
  int dim;
  int it{0};
};

// optimization constraint
struct sphere_data {
  double a, b, c, r;
};

// Sphere constraint
double sphere_constr(unsigned n, const double *X, double *grad, void *sd) {
  // data type conversion
  sphere_data &data_s = *static_cast<sphere_data *>(sd);

  // if constraint gradient is needed
  if (grad) {
    grad[0] = 2 * (X[0] - data_s.a);
    grad[1] = 2 * (X[1] - data_s.b);
    grad[2] = 2 * (X[2] - data_s.c);
  }

  // return vale of constraint
  return pow(X[0] - data_s.a, 2) + pow(X[1] - data_s.b, 2) +
         pow(X[2] - data_s.c, 2) - pow(data_s.r, 2);
}

// Objective function to minimize
double energy_objective(unsigned n, const double *sites, double *grad,
                        void *ed) {
  // get data
  energy_data<SphereDomain> &data_e =
      *static_cast<energy_data<SphereDomain> *>(ed);

  // increase iteration
  data_e.it++;
  LOG << fmt::format("IT. {}", data_e.it);

  // reset sites and recompute diagram
  data_e.vor.vertices().clear();
  new_sites(data_e.verts, sites, data_e.dim, n / data_e.dim, false);
  data_e.vor.polygons().clear();
  data_e.vor.triangles().clear();
  data_e.vor.compute(data_e.domain, data_e.options);

  // calculate energy and record current guess
  double curr_e = energy(data_e.vor, 1, 10, false, data_e.dim);

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

  return curr_e;
}

UT_TEST_SUITE(smooth_test_suite)

UT_TEST_CASE(convergence_comparison) {
  /*** adjustments ***/
  size_t n_sites = 1e3;
  static const int dim = 3;
  bool pre = true;
  std::string alg;
  double xtol = 1e-12;
  double ftol = 1e-12;
  int eval = 100;
  /********************/

  auto irand = [](int min, int max) {
    return min + double(rand()) / (double(RAND_MAX) + 1.0) * (max - min);
  };
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
  Vertices vertices1(dim);
  Vertices vertices2(dim);
  vertices1.reserve(n_sites);
  vertices2.reserve(n_sites);
  coord_t x[dim];
  for (size_t i = 0; i < n_sites; i++) {
    for (int d = 0; d < dim; d++) x[d] = sites[dim * order[i] + d];
    vertices1.add(x);
    vertices2.add(x);
  }
  SphereDomain domain;
  VoronoiDiagram voronoi1(dim, vertices1[0], n_sites);
  VoronoiDiagramOptions options;
  options.n_neighbors = 75;
  options.parallel = true;
  options.verbose = false;
  options.store_mesh = true;
  voronoi1.compute(domain, options);

  // pre-processing with Lloyd Relaxation
  if (pre) {
    for (int i = 0; i < 50; i++) {
      voronoi1.smooth(vertices1, true);
      voronoi1.vertices().clear();
      voronoi1.polygons().clear();
      voronoi1.triangles().clear();
      voronoi1.compute(domain, options);
    }
  }

  std::cout << "Optimization Method\n" << std::endl;

  double e_opt = energy(voronoi1, 1, 10, true, dim);
  std::vector<double> norms(n_sites, 0.);

  // fill optimization data
  energy_data<SphereDomain> data_e = {voronoi1, domain, options, vertices1,
                                      norms,    0,      dim,     0};
  int n = sites.size();
  nlopt::opt opt_sites(nlopt::LD_LBFGS, n);  // change alg.
  opt_sites.set_min_objective(&energy_objective, static_cast<void *>(&data_e));

  alg = opt_sites.get_algorithm_name();
  if (alg.compare("LD_SLSQP") == 0) {
    sphere_data data_s = {0., 0., 0., 1.};
    for (int i = 0; i < n_sites; i++) {
      opt_sites.add_equality_constraint(&sphere_constr, &data_s, 1e-8);
    }
  }

  // optimizaition params
  opt_sites.set_xtol_rel(xtol);
  opt_sites.set_ftol_rel(ftol);

  opt_sites.set_maxeval(eval);

  // set the lower and upper bounds
  std::vector<double> lower_bound(n, -3.);
  std::vector<double> upper_bound(n, 3.);
  opt_sites.set_lower_bounds(lower_bound);
  opt_sites.set_upper_bounds(upper_bound);

  e_opt = 0;
  std::vector<double> X = vertices1.data();
  nlopt::result result = opt_sites.optimize(X, e_opt);
  // UT_ASSERT_EQUALS(result, nlopt::SUCCESS);

  std::cout << "\n\nLloyd Relaxation\n" << std::endl;

  VoronoiDiagram voronoi2(dim, vertices2[0], n_sites);
  voronoi2.compute(domain, options);
  double e_lloyd = energy(voronoi2, 1, 10, true, dim);
  for (int i = 1; i < eval + 1; i++) {
    voronoi2.smooth(vertices2, true);
    voronoi2.vertices().clear();
    voronoi2.polygons().clear();
    voronoi2.triangles().clear();
    voronoi2.compute(domain, options);
    LOG << fmt::format("IT. {}", i);
    e_lloyd = energy(voronoi2, 1, 10, true, dim);
  }
}
UT_TEST_CASE_END(convergence_comparison)

UT_TEST_SUITE_END(smooth_test_suite)