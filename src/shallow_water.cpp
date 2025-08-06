//
//  vortex: Voronoi mesher and fluid simulator for the Earth's oceans and
//  atmosphere.
//
//  Copyright 2023 - 2025 Philip Claude Caplan
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
#include "shallow_water.h"

#include <argparse/argparse.hpp>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>

#include "io.h"
#include "library.h"
#include "math/spmat.h"
#include "math/vec.hpp"
#include "operators.h"

namespace vortex {

template <typename Domain_t>
void ShallowWaterSimulation<Domain_t>::setup() {
  for (size_t k = 0; k < particles_.n(); k++) {
    double h0 = options_.initial_height(particles_[k]);
    if (h0 >= 0) {
      // negative height can be a signal from the case that means
      // the height was imported from a file
      height_[k] = h0 - options_.surface_height(particles_[k]);
    }
  }
  particles_.set_velocity(options_.initial_velocity);
  particles_.set_density([](const double* x) { return 1.0; });
  calculate_initial_conservation_quantities();
}

template <typename Domain_t>
void ShallowWaterSimulation<
    Domain_t>::calculate_initial_conservation_quantities() {
  for (size_t k = 0; k < particles_.n(); k++) {
    volume_[k] = height_[k] * voronoi_.properties()[k].volume;
  }
  initial_mass_ = total_mass();
  initial_momentum_ = total_momentum();
  initial_energy_ = total_energy();
}

template <typename Domain_t>
double ShallowWaterSimulation<Domain_t>::time_step(
    const SimulationOptions& options) {
  reset_timers();
  Timer time_step_timer;
  time_step_timer.start();
  calculate_properties();  // mass, volume, centroids, max displacement
  auto& velocity = particles_.velocity();
  const size_t n = particles_.n();
  int dim = particles_.dim();

  // differential operators are calculated on a unit sphere, so they should be
  // scaled by 1/a since (d/dx)_earth = (d/dx)_unit_sphere / earth_radius
  const double a = earth_.radius;
  const double g = earth_.gravity;  // to calculate gravity wave speeds

  // save target area
  std::vector<double> area(n, 0.0);
  for (size_t k = 0; k < n; k++) {
    area[k] = volume_[k] / height_[k];
  }

  // from the latest Voronoi diagram, determine the maximum time step that
  // can be taken so points always remain in the current cell
  double dt = options.time_step;
  for (size_t k = 0; k < n; k++) {
    vec3d u(velocity[k]);
    double um = length(u);
    // double cm = std::sqrt(g * height_[k]);
    //  TODO this should be done from the advection point (site or centroid)
    //  currently max displacement is measured from the site
    double dt_k = a * max_displacement_[k] / um;  //(um + cm);
    if (dt_k < dt) dt = dt_k;
  }

  // find a time step such that we at least have a Voronoi diagram
  std::vector<double> position(n * dim);
  for (size_t k = 0; k < n; k++) {
    for (int d = 0; d < 3; d++) position[3 * k + d] = particles_(k, d);
  }
  while (true) {
    for (size_t k = 0; k < n; k++) {
      voronoi_.weights()[k] = 0.0;
      vec3d x(particles_[k]);
      vec3d u(velocity[k]);
      vec3d dx = dt * u / a;  // displacement on unit-sphere

      for (int d = 0; d < 3; d++) {
        if (options.advect_from_centroid)
          particles_(k, d) = particles_.centroids()(k, d) + dx[d];
        else
          particles_(k, d) = particles_(k, d) + dx[d];
      }
    }

    // calculate voronoi diagram
    VoronoiDiagramOptions voro_opts;
    voro_opts.verbose = false;
    calculate_power_diagram(domain_, voro_opts);

    size_t n_vanish = 0;
    for (size_t k = 0; k < n; k++) {
      if (voronoi_.properties()[k].volume <= 0) {
        n_vanish++;
      }
    }
    if (n_vanish == 0) break;
    dt *= 0.5;
  }
  // restore particle position
  for (size_t k = 0; k < n; k++) {
    for (int d = 0; d < 3; d++) {
      particles_(k, d) = position[3 * k + d];
    }
  }

  // limit time step so height is always positive
  // TODO is this necessary with semi-implicit time stepping?
  VoronoiOperators<Domain_t> ops(voronoi_);
  std::vector<double> div_u(n);
  ops.calculate_divergence(velocity[0], div_u.data());
  while (true) {
    bool ok = true;
    for (size_t k = 0; k < n; k++) {
      double h = height_[k];
      if (h - dt * h * div_u[k] / a < 1e-6) {
        ok = false;
        break;
      }
    }
    if (ok) break;
    dt = 0.5 * dt;
  }

  // calculate gradient of original total height before advection
  ops.set_boundary_value(0.0);
  std::vector<double> grad_h(3 * n);
  std::vector<double> height(n, 0);
  std::vector<double> hs(n), grad_hs(3 * n);
  for (size_t k = 0; k < n; k++) {
    height[k] = height_[k] + options_.surface_height(particles_[k]);
    hs[k] = options_.surface_height(particles_[k]);
  }
  ops.calculate_gradient(height.data(), grad_h.data());
  ops.calculate_gradient(hs.data(), grad_hs.data());

  // advance the height field according to governing differential equation
  double linear_solver_time = 0;
  if (options_.time_stepping == TimeSteppingScheme::kExplicit) {
    for (size_t k = 0; k < n; k++) {
      ASSERT(height_[k] > 0);
      height_[k] -= dt * height_[k] * div_u[k] / a;
    }
  } else if (options_.time_stepping == TimeSteppingScheme::kSemiImplicit) {
    // build system to solve for h^{n+1}
    std::vector<double> force(3 * n);
    for (size_t i = 0; i < n; i++) {
      vec3d x(particles_[i]);  // on unit sphere
      vec3d c(particles_.centroids()[i]);
      vec3d u(velocity[i]);
      double f = options_.coriolis_parameter(&c[0]);
      vec3d fc = f * cross(c, u);  // coriolis force
      if (options_.constrain) fc = fc + dot(u, u) * c / a;
      for (int d = 0; d < 3; d++) {
        force[3 * i + d] = u[d] - dt * fc[d] - g * dt * grad_hs[3 * i + d] / a;
      }
    }
    std::vector<double> div_f(n);
    ops.calculate_divergence(force.data(), div_f.data());

    Timer timer;
    timer.start();
    h_mat_.clear();
    for (size_t i = 0; i < n; i++) {
      h_mat_(i, i) = 1.0;
      h_rhs_(i) = height_[i] * (1 - dt * div_f[i] / a) / a;
    }
    for (const auto& facet : voronoi_.facets()) {
      if (facet.bj < 0) continue;
      size_t i = facet.bi;
      size_t j = facet.bj;
      if (i >= particles_.n()) continue;
      if (j >= particles_.n()) continue;
      vec3d ri(particles_[i]);
      vec3d rj(particles_[j]);
      vec3d rij = ri - rj;
      const double lij = facet.length;
      const double wi = voronoi_.properties()[i].volume;
      const double wj = voronoi_.properties()[j].volume;
      double xij = std::sqrt(dot(rij, rij));
      double dij = -height_[i] * g * dt * dt * lij / xij / a / a;
      h_mat_(i, j) = dij / wi;
      h_mat_(j, i) = dij / wj;
      h_mat_(i, i) -= h_mat_(i, j);
      h_mat_(j, j) -= h_mat_(j, i);
    }
    SparseSolverOptions opts;
    // opts.tol = 1e-3;
    opts.symmetric = false;
    opts.max_iterations = 50;
    vecd<double> dh(n);
    h_mat_.solve_nl(h_rhs_, h_sol_, opts);
    timer.stop();
    linear_solver_time = timer.seconds();
    for (size_t i = 0; i < n; i++) {
      height_[i] = h_sol_[i];
    }
  }

  if (options_.use_analytic_height) {
    for (size_t k = 0; k < n; k++)
      height_[k] = options_.analytic_height(particles_[k], options.time);
  }

  // update particle positions using the time step and current velocity
  for (size_t k = 0; k < n; k++) {
    vec3d x(particles_[k]);
    vec3d u(velocity[k]);
    vec3d dx = dt * u / a;  // displacement on unit-sphere

    for (int d = 0; d < 3; d++) {
      if (options.advect_from_centroid)
        particles_(k, d) = particles_.centroids()(k, d) + dx[d];
      else
        particles_(k, d) = particles_(k, d) + dx[d];
    }
    if (options_.project_points) project_point<Domain_t>(particles_[k]);
  }

  SimulationConvergence convergence;
  if (options_.use_optimal_transport) {
    // determine the area we need to conserve mass (volume for incompressible)
    double at = 0.0;
    for (size_t k = 0; k < n; k++) {
      area[k] = volume_[k] / height_[k];
      ASSERT(area[k] > 0) << area[k];
      at += area[k];
    }

    // renormalize the areas to make the optimal transport solution tractable
    for (size_t k = 0; k < n; k++) {
      area[k] *= domain_.area() / at;
    }

    // solve semi-discrete optimal transport problem to reproduce cell areas
    convergence = optimize_volumes(domain_, options, area);
    calculate_properties();  // mass, volume, centroids, max displacement

    // update height to satisfy conservation of mass
    for (size_t i = 0; i < n; i++) {
      height_[i] = volume_[i] / voronoi_.properties()[i].volume;
    }
    for (size_t k = 0; k < n; k++) {
      height[k] = height_[k] + options_.surface_height(particles_[k]);
    }
    ops.calculate_gradient(height.data(), grad_h.data());
  } else {
    VoronoiDiagramOptions voro_opts;
    voro_opts.verbose = false;
    calculate_power_diagram(domain_, voro_opts);
    calculate_properties();  // mass, volume, centroids, max displacement

    // update height to satisfy conservation of mass
    for (size_t i = 0; i < n; i++) {
      ASSERT(voronoi_.properties()[i].site == i);
      height_[i] = volume_[i] / voronoi_.properties()[i].volume;
    }
  }

  // compute artificial viscosity
  std::vector<double> viscous(3 * n, 0.0);
  if (options_.add_artificial_viscosity) compute_artificial_viscosity(viscous);
  if (options_.stabilize_pressure_gradient)
    stabilize_pressure_gradient(height, grad_h);

  // update velocity
  if (!options_.use_analytic_velocity) {
    for (size_t i = 0; i < n; i++) {
      vec3d x(particles_[i]);  // on unit sphere
      vec3d c(particles_.centroids()[i]);
      vec3d u(velocity[i]);
      double f = options_.coriolis_parameter(&c[0]);
      vec3d force = f * cross(c, u);  // coriolis force
      if (options_.constrain) force = force + dot(u, u) * c / a;
      vec3d fv(&viscous[3 * i]);  // derivative needs to be scaled by 1/a
      force = force + fv / a;

      for (int d = 0; d < 3; d++) {
        velocity(i, d) -= dt * (g * grad_h[3 * i + d] / a + force[d]);
      }
      if (options_.project_velocity)
        project_velocity<Domain_t>(particles_[i], velocity[i]);
    }
  } else {
    // use the analytic velocity, e.g. for Williamson Case 1
    for (size_t i = 0; i < n; i++) {
      vec3d u = options_.analytic_velocity(particles_[i], options.time);
      for (int d = 0; d < 3; d++) {
        velocity(i, d) = u[d];
      }
    }
  }

  // calculate error
  double h_error = -1;
  double h_total = -1;
  if (options_.has_analytic_height) {
    h_error = 0.0;
    h_total = 0.0;
    for (size_t k = 0; k < n; k++) {
      double ha = options_.analytic_height(particles_[k], options.time + dt);
      double dh = ha - height_[k];
      h_error += dh * dh * voronoi_.properties()[k].volume;
      h_total += ha * ha * voronoi_.properties()[k].volume;
    }
    h_error = std::sqrt(h_error);
    h_total = std::sqrt(h_total);
  }
  double u_error = -1;
  double u_total = -1;
  if (options_.has_analytic_velocity) {
    u_error = 0.0;
    u_total = 0.0;
    for (size_t k = 0; k < n; k++) {
      vec3d ua = options_.analytic_velocity(particles_[k], options.time + dt);
      vec3d uk(particles_.velocity()[k]);
      double du = dot(ua - uk, ua - uk);
      u_error += du * du * voronoi_.properties()[k].volume;
      u_total += dot(ua, ua) * voronoi_.properties()[k].volume;
    }
    u_error = std::sqrt(u_error);
    u_total = std::sqrt(u_total);
  }

  // monitor the simulation
  if (options.iteration % options.print_frequency == 0) print_header();
  double sdpd = options.time / timer_.seconds();
  double area_error = (domain_.area() - total_area());
  double mass_error = (total_mass() - initial_mass_) / initial_mass_;
  double momentum_error =
      (total_momentum() - initial_momentum_) / initial_momentum_;
  double energy_error = (total_energy() - initial_energy_) / initial_energy_;
  std::cout << fmt::format(
      "| {:6d} | {:9s} | {:1.1e} | {:2d}:{:1.1e} | {:+1.1e} | {:+1.1e} | "
      "{:+1.1e} "
      "| {:+1.1e} | {:6.1f} | {:+1.1e} | {:+1.1e} |\n",
      options.iteration, days_hours_minutes(options.time + dt), dt,
      convergence.n_iterations, convergence.error, area_error, mass_error,
      momentum_error, energy_error, sdpd, h_error / h_total, u_error / u_total);
  time_step_timer.stop();

  statistics_.ra.push_back(area_error);
  statistics_.rw.push_back(convergence.error);
  statistics_.rm.push_back(mass_error);
  statistics_.rp.push_back(momentum_error);
  statistics_.re.push_back(energy_error);
  statistics_.sdpd.push_back(sdpd);
  statistics_.h_error.push_back(h_error);
  statistics_.h_total.push_back(h_total);
  statistics_.u_error.push_back(u_error);
  statistics_.u_total.push_back(u_total);
  statistics_.time.push_back(options.time + dt);
  statistics_.voronoi_time.push_back(voronoi_time_);
  statistics_.n_voronoi.push_back(n_voronoi_);
  statistics_.linear_solver_time.push_back(linear_solver_time_ +
                                           linear_solver_time);
  statistics_.time_step_time.push_back(time_step_timer.seconds());

  return dt;
}

template <typename Domain_t>
void ShallowWaterSimulation<Domain_t>::compute_artificial_viscosity(
    std::vector<double>& fv) {
  const size_t n = particles_.n();
  double g = earth_.gravity;
  // double a = earth_.radius;
  double eps = 1e-6;

  fv.resize(3 * n, 0);
  ASSERT(voronoi_.facets().size() > 0);
  for (const auto& facet : voronoi_.facets()) {
    if (facet.bj < 0) continue;
    size_t i = facet.bi;
    size_t j = facet.bj;
    if (i >= particles_.n()) continue;
    if (j >= particles_.n()) continue;
    vec3d ri(particles_[i]);
    vec3d rj(particles_[j]);
    vec3d rij = ri - rj;
    double xij = std::sqrt(dot(rij, rij));

    vec3d ui(particles_.velocity()[i]);
    vec3d uj(particles_.velocity()[j]);
    vec3d uij = ui - uj;

    double hi = height_[i];
    double hj = height_[j];

    double ci = std::sqrt(g * hi);
    double cj = std::sqrt(g * hj);
    double cij = 0.5 * (ci + cj);

    double pi_ij = -cij * dot(uij, rij) / std::sqrt(dot(rij, rij) + eps * eps);

    const vec3d mij(&facet.midpoint.x);
    const double lij = facet.length;
    double wi = voronoi_.properties()[i].volume;
    double wj = voronoi_.properties()[j].volume;

    for (int d = 0; d < 3; d++) {
      fv[3 * i + d] += -lij * (ri[d] - mij[d]) * pi_ij / (xij * wi);
      fv[3 * j + d] += -lij * (rj[d] - mij[d]) * pi_ij / (xij * wj);
    }
  }
}

template <typename Domain_t>
void ShallowWaterSimulation<Domain_t>::stabilize_pressure_gradient(
    const std::vector<double>& h, std::vector<double>& dh) {
  const size_t n = particles_.n();

  // compute laplacian of height field
  std::vector<double> lh(n, 0);
  ASSERT(voronoi_.facets().size() > 0);
  for (const auto& facet : voronoi_.facets()) {
    if (facet.bj < 0) continue;
    size_t i = facet.bi;
    size_t j = facet.bj;
    if (i >= particles_.n()) continue;
    if (j >= particles_.n()) continue;
    vec3d ri(particles_[i]);
    vec3d rj(particles_[j]);
    vec3d rij = ri - rj;
    const double lij = facet.length;
    double xij = std::sqrt(dot(rij, rij));

    double hi = height_[i];
    double hj = height_[j];

    double wi = voronoi_.properties()[i].volume;
    double wj = voronoi_.properties()[j].volume;

    double hij = hi - hj;
    lh[i] -= lij * hij / (xij * wi);
    lh[j] += lij * hij / (xij * wj);
  }

  // stabilize the pressure gradient
  const double d = 2;
  const double fd = 0.5 * (d + 1) / d;
  for (size_t k = 0; k < n; k++) {
    if (lh[k] < 0) continue;
    vec3d xi(particles_[k]);
    vec3d ci(particles_.centroids()[k]);
    vec3d dhs = fd * lh[k] * (ci - xi);
    for (int j = 0; j < 3; j++) {
      dh[3 * k + j] -= dhs[j];
    }
  }
}

template <typename Domain_t>
void ShallowWaterSimulation<Domain_t>::print_header(int n_bars) const {
  std::cout << fmt::format("{:->{}}", "", n_bars) << std::endl;
  std::cout << fmt::format(
      "| {:6s} | {:9s} | {:7s} | {:10s} | {:8s} | {:8s} | {:8s} | {:8s} | "
      "{:6s} | {:8s} | {:8s} |\n",
      "Step", "day:hr:mn", "dt (s)", "Rw", "Ra", "Rm", "Rp", "Re", "SDPD", "Eh",
      "Eu");
  std::cout << fmt::format("{:->{}}", "", n_bars) << std::endl;
}

template <typename Domain_t>
void ShallowWaterSimulation<Domain_t>::save(const std::string& filename) const {
  const auto a = earth_.radius;
  size_t n = particles_.n();
  std::vector<double> w(3 * n, 0.0);
  VoronoiOperators<Domain_t> ops(voronoi_);
  ops.set_boundary_value(0.0);
  ops.calculate_curl(particles_.velocity()[0], w.data());

  size_t n_data = n * 3;
  std::vector<float> data(n_data);
  size_t m = 0;
  for (size_t k = 0; k < n; k++) {
    const double* x = particles_[k];
    if (options_.save_latlon) {
      double theta = atan2(x[1], x[0]);
      double lambda = asin(x[2]);
      data[m++] = theta;
      data[m++] = lambda;
      data[m++] = 0;
    } else {
      for (int d = 0; d < 3; d++) {
        data[m++] = x[d];
      }
    }
  }
  FILE* fid = fopen(filename.c_str(), "wb");
  fprintf(fid, "# vtk DataFile Version 2.0\nvortex vertices\n");
  fprintf(fid, "BINARY\nDATASET UNSTRUCTURED_GRID\nPOINTS %zu float\n", n);
  std::parafor_i(0, n_data,
                 [&data](int tid, size_t k) { io::swap_end(data[k]); });
  fwrite(&data[0], 4, n_data, fid);
  fprintf(fid, "\nCELLS 0 0\nCELL_TYPES 0\n");
  fprintf(fid, "POINT_DATA %zu\n", n);
  fprintf(fid, "FIELD FieldData 3\n");
  std::vector<float> height_data(n), pv_data(n), rv_data(n);
  for (size_t k = 0; k < n; k++) {
    vec3d wk(w.data() + 3 * k);
    vec3d n(particles_[k]);
    double rv = dot(n, wk);
    height_data[k] = height_[k] + options_.surface_height(particles_[k]);
    rv_data[k] = rv / a;
    pv_data[k] =
        (rv / a + options_.coriolis_parameter(particles_[k])) / height_[k];
  }
  std::parafor_i(0, n, [&](int tid, size_t k) {
    io::swap_end(height_data[k]);
    io::swap_end(pv_data[k]);
    io::swap_end(rv_data[k]);
  });

  fprintf(fid, "Height 1 %zu float\n", n);
  fwrite(&height_data[0], 4, n, fid);
  fprintf(fid, "RelativeVorticity 1 %zu float\n", n);
  fwrite(&rv_data[0], 4, n, fid);
  fprintf(fid, "PotentialVorticity 1 %zu float\n", n);
  fwrite(&pv_data[0], 4, n, fid);
  fprintf(fid, "\n");
  fclose(fid);
}

template <typename Domain_t>
void ShallowWaterSimulation<Domain_t>::save_json(
    const std::string& filename) const {
  nlohmann::json data;

  size_t n = particles_.n();
  std::vector<double> x(n), y(n), z(n), h(n);
  for (size_t k = 0; k < n; k++) {
    x[k] = particles_[k][0];
    y[k] = particles_[k][1];
    z[k] = particles_[k][2];
    h[k] = height_[k] + options_.surface_height(particles_[k]);
  }
  data["x"] = x;
  data["y"] = y;
  data["z"] = z;
  data["h"] = h;
  data["w"] = voronoi_.weights();
  data["domain"] = "sphere";

  std::ofstream outfile(filename);
  outfile << std::setw(4) << data << std::endl;
}

template <typename Domain_t>
double ShallowWaterSimulation<Domain_t>::total_area() const {
  double area = 0.0;
  for (size_t i = 0; i < particles_.n(); i++) {
    area += voronoi_.properties()[i].volume;
  }
  return area;
}

template <typename Domain_t>
double ShallowWaterSimulation<Domain_t>::total_mass() const {
  double mass = 0;
  for (size_t k = 0; k < particles_.n(); k++) {
    mass += voronoi_.properties()[k].volume * height_[k];
  }
  return mass;
}

template <typename Domain_t>
double ShallowWaterSimulation<Domain_t>::total_momentum() const {
  double g = earth_.gravity;
  double momentum = 0;
  for (size_t k = 0; k < particles_.n(); k++) {
    double vi = voronoi_.properties()[k].volume;
    vec3d u(particles_.velocity()[k]);
    double hi = height_[k];
    momentum += vi * (length(u) + std::sqrt(g * hi));
  }
  return momentum;
}

template <typename Domain_t>
double ShallowWaterSimulation<Domain_t>::total_energy() const {
  double g = earth_.gravity;
  double energy = 0;
  for (size_t k = 0; k < particles_.n(); k++) {
    double vi = voronoi_.properties()[k].volume;
    vec3d u(particles_.velocity()[k]);
    double hi = height_[k];
    double hs = options_.surface_height(particles_[k]);
    energy += 0.5 * vi * hi * dot(u, u);
    energy += 0.5 * vi * g * ((hi + hs) * (hi + hs) - hs * hs);
  }
  return energy;
}

void run_swe_simulation(const argparse::ArgumentParser& program) {
  static const int dim = 3;
  typedef SphereDomain Domain_t;
  Domain_t domain;

  // create output directory
  std::string output_dir = program.get<std::string>("--output");
  if (output_dir.empty()) {
    output_dir = fmt::format("vortex{}", getpid());
  } else {
    size_t n_removed = std::filesystem::remove_all(output_dir);
    LOG << fmt::format("removed {} files from {}", n_removed, output_dir);
  }
  std::filesystem::create_directories(output_dir);
  std::string prefix = output_dir + "/particles";

  std::string import_height_from =
      program.get<std::string>("--import_height_from");

  // set up the test case
  std::shared_ptr<ShallowWaterOptions> test_case_ptr = nullptr;
  std::string test_case = program.get<std::string>("--case");
  if (test_case == "williamson1") {
    test_case_ptr = std::make_shared<WilliamsonCase1>();
  } else if (test_case == "williamson2") {
    test_case_ptr = std::make_shared<WilliamsonCase2>();
  } else if (test_case == "williamson5") {
    test_case_ptr = std::make_shared<WilliamsonCase5>();
  } else if (test_case == "williamson6") {
    test_case_ptr = std::make_shared<WilliamsonCase6>();
  } else if (test_case == "galewsky") {
    test_case_ptr = std::make_shared<GalewskyCase>();
  } else if (test_case == "galewsky_Set_Initial") {
    test_case_ptr = std::make_shared<GalewskySetInitial>();
  } else {
    LOG << "unknown test case: " << test_case;
  }
  test_case_ptr->use_optimal_transport = true;
  test_case_ptr->add_artificial_viscosity =
      program.get<bool>("--add_artificial_viscosity");
  if (program.get<bool>("--use_explicit_time_stepping"))
    test_case_ptr->time_stepping = TimeSteppingScheme::kExplicit;
  test_case_ptr->project_velocity = true;
  test_case_ptr->stabilize_pressure_gradient = false;

  // initialize particles
  double* sites = nullptr;
  int n_smooth = 0;
  size_t n_sites = 0;
  std::vector<index_t> order;
  std::vector<double> heights;
  std::shared_ptr<Mesh> mesh;
  EarthProperties earth;

  if (import_height_from.empty()) {
    std::string particles = program.get<std::string>("--particles");
    if (size_t pos = particles.find("icosahedron") != std::string::npos) {
      particles.erase(pos - 1, 11);
      mesh = std::make_shared<SubdividedSphere<Icosahedron>>(
          std::atoi(particles.data()));
      n_sites = mesh->vertices().n();
      sites = mesh->vertices()[0];
    }

    // construct a better ordering of the points
    order.resize(n_sites);
    sort_points_on_zcurve(sites, n_sites, dim, order);
  } else {
    mesh = std::make_shared<Mesh>(3);

    // Galewsky jet implementation using Darren Engwirda's initial condition
    // https://github.com/dengwirda/swe-python
    std::ifstream f(import_height_from);
    nlohmann::json json;
    f >> json;

    std::vector<double> xs = json["x"];
    std::vector<double> ys = json["y"];
    std::vector<double> zs = json["z"];
    std::vector<double> hs = json["h"];

    size_t n = xs.size();
    ASSERT(ys.size() == n);
    ASSERT(zs.size() == n);
    ASSERT(hs.size() == n) << fmt::format("|hs| = {}, n = {}", hs.size(), n);

    const double a = earth.radius;
    heights.resize(n);
    std::array<coord_t, 3> coords;
    for (size_t i = 0; i < n; ++i) {
      coords[0] = xs[i] / a;
      coords[1] = ys[i] / a;
      coords[2] = zs[i] / a;
      heights[i] = hs[i];
      ASSERT(heights[i] > 0);
      mesh->vertices().add(coords.data());
    }
    n_sites = mesh->vertices().n();
    sites = mesh->vertices()[0];
    order.resize(n_sites);
    std::iota(order.begin(), order.end(), 0);
  }
  LOG << fmt::format("# sites = {}", n_sites);

  Vertices vertices(dim);
  vertices.reserve(n_sites);
  coord_t x[dim];
  for (size_t i = 0; i < n_sites; i++) {
    for (int d = 0; d < dim; d++) x[d] = sites[dim * order[i] + d];
    vertices.add(x);
  }

  // set up the solver
  ShallowWaterSimulation<Domain_t> solver(domain, n_sites, vertices[0],
                                          vertices.dim(), *test_case_ptr);

  // set up simulation
  SimulationOptions solver_opts;
  solver_opts.save_initial_mesh = true;
  solver_opts.n_smoothing_iterations = n_smooth;
  solver_opts.advect_from_centroid = true;
  solver.initialize(domain, solver_opts);
  solver.setup();

  if (!import_height_from.empty()) {
    for (size_t k = 0; k < n_sites; k++) {
      ASSERT(heights[k] > 0);
      solver.height()[k] = heights[k];
    }
    solver.calculate_initial_conservation_quantities();
  }

  solver.statistics().n_particles = n_sites;
  solver.statistics().name = test_case;

  // step in time
  double days = test_case_ptr->days;
  if (program.present<double>("--days")) days = program.get<double>("--days");
  solver_opts.time_step = program.get<double>("--step");
  solver_opts.verbose = false;
  solver_opts.backtrack = false;
  solver_opts.restart_zero_weights = true;
  solver_opts.max_iter = 3;
  solver_opts.skip_initial_calculation = true;
  solver_opts.neighbor_algorithm = NearestNeighborAlgorithm::kSphereQuadtree;

  Timer timer;
  timer.start();
  int save_every = program.get<int>("--save_every");
  double seconds = 0;
  double hour = 0;
  solver.save(prefix + "0.vtk");
  solver.save_json(fmt::format("{}{}.json", prefix, hour));
  while (seconds < days_to_seconds(days)) {
    double dt = solver.time_step(solver_opts);
    solver_opts.iteration++;
    seconds += dt;
    solver_opts.time = seconds;
    int current_hour = seconds / 3600;
    if (current_hour == hour + 1) {
      if (current_hour % save_every == 0) {
        solver.save(fmt::format("{}{}.vtk", prefix, current_hour));
        solver.save_json(fmt::format("{}{}.json", prefix, current_hour));
      }
      ++hour;
    }
  }
  timer.stop();
  solver.statistics().total_time = timer.seconds();

  std::ofstream outfile(output_dir + "/" +
                        program.get<std::string>("--statistics"));
  outfile << std::setw(4) << solver.statistics().to_json() << std::endl;
}

nlohmann::json ShallowWaterStatistics::to_json() const {
  nlohmann::json data;
  data["n_particles"] = n_particles;
  data["name"] = name;
  data["rw"] = rw;
  data["ra"] = ra;
  data["rm"] = rm;
  data["rp"] = rp;
  data["re"] = re;
  data["time"] = time;
  data["h_error"] = h_error;
  data["h_total"] = h_total;
  data["u_error"] = u_error;
  data["u_total"] = u_total;
  data["sdpd"] = sdpd;
  data["voronoi_time"] = voronoi_time;
  data["n_voronoi"] = n_voronoi;
  data["linear_solver_time"] = linear_solver_time;
  data["time_step_time"] = time_step_time;
  data["total_time"] = total_time;
  return data;
}

template class ShallowWaterSimulation<SphereDomain>;

}  // namespace vortex