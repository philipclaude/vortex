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

#include "io.h"
#include "math/spmat.h"
#include "operators.h"

namespace vortex {

template <typename Domain_t>
void ShallowWaterSimulation<Domain_t>::setup() {
  for (size_t k = 0; k < particles_.n(); k++) {
    height_[k] = options_.initial_height(particles_[k]) -
                 options_.surface_height(particles_[k]);
  }
  particles_.set_velocity(options_.initial_velocity);
  particles_.set_density([](const double* x) { return 1.0; });
  for (size_t k = 0; k < particles_.n(); k++) {
    volume_[k] = height_[k] * voronoi_.properties()[k].volume;
  }
  initial_mass_ = total_mass();
  initial_momentum_ = total_momentum();
  initial_energy_ = total_energy();
}

template <typename Domain_t>
double ShallowWaterSimulation<Domain_t>::forward_euler_step(
    const SimulationOptions& options) {
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
    double cm = std::sqrt(g * height_[k]);
    // TODO this should be done from the advection point (site or centroid)
    // currently max displacement is measured from the site
    double dt_k = a * max_displacement_[k] / (um + cm);
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
  for (size_t k = 0; k < n; k++) {
    height[k] = height_[k] + options_.surface_height(particles_[k]);
  }
  ops.calculate_gradient(height.data(), grad_h.data());

  // advance the height field according to governing differential equation
  for (size_t k = 0; k < n; k++) {
    ASSERT(height_[k] > 0);
    height_[k] -= dt * height_[k] * div_u[k] / a;
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

  // smooth the particles
  for (int iter = 0; iter < 0; iter++) {
    voronoi_.smooth(particles_, true);
    VoronoiDiagramOptions voro_opts;
    voro_opts.verbose = false;
    calculate_power_diagram(domain_, voro_opts);
  }

  SimulationConvergence convergence;
  if (options_.conserve_mass) {
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
  }

  // compute artificial viscosity
  std::vector<double> viscous(3 * n, 0.0);
  if (options_.add_artificial_viscosity) compute_artificial_viscosity(viscous);
  // stabilize_pressure(height, grad_h);

  // update velocity
  if (!options_.use_analytic_velocity) {
    for (size_t i = 0; i < n; i++) {
      vec3d x(particles_[i]);  // on unit sphere
      vec3d c(particles_.centroids()[i]);
      vec3d u(velocity[i]);
      double f = options_.coriolis_parameter(&x[0]);
      vec3d force = f * cross(x, u);  // coriolis force
      if (options_.constrain) force = force + dot(u, u) * x / a;
      vec3d fv(&viscous[3 * i]);  // derivative needs to be scaled by 1/a
      force = force + 1.0 * fv / a;

      vec3d fs = options_.spring_stiffness * (x - c);
      force = force + fs;

      for (int d = 0; d < 3; d++) {
        velocity(i, d) -= dt * (g * grad_h[3 * i + d] / a + force[d]);
      }
      if (options_.project_velocity)
        project_velocity<Domain_t>(particles_[i], velocity[i]);
    }
  } else {
    // use the analytic velocity, e.g. for Williamson Case 1
    for (size_t i = 0; i < n; i++) {
      vec3d u = options_.analytic_velocity(particles_[i]);
      for (int d = 0; d < 3; d++) {
        velocity(i, d) = u[d];
      }
    }
  }

  // monitor the simulation
  if (options.iteration % options.print_frequency == 0) print_header();
  double sdpd = options.time / timer_.seconds();
  double area_error = (domain_.area() - total_area());
  double mass_error = 100 * (total_mass() - initial_mass_) / initial_mass_;
  double momentum_error =
      100 * (total_momentum() - initial_momentum_) / initial_momentum_;
  double energy_error =
      100 * (total_energy() - initial_energy_) / initial_energy_;
  std::cout << fmt::format(
      "| {:6d} | {:9s} | {:1.1e} | {:2d}:{:1.1e} | {:+1.1e} | {:+1.1e} | "
      "{:+1.1e} "
      "| {:+1.1e} | {:6.1f} | \n",
      options.iteration, days_hours_minutes(options.time), dt,
      convergence.n_iterations, convergence.error, area_error, mass_error,
      momentum_error, energy_error, sdpd);
  return dt;
}

template <typename Domain_t>
void ShallowWaterSimulation<Domain_t>::compute_artificial_viscosity(
    std::vector<double>& fv) {
  const size_t n = particles_.n();
  double g = earth_.gravity;
  double a = earth_.radius;
  double eps = 1e-6;

  fv.resize(3 * n, 0);
  ASSERT(voronoi_.facets().size() > 0);
  for (const auto& facet : voronoi_.facets()) {
    if (facet.bj < 0) continue;
    if (facet.bi >= particles_.n()) continue;
    if (facet.bj >= particles_.n()) continue;
    size_t i = facet.bi;
    size_t j = facet.bj;
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
void ShallowWaterSimulation<Domain_t>::stabilize_pressure(
    const std::vector<double>& h, std::vector<double>& dh) {
  const size_t n = particles_.n();

  // compute laplacian of height field
  std::vector<double> lh(n, 0);
  ASSERT(voronoi_.facets().size() > 0);
  for (const auto& facet : voronoi_.facets()) {
    if (facet.bj < 0) continue;
    if (facet.bi >= particles_.n()) continue;
    if (facet.bj >= particles_.n()) continue;
    size_t i = facet.bi;
    size_t j = facet.bj;
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
      "{:6s} |\n",
      "Step", "day:hr:mn", "dt (s)", "Rw", "Ra (%)", "Rm (%)", "Rp (%) ",
      "Re (%)", "SDPD");
  std::cout << fmt::format("{:->{}}", "", n_bars) << std::endl;
}

template <typename Domain_t>
void ShallowWaterSimulation<Domain_t>::save(const std::string& filename) const {
  size_t n = particles_.n();
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
      //(height_[k] + options_.surface_height(particles_[k])) / 1000;
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
  fprintf(fid, "FIELD FieldData 1\n");
  fprintf(fid, "height 1 %zu float\n", n);
  std::vector<float> height_data(n);
  for (size_t k = 0; k < n; k++) {
    height_data[k] = height_[k] + options_.surface_height(particles_[k]);
    // vec3d u(particles_.velocity()[k]);
    //  height_data[k] = length(u);
  }
  std::parafor_i(0, n, [&height_data](int tid, size_t k) {
    io::swap_end(height_data[k]);
  });
  fwrite(&height_data[0], 4, n, fid);
  fprintf(fid, "\n");
  fclose(fid);
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
  for (int k = 0; k < particles_.n(); k++) {
    mass += voronoi_.properties()[k].volume * height_[k];
  }
  return mass;
}

template <typename Domain_t>
double ShallowWaterSimulation<Domain_t>::total_momentum() const {
  double g = earth_.gravity;
  double momentum = 0;
  for (int k = 0; k < particles_.n(); k++) {
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
  for (int k = 0; k < particles_.n(); k++) {
    double vi = voronoi_.properties()[k].volume;
    vec3d u(particles_.velocity()[k]);
    double hi = height_[k];
    double hs = options_.surface_height(particles_[k]);
    energy += 0.5 * vi * hi * dot(u, u);
    energy += 0.5 * vi * g * ((hi + hs) * (hi + hs) - hs * hs);
  }
  return energy;
}

template class ShallowWaterSimulation<SphereDomain>;

}  // namespace vortex