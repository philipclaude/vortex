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
#pragma once

#include <vector>

#include "math/spmat.h"
#include "voronoi.h"

namespace vortex {

struct SimulationOptions {
  double volume_grad_tol{1e-8};
  int max_iter{25};
  int n_neighbors{80};
  bool verbose{true};
  double time_step{1e-3};
  double time{0};
  int iteration{0};
  bool backtrack{false};
  double epsilon{1};
  bool save_initial_mesh{false};
  bool backtrack_check_gradient{false};
  bool backtrack_check_zero_volume{false};
  bool restart_zero_weights{false};
  int print_frequency{50};
  bool reflection_boundary_condition{false};
  bool advect_from_centroid{true};
};

struct SimulationConvergence {
  double error{1};
  bool converged{false};
  int n_iterations{0};
  double min_volume{0};
};

struct FluidProperties {
  double viscosity{0};
};

struct Particle {
  double mass;
  vec3d velocity;
  vec3d position;
  vec3d centroid;
};

template <typename Domain_t>
bool has_boundary();

template <typename Domain_t>
void project_point(double* x);

template <typename Domain_t>
void project_velocity(double* v);

class Particles : public Vertices {
 public:
  Particles(size_t np, const coord_t* xp, int dim)
      : Vertices(4),
        volume_(np),
        mass_(np),
        density_(np),
        pressure_(np),
        centroids_(3),
        velocity_(3) {
    create(np, xp, dim);
  }

  template <typename VelocityFunction>
  void set_velocity(const VelocityFunction& fn) {
    velocity_.clear();
    velocity_.reserve(n());
    for (size_t k = 0; k < n(); k++) {
      auto vk = fn((*this)[k]);
      velocity_.add(&vk[0]);
    }
  }

  template <typename DensityFunction>
  void set_density(const DensityFunction& fn) {
    for (size_t k = 0; k < n(); k++) density_[k] = fn((*this)[k]);
  }

  void create(size_t np, const coord_t* xp, int dim);
  void save(const std::string& name) const;

  const auto& volume() const { return volume_; }
  auto& volume() { return volume_; }
  const auto& mass() const { return mass_; }
  auto& mass() { return mass_; }
  const auto& pressure() const { return pressure_; }
  auto& pressure() { return pressure_; }
  const auto& density() const { return density_; }
  auto& density() { return density_; }

  const auto& centroids() const { return centroids_; }
  auto& centroids() { return centroids_; }

  const auto& velocity() const { return velocity_; }
  auto& velocity() { return velocity_; }

  Particle particle(size_t k) const {
    Particle p;
    for (int d = 0; d < 3; d++) {
      p.velocity[d] = velocity_(k, d);
      p.position[d] = (*this)(k, d);
      p.centroid[d] = centroids_(k, d);
    }
    p.mass = mass_[k];
    return p;
  }

 private:
  vecd<double> volume_;
  vecd<double> mass_;
  vecd<double> density_;
  vecd<double> pressure_;
  Vertices centroids_;
  array2d<double> velocity_;
};

class ParticleSimulation {
 public:
  ParticleSimulation(size_t np, const coord_t* xp, int dim)
      : particles_(np, xp, dim),
        voronoi_(particles_.dim(), particles_[0], np),
        max_displacement_(np),
        hessian_(np, np),
        dw_(np),
        gradient_(np) {}

  template <typename Domain_t>
  void initialize(const Domain_t& domain, SimulationOptions sim_opts) {
    // we need to calculate the voronoi diagram to obtain initial volumes,
    // masses and centroids
    VoronoiDiagramOptions voro_opts;
    voro_opts.store_mesh = sim_opts.save_initial_mesh;
    voro_opts.store_facet_data = true;
    voro_opts.allow_reattempt = false;
    voro_opts.n_neighbors = sim_opts.n_neighbors;
    voro_opts.verbose = false;
    voronoi_.weights().resize(particles_.n(), 0);
    voronoi_.compute(domain, voro_opts);

    // save the initial volume and mass of each particles
    for (size_t k = 0; k < particles_.n(); k++) {
      particles_.volume()[k] = voronoi_.properties()[k].volume;
      particles_.mass()[k] = particles_.density()[k] * particles_.volume()[k];
      ASSERT(particles_.mass()[k] > 0);
    }
  }

  template <typename Domain_t>
  void calculate_power_diagram(const Domain_t& domain,
                               VoronoiDiagramOptions opts) {
    lift_sites(particles_, voronoi_.weights());
    voronoi_.compute(domain, opts);
  }

  void compute_gradient(const std::vector<double>& target_volume) {
    for (size_t k = 0; k < particles_.n(); k++)
      gradient_[k] = target_volume[k] - voronoi_.properties()[k].volume;
  }

  template <typename Domain_t>
  SimulationConvergence optimize_volumes(
      const Domain_t& domain, SimulationOptions sim_opts,
      const std::vector<double>& target_volume) {
    // calculate the voronoi diagram
    VoronoiDiagramOptions voro_opts;
    voro_opts.store_mesh = false;
    voro_opts.store_facet_data = true;
    voro_opts.allow_reattempt = false;
    voro_opts.n_neighbors = sim_opts.n_neighbors;
    voro_opts.verbose = false;

    // utility to get the minimum volume in the voronoi diagram
    auto min_volume = [this]() -> double {
      auto props = *std::min_element(
          voronoi_.properties().begin(), voronoi_.properties().end(),
          [](const auto& pa, const auto& pb) { return pa.volume < pb.volume; });
      return props.volume;
    };

  start_optimization:
    if (sim_opts.restart_zero_weights)
      for (size_t k = 0; k < voronoi_.weights().size(); k++)
        voronoi_.weights()[k] = 0.0;
    const auto nu_min =
        *std::min_element(target_volume.begin(), target_volume.end());

    // iterate until converged to the desired mass
    SimulationConvergence convergence;
    convergence.error = 1;
    convergence.converged = false;
    double v_min = std::numeric_limits<double>::max();
    for (convergence.n_iterations = 0;
         convergence.n_iterations < sim_opts.max_iter;
         convergence.n_iterations++) {
      // calculate the power diagram with the current weights
      calculate_power_diagram(domain, voro_opts);

      // calculate the minimum volume of the voronoi diagram (zero weights)
      if (convergence.n_iterations == 0) {
        v_min = min_volume();
        if (sim_opts.restart_zero_weights)
          ASSERT(v_min > 0)
              << fmt::format("Voronoi diagram has an empty cell.");
      }

      // calculate gradient and check for convergence
      compute_gradient(target_volume);
      convergence.error = length(gradient_);
      if (sim_opts.verbose)
        LOG << fmt::format("iter[{:3}]: error = {:.3e}",
                           convergence.n_iterations, convergence.error);
      if (convergence.error < sim_opts.volume_grad_tol) {
        convergence.converged = true;
        break;
      }

      // solve hessian_ * dw_ = gradient_ for dw_ (search direction)
      compute_search_direction(target_volume);

      // backtrack to make sure cells always have an acceptable volume
      const double initial_gradient_norm = convergence.error;
      const double a0 = 0.5 * std::min({v_min, nu_min});
      std::vector<double> initial_weights(voronoi_.weights().begin(),
                                          voronoi_.weights().end());
      double alpha = 1.0;
      while (true) {
        // update weights using the search direction
        for (size_t k = 0; k < voronoi_.weights().size(); k++)
          voronoi_.weights()[k] = initial_weights[k] - alpha * dw_[k];
        if (!sim_opts.backtrack) break;  // option to ignore volume/grad check

        // get the minimum volume and gradient norm with these weights
        calculate_power_diagram(domain, voro_opts);
        compute_gradient(target_volume);
        const auto gradient_norm = length(gradient_);
        const auto v_min_alpha = min_volume();

        // the third conditional is the full check in Algorithm 1 of
        // Kitagawa (2019) https://arxiv.org/pdf/1603.05579
        if (!sim_opts.backtrack_check_gradient &&
            sim_opts.backtrack_check_zero_volume && v_min_alpha > 0)
          break;
        if (!sim_opts.backtrack_check_gradient && v_min_alpha > a0) break;
        if (v_min_alpha > a0 &&
            gradient_norm <= (1 - 0.5 * alpha) * initial_gradient_norm)
          break;
        alpha /= 2.0;
        if (alpha < 1e-8) {
          LOG << fmt::format("[ERROR] backtracking: min vol = {}", v_min_alpha);
          break;
        }
      }
    }
    convergence.min_volume = min_volume();
    if (convergence.min_volume <= 0) {
      // restart the optimization with backtracking enabled (if not already)
      if (!sim_opts.backtrack) {
        sim_opts.backtrack = true;
        sim_opts.verbose = true;
        for (size_t k = 0; k < voronoi_.weights().size(); k++)
          voronoi_.weights()[k] = 0.0;
        goto start_optimization;
      } else {
        NOT_POSSIBLE;  // this is a fatal error in the simulation
      }
    }
    return convergence;
  }

  void compute_search_direction(const std::vector<double>& target_volumes);

  // store mass, volume, centroids, etc. based on the current Voronoi diagram
  void calculate_properties();

  template <typename ForceFunction>
  double max_time_step(const ForceFunction& fext, SimulationOptions opts) {
    auto len = [](vec3d x) -> double {
      return std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    };
    double eps = opts.epsilon;  // inverse of spring constant
    double dt = opts.time_step;
    for (size_t k = 0; k < particles_.n(); k++) {
      vec3d f = fext(particles_.particle(k));  // initialize to external force
      vec3d v(particles_.velocity()[k]);
      for (int d = 0; d < 3; d++)
        f[d] += (particles_.centroids()(k, d) - particles_(k, d)) / (eps * eps);

      double f0 = len(f);
      double v0 = len(v);
      double dx = max_displacement_[k];
      double m = particles_.mass()[k];
      if (v0 < 1e-12) continue;  // any step is allowed
      if (f0 < 1e-12)
        dt = dx / v0;
      else
        dt = 0.5 * m * (std::sqrt(v0 * v0 + 4 * dx * f0 / m) - v0) / f0;
    }
    return 0.95 * dt;  // 95% of permissible time step, TODO: make user option
  }

  auto& particles() { return particles_; }
  const auto& voronoi() const { return voronoi_; }
  auto& voronoi() { return voronoi_; }

 protected:
  Particles particles_;
  VoronoiDiagram voronoi_;
  vecd<double> max_displacement_;

 private:
  spmat<double> hessian_;
  vecd<double> dw_;
  vecd<double> gradient_;
};

template <typename Domain_t>
class SpringParticles : public ParticleSimulation {
 public:
  SpringParticles(const Domain_t& domain, int np, const coord_t* xp, int dim)
      : ParticleSimulation(np, xp, dim), domain_(domain) {
    timer_.start();
  }

  // perform one time step in the simulation
  template <typename ForceFunction>
  void step(const ForceFunction& fext, FluidProperties properties,
            SimulationOptions& options) {
    calculate_properties();  // mass, volume, centroids

    // get the maximum time step
    double dt = options.time_step;
    double dt_max = max_time_step(fext, options);
    if (dt_max < options.time_step) dt = dt_max;

    auto& velocity = particles_.velocity();
    double eps = options.epsilon;  // inverse of spring constant
    double dx_max = 0, v_max = 0;
    double dx_tot = 0;
    for (size_t k = 0; k < particles_.n(); k++) {
      vec3d f = fext(particles_.particle(k));  // initialize to external force

      // update velocity
      double v_k = 0;
      for (int d = 0; d < 3; d++) {
        // add spring force
        f[d] += (particles_.centroids()(k, d) - particles_(k, d)) / (eps * eps);

        // update velocity
        velocity(k, d) += dt * f[d] / particles_.mass()[k];
        v_k += velocity(k, d) * velocity(k, d);
      }
      // project_velocity<Domain_t>(velocity[k]);
      v_k = std::sqrt(v_k);
      if (v_k > v_max) v_max = v_k;

      // update position
      double dx_k = 0;
      for (int d = 0; d < 3; d++) {
        double dx_d = dt * velocity(k, d);
        dx_k += dx_d * dx_d;
        if (options.advect_from_centroid)
          particles_(k, d) = particles_.centroids()(k, d) + dx_d;
        else
          particles_(k, d) = particles_(k, d) + dx_d;
      }
      project_point<Domain_t>(particles_[k]);
      dx_k = std::sqrt(dx_k);
      if (dx_k > dx_max) dx_max = dx_k;
      dx_tot += dx_k;

      // apply boundary conditions
      // TODO(philip) more general boundary condition treatment
      // for now just reflect off the boundary (assuming a square)
      if (options.reflection_boundary_condition && has_boundary<Domain_t>()) {
        for (int d = 0; d < 2; d++) {
          if (particles_(k, d) < 0.0) {
            particles_(k, d) = 0.0;
            velocity(k, d) *= -1.0;
          }
          if (particles_(k, d) > 1.0) {
            particles_(k, d) = 1.0;
            velocity(k, d) *= -1.0;
          }
        }
      }
    }
    dx_tot = std::sqrt(dx_tot);

    std::vector<double> target_volume(particles_.n());
    for (size_t k = 0; k < particles_.n(); k++) {
      target_volume[k] = particles_.mass()[k] / particles_.density()[k];
      ASSERT(target_volume[k] > 0);
    }

    options.max_iter = 30;
    auto convergence = optimize_volumes(domain_, options, target_volume);
    if (options.iteration % options.print_frequency == 0) {
      timer_.stop();
      simulation_rate_ = double(options.print_frequency / timer_.seconds());
      print_footer();
      print_header();
      timer_.start();
    }
    std::cout << fmt::format(
        "| {:6d} | {:1.3e} | {:1.3e} | {:1.3e} | {:1.3e} | {:3d} | {:1.3e} | "
        "{:1.3e} | {:1.3e} |\n",
        options.iteration, options.time, dt_max, dt, convergence.error,
        convergence.n_iterations, dx_max, dx_tot, v_max);
    options.time += dt;
  }

  void print_header(int n_bars = 100) const {
    std::cout << fmt::format(
        "| {:6s} | {:9s} | {:9s} | {:9s} | {:9s} | {:3s} | {:9s} | {:9s} | "
        "{:9s} | @{:3.1f} steps / sec.\n",
        "step", "time", "dt_max", "dt", "|rm|", "nm", "max(dx)", "|dx|",
        "max(v)", simulation_rate_);
    std::cout << fmt::format("{:->{}}", "", n_bars) << std::endl;
  }

  void print_footer(int n_bars = 100) const {
    std::cout << fmt::format("{:->{}}", "", n_bars) << std::endl;
  }

 private:
  const Domain_t& domain_;
  Timer timer_;
  double simulation_rate_{0};  // # time steps per second (of wall-clock time)
};

}  // namespace vortex