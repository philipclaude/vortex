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
  vecd<double> mass_;     // ideally constant during the simulation
  vecd<double> density_;  // constant for incompressible fluid
  vecd<double> pressure_;
  Vertices centroids_;
  array2d<double> velocity_;
};

class ParticleSimulation {
 public:
  ParticleSimulation(size_t np, const coord_t* xp, int dim)
      : particles_(np, xp, dim),
        voronoi_(particles_.dim(), particles_[0], np),
        hessian_(np, np),
        dw_(np),
        gradient_(np) {}

  template <typename Domain_t>
  void initialize(const Domain_t& domain, SimulationOptions sim_opts) {
    // we need to calculate the voronoi diagram to obtain initial volumes,
    // masses and centroids
    VoronoiDiagramOptions voro_opts;
    voro_opts.store_mesh = false;
    voro_opts.store_facet_data = true;
    voro_opts.allow_reattempt = false;
    voro_opts.n_neighbors = sim_opts.n_neighbors;
    voro_opts.verbose = false;
    voronoi_.weights().resize(particles_.n(), 0);
    voronoi_.compute(domain, voro_opts);
    for (size_t k = 0; k < particles_.n(); k++) {
      double m = particles_.density()[k] * voronoi_.properties()[k].volume;
      particles_.mass()[k] = m;
      ASSERT(particles_.mass()[k] > 0);
    }
  }

  void compute_gradient(const std::vector<double>& target_volumes) {
    for (size_t k = 0; k < particles_.n(); k++)
      gradient_[k] = target_volumes[k] - voronoi_.properties()[k].volume;
  }

  template <typename Domain_t>
  SimulationConvergence optimize_volumes(
      const Domain_t& domain, SimulationOptions sim_opts,
      const std::vector<double>& target_volumes) {
    // calculate the voronoi diagram
    VoronoiDiagramOptions voro_opts;
    voro_opts.store_mesh = false;
    voro_opts.store_facet_data = true;
    voro_opts.allow_reattempt = false;
    voro_opts.n_neighbors = sim_opts.n_neighbors;
    voro_opts.verbose = false;
    // voro_opts.parallel = false;

    // utility to get the minimum volume in the voronoi diagram
    auto min_volume = [this]() {
      auto props = *std::min_element(
          voronoi_.properties().begin(), voronoi_.properties().end(),
          [](const auto& pa, const auto& pb) { return pa.volume < pb.volume; });
      return props.volume;
    };

    // set initial guess
    voronoi_.weights().resize(particles_.n(), 0.0);
    const double nu_min =
        *std::min_element(target_volumes.begin(), target_volumes.end());

    // iterate until converged to the desired mass
    SimulationConvergence convergence;
    convergence.error = 1;
    convergence.converged = false;
    double v_min = std::numeric_limits<double>::max();
    for (convergence.n_iterations = 0;
         convergence.n_iterations < sim_opts.max_iter;
         convergence.n_iterations++) {
      // calculate the power diagram with the current weights
      lift_sites(particles_, voronoi_.weights());
      voronoi_.compute(domain, voro_opts);

      // calculate the minimum volume of the voronoi diagram (zero weights)
      if (convergence.n_iterations == 0) v_min = min_volume();

      // calculate gradient and check for convergence
      compute_gradient(target_volumes);
      convergence.error = length(gradient_);
      if (sim_opts.verbose)
        LOG << fmt::format("iter[{:3}]: error = {:.3e}",
                           convergence.n_iterations, convergence.error);
      if (convergence.error < sim_opts.volume_grad_tol) {
        convergence.converged = true;
        break;
      }

      // compute the search direction
      compute_search_direction(target_volumes);

      // backtrack to make sure cells always have a positive volume
      const double initial_gradient_norm = convergence.error;
      const double a0 = 0.5 * std::min({v_min, nu_min});
      std::vector<double> initial_weights(voronoi_.weights().begin(),
                                          voronoi_.weights().end());
      double alpha = 1.0;
      while (true) {
        for (size_t k = 0; k < voronoi_.weights().size(); k++)
          voronoi_.weights()[k] = initial_weights[k] - alpha * dw_[k];
        if (!sim_opts.backtrack) break;
        voronoi_.compute(domain, voro_opts);
        compute_gradient(target_volumes);
        const double gradient_norm = length(gradient_);

        // get the minimum volume with this step size
        const double v_min_alpha = min_volume();

        if (v_min_alpha > a0 &&
            gradient_norm <= (1 - 0.5 * alpha) * initial_gradient_norm)
          break;
        alpha /= 2.0;
      }
    }
    convergence.min_volume = min_volume();
    return convergence;
  }

  void compute_search_direction(const std::vector<double>& target_volumes);

  // store mass, volume, centroids, etc. based on the current Voronoi diagram
  void calculate_properties();

  auto& particles() { return particles_; }
  const auto& voronoi() const { return voronoi_; }
  auto& voronoi() { return voronoi_; }

 protected:
  Particles particles_;
  VoronoiDiagram voronoi_;

 private:
  spmat<double> hessian_;
  vecd<double> dw_;
  vecd<double> gradient_;
};

template <typename Domain_t>
class SpringParticles : public ParticleSimulation {
 public:
  SpringParticles(const Domain_t& domain, int np, const coord_t* xp, int dim)
      : ParticleSimulation(np, xp, dim), domain_(domain) {}

  // perform one time step in the simulation
  template <typename ForceFunction>
  void step(const ForceFunction& fext, FluidProperties properties,
            SimulationOptions options) {
    calculate_properties();  // mass, volume, centroids

    auto& velocity = particles_.velocity();
    double eps = 4e-3;  // spring constant
    double dt = options.time_step;
    for (size_t k = 0; k < particles_.n(); k++) {
      vec3d f = fext(particles_.particle(k));  // initialize to external force
      for (int d = 0; d < 3; d++) {
        // add spring force
        f[d] += (particles_.centroids()(k, d) - particles_(k, d)) / (eps * eps);

        // update velocity and position
        velocity(k, d) += dt * f[d] / particles_.mass()[k];
        particles_(k, d) += dt * velocity(k, d);
      }

      // apply boundary conditions
      // TODO(philip) more general boundary condition treatment
      // for now just reflect off the boundary (assuming a square)
      if (has_boundary<Domain_t>()) {
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

      // project to the domain
      project_point<Domain_t>(particles_[k]);
      // project_velocity<Domain_t>(velocity[k]);
    }

    // particles_.print();

    std::vector<double> target_volume(particles_.n());
    for (size_t k = 0; k < particles_.n(); k++) {
      target_volume[k] = particles_.mass()[k] / particles_.density()[k];
      ASSERT(target_volume[k] > 0);
    }

    options.max_iter = 30;
    auto convergence = optimize_volumes(domain_, options, target_volume);
    std::string status = convergence.converged ? "converged" : "unconverged";
    if (options.iteration % 10 == 0)
      LOG << "_______________________________________________";
    LOG << fmt::format("{:3d} | {:1.4e} | {:2d} | {:.3e} | {}",
                       options.iteration, options.time,
                       convergence.n_iterations, convergence.error, status);
  }

 private:
  const Domain_t& domain_;
};

template <typename Domain_t>
class PowerParticles : public ParticleSimulation {
 public:
  PowerParticles(const Domain_t& domain, int np, const coord_t* xp, int dim)
      : ParticleSimulation(np, xp, dim),
        domain_(domain),
        vstar_(3),
        grad_{spmat<double>(np, np), spmat<double>(np, np),
              spmat<double>(np, np)},
        div_{spmat<double>(np, np), spmat<double>(np, np),
             spmat<double>(np, np)},
        laplacian_(np, np),
        aux_matrix_(np, np),
        aux_vector_(np),
        aux_rhs_(np) {
    reserve();
  }

  void reserve() {
    vstar_.reserve(particles_.n());
    vec3d uvw;
    for (size_t k = 0; k < particles_.n(); k++) vstar_.add(&uvw[0]);
  }

  // perform one time step in the simulation
  template <typename ForceFunction>
  void step(const ForceFunction& fext, FluidProperties properties,
            SimulationOptions options) {
    calculate_properties();  // mass, volume, centroids
    build_operators();       // divergence, gradient, Laplacian
    apply_external_forces(options, properties, fext);  // calculate v*
    compute_pressure(options.time_step);               // calculate p
    apply_internal_forces(options.time_step);  // apply vn = v* - grad(p)
    advect(options.time_step);                 // apply x = x + dt.vn

    std::vector<double> target_volume(particles_.n());
    for (size_t k = 0; k < particles_.n(); k++)
      target_volume[k] = particles_.mass()[k] / particles_.density()[k];

    options.max_iter = 30;
    auto convergence = optimize_volumes(domain_, options, target_volume);
    std::string status = convergence.converged ? "converged" : "unconverged";
    if (options.iteration % 10 == 0)
      LOG << "_______________________________________________";
    LOG << fmt::format("{:3d} | {:1.4e} | {:2d} | {:.3e} | {}",
                       options.iteration, options.time,
                       convergence.n_iterations, convergence.error, status);
  }

  // compute fractional velocity (v*) according to external forces
  template <typename ForceFunction>
  void apply_external_forces(SimulationOptions options,
                             FluidProperties properties,
                             const ForceFunction& fext) {
    // evaluate the external forces
    double dt = options.time_step;
    if (properties.viscosity > 0) {
      // TODO(philip): use OpenNL multiple RHS feature
      for (int d = 0; d < 3; d++) {
        // build rhs vector for this dimension
        for (size_t k = 0; k < particles_.n(); k++) {
          auto fk = fext(particles_.particle(k));
          aux_rhs_[k] =
              particles_.volume()[k] *
              (particles_.velocity()[k][d] + dt * fk[d] / particles_.mass()[k]);
        }

        // build the matrix for this dimension
        aux_matrix_.clear();
        const auto& rows = laplacian_.rows();
        ASSERT(rows.size() == particles_.n());
        for (size_t row = 0; row < particles_.n(); row++) {
          for (const auto& [col, value] : rows[row]) {
            aux_matrix_(row, col) = -properties.viscosity * dt * value;
          }
          aux_matrix_(row, row) += particles_.volume()[row];
        }

        // solve for the velocity update
        aux_matrix_.solve_nl(aux_rhs_, aux_vector_, 1e-10, false);
        for (size_t k = 0; k < particles_.n(); k++)
          for (int d = 0; d < 3; d++) vstar_[k][d] = aux_vector_[k];
      }
    } else {
      for (size_t k = 0; k < particles_.n(); k++) {
        auto fk = fext(particles_.particle(k));
        for (int d = 0; d < 3; d++)
          vstar_[k][d] =
              particles_.velocity()[k][d] + dt * fk[d] / particles_.mass()[k];
      }
    }
  }

  // solve the pressure projection equation
  void compute_pressure(double dt);

  // add the contribution of internal pressure forces
  void apply_internal_forces(double dt);

  // advect the particles
  void advect(double dt);

  // populates the gradient, divergence and Laplacian matrices
  void build_operators();

 private:
  const Domain_t& domain_;
  array2d<double> vstar_;
  spmat<double> grad_[3];
  spmat<double> div_[3];
  spmat<double> laplacian_;

  // for solving arbitrary systems of equations
  spmat<double> aux_matrix_;
  vecd<double> aux_vector_;
  vecd<double> aux_rhs_;
};

}  // namespace vortex