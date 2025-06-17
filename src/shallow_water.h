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
#pragma once

#include <nlohmann/json_fwd.hpp>

#include "particles.h"

namespace argparse {
class ArgumentParser;
}

namespace vortex {

struct EarthProperties {
  const double radius = 6371220;
  const double gravity = 9.80616;
  const double angular_velocity = 7.292e-5;
};

static inline std::string days_hours_minutes(double s) {
  int d = s / (24 * 3600);
  s -= d * 24 * 3600;
  int h = s / 3600;
  s -= h * 3600;
  int m = s / 60;
  return fmt::format("{:03d}:{:02d}:{:02d}", d, h, m);
}

static inline int days_to_seconds(double d) { return d * 24 * 3600; }

typedef std::function<double(const coord_t*)> ScalarFunction;
typedef std::function<double(const coord_t*, double)> ScalarFunctionXT;
typedef std::function<vec3d(const coord_t*)> VectorFunction;

enum TimeSteppingScheme : uint8_t { kExplicit, kSemiImplicit };

struct ShallowWaterOptions {
  bool project_points{false};
  bool project_velocity{true};
  bool advect_from_centroid{true};
  bool use_analytic_velocity{false};
  bool use_analytic_height{false};
  bool use_optimal_transport{true};
  bool smoothing_iterations{0};
  bool constrain{true};
  bool add_artificial_viscosity{false};
  bool stabilize_pressure_gradient{false};
  TimeSteppingScheme time_stepping{TimeSteppingScheme::kSemiImplicit};
  ScalarFunction surface_height;
  ScalarFunction initial_height;
  ScalarFunctionXT analytic_height;
  VectorFunction surface_height_gradient;
  bool has_analytic_height{false};
  double days{1};

  VectorFunction initial_velocity;
  VectorFunction analytic_velocity;

  ScalarFunction coriolis_parameter;

  bool save_latlon{true};
  EarthProperties earth;
};

struct ShallowWaterStatistics {
  std::vector<double> rw;       // weight residual
  std::vector<double> ra;       // area residual
  std::vector<double> rm;       // mass residual
  std::vector<double> rp;       // momentum residual
  std::vector<double> re;       // energy residual
  std::vector<double> time;     // current time
  std::vector<double> h_error;  // error in analytic solution (if available)
  std::vector<double> sdpd;     // simulated day per day
  std::vector<double> voronoi_time;
  std::vector<double> linear_solver_time;
  std::vector<double> time_step_time;
  std::vector<int> n_voronoi;
  std::string name;
  size_t n_particles{0};
  double total_time{0};
  nlohmann::json to_json() const;
};

template <typename Domain_t>
class ShallowWaterSimulation : public ParticleSimulation {
 public:
  ShallowWaterSimulation(const Domain_t& domain, int np, const coord_t* xp,
                         int dim, const ShallowWaterOptions& test_case)
      : ParticleSimulation(np, xp, dim),
        domain_(domain),
        height_(np, 0.0),
        volume_(np, 0.0),
        options_(test_case),
        h_mat_(np, np),
        h_rhs_(np),
        h_sol_(np) {
    timer_.start();
    for (int i = 0; i < np; i++) particles_.density()[i] = 1.0;
  }

  void setup();

  double time_step(const SimulationOptions& options);
  void compute_artificial_viscosity(std::vector<double>& fv);
  void stabilize_pressure_gradient(const std::vector<double>& h,
                                   std::vector<double>& dh);

  void print_header(int n_bars = 108) const;
  void save(const std::string& filename) const;

  double total_area() const;
  double total_mass() const;
  double total_momentum() const;
  double total_energy() const;

  const auto& height() const { return height_; }

  void save_json(const std::string& filename) const;
  const auto& statistics() const { return statistics_; }
  auto& statistics() { return statistics_; }

 private:
  const Domain_t& domain_;
  Timer timer_;
  std::vector<double> height_;
  EarthProperties earth_;
  std::vector<double> volume_;
  const ShallowWaterOptions& options_;

  double initial_mass_;
  double initial_momentum_;
  double initial_energy_;
  ShallowWaterStatistics statistics_;

  spmat<double> h_mat_;
  vecd<double> h_rhs_;
  vecd<double> h_sol_;
};

struct WilliamsonCase1 : ShallowWaterOptions {
  WilliamsonCase1();
};

struct WilliamsonCase2 : ShallowWaterOptions {
  WilliamsonCase2();
};

struct WilliamsonCase5 : ShallowWaterOptions {
  WilliamsonCase5();
};

struct WilliamsonCase6 : ShallowWaterOptions {
  WilliamsonCase6();
};

void run_swe_simulation(const argparse::ArgumentParser& program);

}  // namespace vortex