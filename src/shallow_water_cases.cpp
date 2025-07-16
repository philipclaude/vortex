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

namespace vortex {

namespace {
void xyz_to_latlon(const double* x, double& lat, double& lon) {
  lon = atan2(x[1], x[0]);  // atan2 returns in [-pi, pi]
  double z = x[2];
  ASSERT(std::fabs(x[2]) <= 1);
  lat = asin(z);  // asin returns in [-pi/2, pi/2]
}
}  // namespace

WilliamsonCase1::WilliamsonCase1() {
  const double omega = earth.angular_velocity;
  const double a = earth.radius;
  const double R = a / 3;
  const double h0 = 1000;
  const double lc = 3 * M_PI / 2;
  const double tc = 0;
  const double twelve_days = 12 * 24 * 3600;
  const double u0 = 2 * M_PI * a / twelve_days;
  const double hm = 1000;

  surface_height = [](const double* x) -> double { return 0.0; };
  initial_height = [a, R, h0, lc, tc, hm](const double* x) -> double {
    double lambda, theta;
    xyz_to_latlon(x, theta, lambda);

    double dl = lambda - lc;
    // double dt = theta - tc;
    double r =
        a * std::acos(sin(tc) * sin(theta) + cos(tc) * cos(theta) * cos(dl));
    if (r >= R) return hm;
    return 0.5 * h0 * (1 + cos(M_PI * r / R)) + hm;
  };
  initial_velocity = [u0](const double* x) -> vec3d {
    double lambda, theta;
    xyz_to_latlon(x, theta, lambda);
    double us = u0 * sin(theta) * cos(lambda);
    double vs = -u0 * sin(lambda);
    double ux = -us * sin(lambda) - vs * sin(theta) * cos(lambda);
    double uy = us * cos(lambda) - vs * sin(theta) * sin(lambda);
    double uz = vs * cos(theta);
    return {ux, uy, uz};
  };
  coriolis_parameter = [&](const double* x) -> double {
    return 2.0 * omega * x[2];
  };
  analytic_velocity = [&](const double* x, double time) -> vec3d {
    return initial_velocity(x);
  };
  use_analytic_velocity = true;

  analytic_height = [this, twelve_days](const double* x,
                                        double time) -> double {
    if (time != twelve_days) return 1e20;
    return initial_height(x);
  };
  has_analytic_height = true;
  time_stepping = TimeSteppingScheme::kExplicit;
  days = 12;
}

WilliamsonCase2::WilliamsonCase2() {
  const double omega = earth.angular_velocity;
  const double a = earth.radius;
  const double g = earth.gravity;

  double h0 = 2.94e4 / g;
  double u0 = 2 * M_PI * a / (12 * 24 * 3600);
  surface_height = [](const double* x) -> double { return 0.0; };

  initial_velocity = [u0](const double* x) -> vec3d {
    double lambda, theta;
    xyz_to_latlon(x, theta, lambda);
    double us = u0 * sin(theta) * cos(lambda);
    double vs = -u0 * sin(lambda);
    double ux = -us * sin(lambda) - vs * sin(theta) * cos(lambda);
    double uy = us * cos(lambda) - vs * sin(theta) * sin(lambda);
    double uz = vs * cos(theta);
    return {ux, uy, uz};
  };

  coriolis_parameter = [omega](const double* x) -> double {
    return -2.0 * omega * x[0];
  };

  analytic_height = [h0, omega, a, u0, g](const double* x,
                                          double time) -> double {
    return h0 - (omega * a * u0 + 0.5 * u0 * u0) * x[0] * x[0] / g;
  };
  initial_height = [this](const double* x) -> double {
    return analytic_height(x, 0);
  };
  has_analytic_height = true;
  has_analytic_velocity = true;
  analytic_velocity = [this](const double* x, double time) -> vec3d {
    return initial_velocity(x);
  };
  days = 12;
}

WilliamsonCase5::WilliamsonCase5() {
  const double omega = earth.angular_velocity;
  const double a = earth.radius;
  const double g = earth.gravity;

  const double u0 = 20.0;
  initial_velocity = [u0](const double* x) -> vec3d {
    return {-u0 * x[1], u0 * x[0], 0.0};
  };

  const double h0 = 5960;
  const double hs0 = 2000;
  const double R = M_PI / 9.0;
  const double lambda_c = -M_PI / 2;
  const double theta_c = M_PI / 6;

  surface_height = [hs0, R, lambda_c, theta_c](const double* x) -> double {
    double lambda, theta;
    xyz_to_latlon(x, theta, lambda);
    double d_lambda = lambda - lambda_c;
    double d_theta = theta - theta_c;
    double d = std::min(R * R, d_lambda * d_lambda + d_theta * d_theta);
    double hs = hs0 * (1 - std::pow(d, 0.5) / R);
    // return hs0 * std::exp(-2.8 * 2.8 * d / (R * R));
    return hs;
  };
  initial_height = [a, omega, g, u0, h0](const double* x) -> double {
    double h = h0 - (omega * a * u0 + 0.5 * u0 * u0) * x[2] * x[2] / g;
    return h;
  };

  coriolis_parameter = [omega](const double* x) -> double {
    return 2.0 * omega * x[2];
  };
  double hmin = h0 - (omega * a * u0 + 0.5 * u0 * u0) / g;
  LOG << fmt::format("hmin = {}, hmax = {}", hmin, h0);
  days = 15;
}

WilliamsonCase6::WilliamsonCase6() {
  const double omega = earth.angular_velocity;
  const double a = earth.radius;
  const double g = earth.gravity;

  const double w = 7.848e-6;
  const double K = 7.848e-6;
  const double h0 = 8000;
  double m = 4;

  auto A = [w, omega, K, m](double t) -> double {
    // the last term is expanded to avoid dividing by cos(t)^2
    double ct = cos(t);
    double ct2 = ct * ct;
    return 0.5 * w * (2 * omega + w) * ct2 +
           0.25 * K * K * pow(ct, 2 * m) *
               ((m + 1) * ct2 + (2 * m * m - m - 2)) -
           0.5 * K * K * m * m * pow(ct, 2 * (m - 1));
  };

  auto B = [w, omega, K, m](double t) -> double {
    double ct = cos(t);
    double ct2 = ct * ct;
    return 2 * (omega + w) * K * pow(ct, m) *
           ((m * m + 2 * m + 2) - pow(m + 1, 2) * ct2) / ((m + 1) * (m + 2));
  };
  auto C = [K, m](double t) -> double {
    double ct = cos(t);
    double ct2 = ct * ct;
    return 0.25 * K * K * pow(ct, 2 * m) * ((m + 1) * ct2 - (m + 2));
  };

  initial_velocity = [a, w, K, m](const double* x) -> vec3d {
    double l, t;
    xyz_to_latlon(x, t, l);
    double ct = cos(t);
    double st = sin(t);
    double ct2 = ct * ct;
    double st2 = st * st;
    double us =
        a * w * ct + a * K * pow(ct, m - 1) * (m * st2 - ct2) * cos(m * l);
    double vs = -a * K * m * pow(ct, m - 1) * st * sin(m * l);
    double ux = -us * sin(l) - vs * st * cos(l);
    double uy = us * cos(l) - vs * st * sin(l);
    double uz = vs * ct;
    return {ux, uy, uz};
  };

  initial_height = [A, B, C, h0, a, m, g](const double* x) -> double {
    double l, t;
    xyz_to_latlon(x, t, l);
    return h0 + a * a * (A(t) + B(t) * cos(m * l) + C(t) * cos(2 * m * l)) / g;
  };

  surface_height = [](const double* x) -> double { return 0.0; };
  coriolis_parameter = [omega](const double* x) -> double {
    return 2.0 * omega * x[2];
  };

  double nu = (m * (3 + m) * w - 2 * omega) / ((1 + m) * (2 + m));
  analytic_height = [A, B, C, h0, a, m, g, nu](const double* x,
                                               double time) -> double {
    double l, t;
    xyz_to_latlon(x, t, l);
    l -= nu * time;
    return h0 + a * a * (A(t) + B(t) * cos(m * l) + C(t) * cos(2 * m * l)) / g;
  };
  has_analytic_height = true;
  days = 15;
}

GalewskyCase::GalewskyCase() {
  const double omega = earth.angular_velocity;
  const double a = earth.radius;
  const double g = earth.gravity;

  const double umax = 80.0;
  const double phi0 = M_PI / 7.0;
  const double phi1 = M_PI / 2.0 - M_PI / 7.0;
  const double en = std::exp(-4.0 / std::pow(phi1 - phi0, 2.0));
  const double u0 = 20.0;
  const double h0 = 5960;

  surface_height = [](const double* x) -> double { return 0.0; };

  coriolis_parameter = [omega](const double* x) -> double {
    return 2.0 * omega * x[2];
  };

  // This initial height is from Williamson, test case 5
  initial_height = [a, omega, g, u0, h0](const double* x) -> double {
    double h = h0 - (omega * a * u0 + 0.5 * u0 * u0) * x[2] * x[2] / g;
    return h;
  };

  initial_velocity = [umax, phi0, phi1, en](const double* x) -> vec3d {
    double lat, lon;
    xyz_to_latlon(x, lat, lon);

    double u_zonal = 0.0;
    if (lat >= phi0 && lat <= phi1) {
      double denom = (lat - phi0) * (lat - phi1);
      u_zonal = umax * std::exp(1.0 / denom) / en;
    }

    double us = u_zonal;
    double vs = 0.0;
    double ux = -us * sin(lon) - vs * sin(lat) * cos(lon);
    double uy = us * cos(lon) - vs * sin(lat) * sin(lon);
    double uz = vs * cos(lat);

    return {ux, uy, uz};
  };
  days = 6;
}

}  // namespace vortex
