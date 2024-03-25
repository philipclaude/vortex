//
//  vortex: Voronoi mesher and fluid simulator for the Earth's oceans and
//  atmosphere.
//
//  Copyright 2023 - 2024 Philip Claude Caplan & Tobias W Pouler
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

#include <Predicates_psm.h>

#include <cmath>

#include "graphics.h"
#include "halfedges.h"
#include "io.h"
#include "library.h"
#include "log.h"
#include "mesh.h"
#include "mesher.h"
#include "numerics.h"
#include "texture.h"
#include "util.h"
#include "voronoi.h"
#include "tester.h"

#ifndef VORTEX_SOURCE_DIR
#define VORTEX_SOURCE_DIR "./"
#endif

using namespace vortex;

UT_TEST_SUITE(oceans_test_suite)

UT_TEST_CASE(test_random_oceans) {
    auto irand = [](int min, int max) {
        return min + double(rand()) / (double(RAND_MAX) + 1.0) * (max - min);
    };
    size_t n_points = 1000000;
    int n_smooth = 10;
    int count = 0;
    Vertices sample(3);

    sample.reserve(n_points);
    std::string tex_file =
        std::string(VORTEX_SOURCE_DIR) + "/../data/oceans_2048.png";
    TextureOptions tex_opts;
    tex_opts.format = TextureFormat::kGrayscale;
    Texture texture(tex_file, tex_opts);
    texture.make_binary(10, 10, 255);

    vec3d xx, uv;
    while (sample.n() < n_points) {
      coord_t theta = 2.0 * M_PI * irand(0, 1);
      coord_t phi = acos(2.0 * irand(0, 1) - 1.0);
      xx[0] = cos(theta) * sin(phi);
      xx[1] = sin(theta) * sin(phi);
      xx[2] = cos(phi);

      // calculate (u, v) consistent with other algorithms
      sphere_params(xx, uv);
      double t;
      texture.sample(uv[0], uv[1], &t);
      if (t < 50){ 
        count++;
        continue;
      }
      sample.add(&xx[0]);
    }
    LOG << fmt::format("land = {}, total = {}", count, n_points);

    int dim  = 3;

    n_points = sample.n();
    std::vector<index_t> order(n_points);
    sort_points_on_zcurve(sample[0], n_points, dim, order);
    Vertices points(dim);
    points.reserve(n_points);
    coord_t x[dim];
    for (size_t i = 0; i < n_points; i++) {
        for (int d = 0; d < dim; d++) x[d] = sample[order[i]][d];
        points.add(x);
    }

    VoronoiDiagram voronoi(dim, points[0], n_points);
    VoronoiDiagramOptions options;
    SphereDomain domain;

    options.n_neighbors = 75;
    options.allow_reattempt = false;
    options.parallel = true;

    auto& weights = voronoi.weights();
    weights.resize(n_points, 0.0);

    int n_iter = n_smooth;
    for (int iter = 1; iter <= n_iter; ++iter) {
      options.store_mesh = iter == n_iter;
      options.verbose = (iter == 1 || iter == n_iter - 1);
      voronoi.vertices().clear();
      //voronoi.vertices().set_dim(3);
      voronoi.polygons().clear();
      voronoi.triangles().clear();
      voronoi.compute(domain, options);

      // move each site to the centroid of the corresponding cell
      voronoi.smooth(points, true);
      auto props = voronoi.analyze();
      LOG << fmt::format("iter = {}, area = {}", iter, props.area);
    }

    int test = 0;
    for(int j = 0; j < n_points; j++){
        vec3d p, q;
        p[0] = points[j][0];
        p[1] = points[j][1];
        p[2] = points[j][2];
        sphere_params(p, q);
        double t;
        texture.sample(q[0], q[1], &t);
        if (t < 50){ 
            test++;
        }
    }
    LOG << fmt::format("land_points = {}, total = {}", test, n_points);
}
UT_TEST_CASE_END(test_random_oceans)

UT_TEST_SUITE_END(oceans_test_suite)