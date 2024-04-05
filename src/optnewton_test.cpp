//
//  vortex: Semi-discreet Optimal Transport Problem on a Sphere
//
//  Copyright 2023 - 2024 Philip Claude Caplan, Otis Milliken
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
#include "voronoi.h"

#include <Predicates_psm.h>

#include <cmath>

#include "graphics.h"
#include "io.h"
#include "library.h"
#include <nlopt.hpp>
#include "math/vec.hpp"
#include "tester.h"
#include "util.h"
#include <iostream>
#include <fstream>
#include "optimaltransport.h"
#include "spmat.h"

using namespace vortex;

UT_TEST_SUITE(optimaltransportnewtonsphere_test_suite)

UT_TEST_CASE(test_optimaltransportnewtonsphere)
{
    int n_iter = 10;
    size_t n_sites = 1000;

    auto irand = [](int min, int max)
    {
        return min + double(rand()) / (double(RAND_MAX) + 1.0) * (max - min);
    };
    static const int dim = 4; // number of dimensions

    std::vector<coord_t> sites(n_sites * dim, 0.0);

    // generate random sites on the sphere
    for (size_t k = 0; k < n_sites; k++)
    {
        coord_t theta = 2.0 * M_PI * irand(0, 1);
        coord_t phi = acos(2.0 * irand(0, 1) - 1.0);
        sites[k * dim + 0] = cos(theta) * sin(phi);
        sites[k * dim + 1] = sin(theta) * sin(phi);
        sites[k * dim + 2] = cos(phi);
        if (dim > 3)
            sites[k * dim + 3] = 0.0;
    }

    // ideal cell size
    std::vector<double> cell_sizes(n_sites, 4 * M_PI / n_sites);

    std::vector<index_t> order(n_sites);
    sort_points_on_zcurve(sites.data(), n_sites, dim, order);

    Vertices vertices(dim);
    vertices.reserve(n_sites);
    coord_t vertex[dim];
    // add sites in z-curve order to Vertices list
    for (size_t i = 0; i < n_sites; i++)
    {
        for (int d = 0; d < dim; d++)
            vertex[d] = sites[dim * order[i] + d];
        vertices.add(vertex);
    }

    SphereDomain domain;
    VoronoiDiagram voronoi(dim, vertices[0], n_sites);
    VoronoiDiagramOptions options;
    options.n_neighbors = 75;
    options.allow_reattempt = false;
    options.parallel = true;

    for (int iter = 1; iter <= n_iter; ++iter)
    {
        options.store_mesh = true;
        options.verbose = (iter == 1 || iter == n_iter - 1);
        voronoi.vertices().clear();
        voronoi.polygons().clear();
        voronoi.triangles().clear();
        voronoi.compute(domain, options);

        // move each site to the centroid of the corresponding cell
        voronoi.smooth(vertices, true);
    }

    auto &weights = voronoi.weights();
    weights.resize(n_sites, 0.0);

    for (size_t k = 0; k < n_sites; k++)
    {
        weights[k] = 0;
    }

    spmat<double> hessian(voronoi.polygons().n(), voronoi.polygons().n());
    std::vector<double> de_dw(voronoi.polygons().n());
    std::vector<double> cell_size(n_sites, 4 * M_PI / n_sites);
    int max_iter = 10000;
    int iter = 0;
    while (iter < max_iter)
    {
        lift_sites(vertices, voronoi.weights());

        domain.set_initialization_fraction(0.7);
        voronoi.vertices().clear();
        voronoi.polygons().clear();
        voronoi.triangles().clear();
        voronoi.compute(domain, options);

        voronoi.merge();

        build_hessian(voronoi, hessian, de_dw, cell_size);
        vecd<double> search_direction(voronoi.polygons().n());

        hessian.solve_nl(de_dw, search_direction);

        for (size_t k = 0; k < n_sites; k++)
        {
            weights[k] = weights[k] - search_direction[k];
        }
    }
}
UT_TEST_CASE_END(test_optimaltransportnewtonsphere)

UT_TEST_SUITE_END(optimaltransportnewtonsphere_test_suite)