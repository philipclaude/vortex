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

using namespace vortex;

UT_TEST_SUITE(optimaltransportgradient_test_suite)

UT_TEST_CASE(test_optimaltransportgradient)
{
    int n_iter = 50;
    size_t n_sites = 5000;
    int neighbors = 100;
    bool output_converge = false;

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

    std::vector<index_t> order(n_sites);
    sort_points_on_zcurve(sites.data(), n_sites, dim, order);

    Vertices vertices(dim);
    vertices.reserve(n_sites);
    coord_t x[dim];
    // add sites in z-curve order to Vertices list
    for (size_t i = 0; i < n_sites; i++)
    {
        for (int d = 0; d < dim; d++)
            x[d] = sites[dim * order[i] + d];
        vertices.add(x);
    }

    SphereDomain domain;
    VoronoiDiagram voronoi(dim, vertices[0], n_sites);
    VoronoiDiagramOptions options;
    options.n_neighbors = neighbors;
    options.allow_reattempt = false;
    options.parallel = true;

    auto &weights = voronoi.weights();
    weights.resize(n_sites, 0.0);

    std::vector<double> cell_size(n_sites, 4 * M_PI / n_sites);

    for (int iter = 1; iter <= n_iter; ++iter)
    {
        options.store_mesh = iter == n_iter;
        options.verbose = (iter == 1 || iter == n_iter - 1);
        voronoi.vertices().clear();
        voronoi.polygons().clear();
        voronoi.triangles().clear();
        voronoi.compute(domain, options);

        // move each site to the centroid of the corresponding cell
        voronoi.smooth(vertices, true);
    }
    // start weights as 0
    for (size_t k = 0; k < n_sites; k++)
    {
        weights[k] = 0;
    }

    double error = 1.0;
    double delta = 0.55;
    const double tol = 1e-8;
    int iter = 0;
    int max_iter = 10000;

    std::string hyphen = "_";
    std::string file_path = "../../data_test/gradient/runtime/cg" + hyphen + std::to_string(n_sites) + hyphen + std::to_string(delta) + ".txt";
    std::ofstream outputFile(file_path);

    Timer timer;
    timer.start();

    while (error >= tol && iter < max_iter)
    {
        iter++;

        std::vector<double> de_dw(n_sites);

        lift_sites(vertices, voronoi.weights());
        domain.set_initialization_fraction(0.7);
        voronoi.vertices().clear();
        voronoi.polygons().clear();
        voronoi.triangles().clear();
        voronoi.compute(domain, options);

        for (size_t k = 0; k < n_sites; k++)
        {
            int group_index = voronoi.polygons().group(k);
            de_dw[group_index] = voronoi.properties()[group_index].mass - cell_size[group_index];
            weights[group_index] = weights[group_index] - delta * (de_dw[group_index]);
        }

        error = calc_gradient_norm(de_dw);

        if (output_converge)
        {
            if (outputFile.is_open())
            {
                outputFile << "iter: " << iter << " error: " << error << std::endl;
            }
            else
            {
                break;
            }
        }
    }

    timer.stop();

    if (!output_converge)
    {
        if (outputFile.is_open())
        {
            outputFile << "Number Sites: " << n_sites << " Neighbors: " << neighbors << std::endl;
            outputFile << "Time: " << timer.seconds() << std::endl;
            outputFile << "Iterations: " << iter << " Error: " << error << std::endl;
        }
        else
        {
            std::cout << "Error opening file" << std::endl;
        }
    }
    outputFile.close();
}
UT_TEST_CASE_END(test_optimaltransportgradient)
UT_TEST_SUITE_END(optimaltransportgradient_test_suite)