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

UT_TEST_SUITE(optimaltransportsphere_test_suite)

UT_TEST_CASE(test_optimaltransportsphere)
{
    int n_iter = 50;
    size_t n_sites = 5000;
    int neighbors = 100;
    bool converge = true;

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
    options.n_neighbors = neighbors;
    options.allow_reattempt = false;
    options.parallel = true;
    auto &weights = voronoi.weights();
    weights.resize(n_sites, 0.0);

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

    std::string hyphen = "_";
    std::string file_path = "";
    if (converge)
    {
        file_path = "../../data_test/quasi/energy/csqq" + hyphen + std::to_string(n_sites) + hyphen + std::to_string(neighbors) + ".txt";
    }
    else
    {
        file_path = "../../data_test/quasi/energy/runtime/csqq" + hyphen + std::to_string(n_sites) + hyphen + std::to_string(neighbors) + ".txt";
    }
    std::ofstream outputFile(file_path);

    double energy_time = 0.0;
    Timer timer;
    timer.start();
    nlopt_data<SphereDomain> data = {voronoi, domain, vertices, options, cell_sizes, 0, 0.0, 0.0, converge, outputFile, energy_time};
    std::vector<double> x(n_sites, 0.0);

    nlopt::opt opt(nlopt::LD_LBFGS, n_sites);

    opt.set_min_objective(objective_func<SphereDomain>, static_cast<void *>(&data));

    // set some optimization parameters
    opt.set_xtol_rel(0);
    opt.set_ftol_abs(1e-16);
    // opt.set_stopval(1e-8);
    opt.set_maxeval(10000);

    // set the lower and upper bounds on the weights
    std::vector<double> lower_bound(n_sites, 0.0);
    opt.set_lower_bounds(lower_bound);

    double f_opt;
    try
    {
        auto result = opt.optimize(x, f_opt);
        timer.stop();
        if (!converge)
        {
            if (outputFile.is_open())
            {
                outputFile << "Number Sites: " << n_sites << " Neighbors: " << neighbors << std::endl;
                outputFile << "Success" << std::endl;
                outputFile << "Time: " << timer.seconds() << std::endl;
                outputFile << "Energy: " << f_opt << std::endl;
                outputFile << "Average Energy Calculation: " << (energy_time / (double)data.iter) << std::endl;
                outputFile << "Error: " << data.error << " Iterations: " << data.iter << std::endl;
                outputFile << "result: " << result << std::endl;
            }
            else
            {
                std::cout << "Error opening file" << std::endl;
            }
        }
    }
    catch (std::exception &e)
    {
        timer.stop();
        LOG << f_opt;

        if (!converge)
        {
            if (outputFile.is_open())
            {
                outputFile << "Number Sites: " << n_sites << " Neighbors: " << neighbors << std::endl;
                outputFile << "Fail" << std::endl;
                outputFile << "Time: " << timer.seconds() << std::endl;
                outputFile << "Average Energy Calculation: " << (energy_time / (double)data.iter) << std::endl;
                outputFile << "Error: " << data.error << "Iterations: " << data.iter << std::endl;
            }
            else
            {
                std::cout << "Error opening file" << std::endl;
            }
        }
        std::cout << e.what() << std::endl;
    }
    outputFile.close();
}
UT_TEST_CASE_END(test_optimaltransportsphere)

UT_TEST_SUITE_END(optimaltransportsphere_test_suite)