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
#include "math/vec.hpp"
#include "tester.h"
#include "util.h"
#include <iostream>
#include <fstream>

using namespace vortex;

UT_TEST_SUITE(optimaltransport_test_suite)

double calc_triangle_area(vec3d &p0, vec3d &p1, vec3d &p2)
{
    vec3d edge1 = p1 - p0;
    vec3d edge2 = p2 - p0;

    vec3d crossProduct = cross(edge1, edge2);
    return 0.5 * length(crossProduct);
};

double calc_energy(VoronoiDiagram &voronoi, std::vector<double> &weights, double &idealWeight)
{
    const Topology<Polygon> &polygons = voronoi.polygons();

    double energy = 0.0;
    double second_sum = 0.0;
    for (int k = 0; k < polygons.n(); k++)
    {
        double cell_area = 0.0;
        // get the site of the polygon (assumes site's within polygon)
        int group_index = polygons.group(k);
        vec3d site(voronoi.vertices().group(group_index));

        int p0_index = polygons(k, 0);
        vec3d p0(voronoi.vertices()[p0_index]);

        for (int j = 1; j < voronoi.polygons().length(k) - 1; j++)
        {
            // get two vertices of the polygon
            int p1_index = polygons(k, j);
            int p2_index = polygons(k, j + 1);
            vec3d p1(voronoi.vertices()[p1_index]);
            vec3d p2(voronoi.vertices()[p2_index]);

            // calculate area of the triangle
            double t_area = calc_triangle_area(p0, p1, p2);

            // taken directly from Levy Geogram
            double cur_f = 0.0;
            for (index_t c = 0; c < 3; c++)
            {
                double u0 = site[c] - p0[c];
                double u1 = site[c] - p1[c];
                double u2 = site[c] - p2[c];
                cur_f += u0 * u0;
                cur_f += u1 * (u0 + u1);
                cur_f += u2 * (u0 + u1 + u2);
            }

            cell_area += t_area;
            energy += t_area * cur_f / 6.0;
        }
        energy += (-weights[k] * cell_area);
        second_sum += idealWeight * weights[k];
    }
    return energy + second_sum;
};

UT_TEST_CASE(test_optimaltransport)
{
    int n_iter = 10;
    double delta = 10e-3;
    double cell_size_tol = 10e-6;
    size_t n_sites = 10000;

    auto irand = [](int min, int max)
    {
        return min + double(rand()) / (double(RAND_MAX) + 1.0) * (max - min);
    };
    const double tol = 1e-12;
    static const int dim = 4; // number of dimensions

    double cell_size = 4 * M_PI / n_sites;

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
    options.n_neighbors = 75;
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
        auto props = voronoi.analyze();
        LOG << fmt::format("area = {}", props.area);
        LOG << fmt::format("iter = {}, area = {}", iter, calc_energy(voronoi, weights, cell_size));
    }
    auto props = voronoi.analyze();

    // // start weights as 0
    // for (size_t k = 0; k < n_sites; k++)
    // {
    //     weights[k] = 0;
    // }

    // // lift vertices to 4D
    // double error;
    // std::ofstream outFile;
    // outFile.open("../../data_test/hundred_convergence_error.txt");

    // do
    // {
    //     double sum = 0.0;
    //     lift_sites(vertices, voronoi.weights());

    //     domain.set_initialization_fraction(0.7);
    //     voronoi.vertices().clear();
    //     voronoi.polygons().clear();
    //     voronoi.triangles().clear();
    //     voronoi.compute(domain, options);

    //     for (size_t k = 0; k < n_sites; k++)
    //     {
    //         weights[k] = weights[k] - delta * (voronoi.properties()[k].mass - cell_size);
    //         sum = sum + pow((voronoi.properties()[k].mass - cell_size), 2);
    //     }
    //     error = pow((sum / n_sites), 0.5);
    //     outFile << std::to_string(error) << std::endl;
    // } while (error > cell_size_tol);

    // outFile.close();
    // auto props = voronoi.analyze();
    // UT_ASSERT_NEAR(props.area, 4 * M_PI, tol);

    // voronoi.merge();

    // LOG << fmt::format("writing {} polygons", voronoi.polygons().n());
    // if (voronoi.polygons().n() > 0)
    // {
    //     meshb::write(voronoi, "hundred_power_diagram.meshb");
    // }
}
UT_TEST_CASE_END(test_optimaltransport)

UT_TEST_SUITE_END(optimaltransport_test_suite)