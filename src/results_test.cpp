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
// This file was used to test polygon counts and surface area values for all the Voronoi diagrams that I generated. (throw away)
//
#include <fmt/format.h>

#include <argparse/argparse.hpp>
#include <queue>
#include <unordered_set>

#include <Predicates_psm.h>

#include <cmath>

#include "elements.h"
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
#include "tester.h"
#include "stlext.h"
#include "triangulate.h"
#include "trees/kdtree.h"
#include "voronoi.h"
#include "voronoi_polygon.hpp"

#include <cstdlib>

#ifndef VORTEX_SOURCE_DIR
#define VORTEX_SOURCE_DIR "./"
#endif

using namespace vortex;

UT_TEST_SUITE(results_test_suite)

double calculateArea(const coord_t* p1, const coord_t* p2, const coord_t* p3) {
    // Calculate vectors
    double v1x = p2[0] - p1[0];
    double v1y = p2[1] - p1[1];
    double v1z = p2[2] - p1[2];

    double v2x = p3[0] - p1[0];
    double v2y = p3[1] - p1[1];
    double v2z = p3[2] - p1[2];

    // Calculate cross product
    double crossX = v1y * v2z - v1z * v2y;
    double crossY = v1z * v2x - v1x * v2z;
    double crossZ = v1x * v2y - v1y * v2x;

    // Calculate magnitude of cross product
    double magnitude = sqrt(crossX * crossX + crossY * crossY + crossZ * crossZ);

    // Calculate area
    double area = magnitude / 2.0;

    return area;
}

UT_TEST_CASE(results){
// read in mesh to be edited
    
    
    static const int dim = 3;
    Mesh mesh(dim);
    meshb::read("../build/release/triangle_01_1M.meshb", mesh);
    //meshb::read("../build/release/example6.meshb", mesh);
    //meshb::read("../data/A_10M.meshb", mesh);
    // size_t n_sites = mesh.polygons().n();
    size_t n_sites = 0;
    for(int p = 0; p < mesh.polygons().n(); p++){
        if(mesh.polygons().group(p) > n_sites){
            n_sites = mesh.polygons().group(p);
        }
    }
    std::cout << n_sites;
    
    
    //system("../build/release/bin/vortex mesh ../build/release/../../data/oceans_2048.png --hmin 0.005 --output earth.meshb");
    //system("../build/release/bin/vortex extract ../build/release/earth.meshb --oceans water.meshb --continents land.meshb");
    //system("../build/release/bin/vortex voronoi --domain ../build/release/water.meshb --points random --n_points 10000 --n_smooth 10 --output example6.meshb");

    
    PolygonTriangulation triangulation(mesh.vertices(), mesh.polygons());
    triangulation.triangulate(TangentSpaceType::kSphere, 0, mesh.polygons().n());

    mesh.triangles().clear();
    for (size_t k = 0; k < triangulation.n(); k++) {
        mesh.triangles().add(triangulation.triangle(k));
        mesh.triangles().set_group(k, triangulation.group(k));
    }
    double area = 0.0;
    for (size_t k = 0; k < mesh.triangles().n(); k++) {
        area += calculateArea(mesh.vertices()[mesh.triangles()[k][0]], mesh.vertices()[mesh.triangles()[k][1]], mesh.vertices()[mesh.triangles()[k][2]]);
    }
    std::cout << area;
    
}
UT_TEST_CASE_END(results)

UT_TEST_SUITE_END(results_test_suite)