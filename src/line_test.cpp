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

UT_TEST_SUITE(line_test_suite)

// Test to see the order of the line mesh
UT_TEST_CASE(water_mesh) {
    Mesh mesh(3);
    meshb::read("../build/release/water.meshb", mesh);

    //size_t l = mesh.lines().n();
    
    for(int i = 5159; i < 5257; i++){
        index_t L1F = mesh.lines()[i][0]; 
        //index_t L1S = mesh.lines()[i][1];
        //std::cout << i << ": (" << mesh.vertices()[L1F][0] << ", " << mesh.vertices()[L1F][1] << ") & (" << mesh.vertices()[L1S][0] << ", " << mesh.vertices()[L1S][1] << ")" << std::endl;
        //std::cout << i << ": (" << mesh.vertices()[L1F][0] << ", " << mesh.vertices()[L1F][1] << ", " << mesh.vertices()[L1F][2] << ") & (" << mesh.vertices()[L1S][0] << ", " << mesh.vertices()[L1S][1] << ", " << mesh.vertices()[L1S][2] << ")" << std::endl;
        std::cout << "(" << mesh.vertices()[L1F][0] << ", " << mesh.vertices()[L1F][1] << ", " << mesh.vertices()[L1F][2] << ")" << std::endl;
    }

    /*
    index_t first = mesh.lines()[5356][0];
    //2690, 3091, 4343, 4916, 5159, 5257, 5356, 6021

    index_t c = 5356;
    index_t second = mesh.lines()[c][1];
    while(first != second){
        second = mesh.lines()[c][1];
        c++;
    }
    std::cout << c;
    */
}
UT_TEST_CASE_END(water_mesh)


UT_TEST_SUITE_END(line_test_suite)