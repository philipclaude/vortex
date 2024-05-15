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

UT_TEST_SUITE(plane_test_suite)

// Returns true if point is on the right side of the plane between edgeStart and edgeEnd (edge runs upward)
int planeSide(const coord_t* edgeStart, const coord_t* edgeEnd, const coord_t* point) {

    vec3d e1 = edgeStart; vec3d e2 = edgeEnd;
    vec3d e = e2 - e1;
    vec3d mid = (e1 + e2)/2.0;

    vec3d ns = normalize(mid);

    vec3d nplane = normalize(cross(e,ns));

    vec3d p = point;
    double r = dot(nplane,(p-mid));
    
    if(r < 0){
        return 0;
    }else if(r > 0){
        return 1;
    }else{
        return 2;
    }
}

// Find the intersection point between two 3D line segments (geodesics) on a unit sphere
bool findIntersection(const coord_t* edge1Start, const coord_t* edge1End, const coord_t* edge2Start, const coord_t* edge2End, vec3d& x) {
    vec3d e1S = edge1Start; vec3d e1E = edge1End;vec3d e2S = edge2Start; vec3d e2E = edge2End;
    // Create plane of intersecting edge
    vec3d i = e1E - e1S;
    vec3d mid = (e1S + e1E)/2.0;
    vec3d ns = normalize(mid);
    vec3d nplane = normalize(cross(i,ns));

    vec3d e = e2E - e2S;
    double t = dot(nplane,(mid-e2S)) / dot(nplane,e);

    vec3d r = e2S + t*e;

    x = normalize(r);

    const double epsilon = 1e-9;

    if(t < 0.0 - epsilon || t > 1.0 + epsilon){
        return false;
    }else{
        return true;
    }
}

UT_TEST_CASE(plane_side) {

    coord_t a[3] = {-0.11672606,-0.6245624,0.7722026};
    coord_t b[3] = {-0.11159379,-0.84252536,0.52697057};
    coord_t c[3] = {0.1848934,-0.6997241,0.6900728};
    coord_t d[3] = {-0.36002824,-0.7606625,0.5401596};

    //coord_t e[3] = {0.23325026, 0.8491885, 0.47378588};

    //bool x = planeSide(a,b,c)==planeSide(a,b,d);
    //int y = planeSide(a,b,c);

    //bool t = insideCheck(a,b,c);
    //std::cout << t << std::endl;

    //std::cout << x << std::endl;
    //std::cout << y << std::endl;

    vec3d r;
    bool s = findIntersection(a,b,c,d,r);

    std::cout << "(" << r[0] << ", " << r[1] << ", " << r[2] << ")" << std::endl;
    LOG << fmt::format("s = {}",s);
}
UT_TEST_CASE_END(plane_side)


UT_TEST_SUITE_END(plane_test_suite)