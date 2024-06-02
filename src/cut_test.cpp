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
// This was the main code file used to actually complete the boundary cutting process. Other than the handling of intersected polygons without 
// boundary vertices inside and the check for more than two (four) intersection points (Frank's issue), this file should contain re-usable 
// methods for Vortex. Please refer to the Methods section of my thesis for more information about the functions and clipping process.
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
#include "trees/kdtree.h"
#include "voronoi.h"
#include "voronoi_polygon.hpp"

#ifndef VORTEX_SOURCE_DIR
#define VORTEX_SOURCE_DIR "./"
#endif

using namespace vortex;

UT_TEST_SUITE(cut_test_suite)

// KD tree for nearest neighbors to a single point
template <int dim>
std::shared_ptr<trees::KdTreeNd<coord_t, index_t>> get_nearest_neighbor(
    const coord_t* p, uint64_t np, const coord_t* q, uint64_t nq,
    std::vector<index_t>& nn, std::shared_ptr<trees::KdTreeNd<coord_t, index_t>> ptree = nullptr) {
  Timer timer;
  trees::KdTreeOptions kdtree_opts;
  kdtree_opts.max_dim = dim;
  using kdtree_t = trees::KdTree<dim, coord_t, index_t>;
  if (!ptree) {
    timer.start();
    ptree = std::make_shared<kdtree_t>(p, np, kdtree_opts);
    timer.stop();  
    LOG << "kdtree created in " << timer.seconds() << " s.";
  }

  auto* tree = static_cast<kdtree_t*>(ptree.get());
  timer.start();
  std::parafor_i(0, nq,
                 [&](int tid, int k) { nn[k] = tree->nearest(&q[k * dim]); });
  timer.stop();
  LOG << "nearest neighbors computed in " << timer.seconds() << " s.";
  return ptree;
}

// Returns 1 if point is on the right side of the plane between edgeStart and edgeEnd (edge runs upward), 0 if on left, and 2 otherwise
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

UT_TEST_CASE(voronoi_cut) {
    Timer timer;
    timer.start();
    
    // read in mesh to be edited
    static const int dim = 3;
    Mesh o_mesh(dim);
    meshb::read("../build/release/analytic_10M.meshb", o_mesh);
    //meshb::read("../build/release/example7.meshb", o_mesh);
    //meshb::read("../build/release/example4.meshb", o_mesh);
    size_t n_sites = o_mesh.polygons().n();
    
    // create boundary object (l_mesh) and place lines on original mesh (o_mesh)
    Mesh l_mesh(dim);
    meshb::read("../build/release/water_0025.meshb", l_mesh);
    index_t current = o_mesh.vertices().n();
    index_t l = 0;
    o_mesh.lines().reserve(l_mesh.lines().n());
    while(l < l_mesh.lines().n()){
        int g = l_mesh.lines().group(l); 
        index_t c = 0;
        int p = l_mesh.lines()[l][0];
        coord_t v[dim] = {l_mesh.vertices()[p][0],l_mesh.vertices()[p][1],l_mesh.vertices()[p][2]};
        o_mesh.vertices().add(v);
        while(l + 1 < l_mesh.lines().n() && g == l_mesh.lines().group(l+1)){
            p = l_mesh.lines()[l][1];
            v[0] = l_mesh.vertices()[p][0]; v[1] = l_mesh.vertices()[p][1]; v[2] = l_mesh.vertices()[p][2];
            o_mesh.vertices().add(v);
            int a = current+c; int b = current+c+1;
            int line[2] = {a,b};
            size_t id = o_mesh.lines().n();
            o_mesh.lines().add(line);
            o_mesh.lines().set_group(id, l_mesh.lines().group(l));
            l++; c++;
        }
        int a = current+c; int b = current;
        int line[2] = {a,b};
        size_t id = o_mesh.lines().n();
        o_mesh.lines().add(line);
        o_mesh.lines().set_group(id, l_mesh.lines().group(l));
        l++; current+=(c+1);
    }    

    // new mesh object
    Mesh mesh(dim);
    mesh.vertices().reserve(o_mesh.vertices().n());
    mesh.polygons().reserve(n_sites);
    
    // texture options
    std::string tex_file =
        std::string(VORTEX_SOURCE_DIR) + "/../data/oceans_2048.png";
    TextureOptions tex_opts;
    tex_opts.format = TextureFormat::kGrayscale;
    Texture texture(tex_file, tex_opts);
    texture.make_binary(10, 10, 255);
    timer.stop();
    LOG << fmt::format("All preprocessing in {} seconds", timer.seconds());

    timer.start();
    // KD tree
    size_t n_boundaries = o_mesh.lines().n();
    std::vector<index_t> VNN(o_mesh.lines().n());
    std::shared_ptr<trees::KdTreeNd<coord_t, index_t>> tree{nullptr};
    get_nearest_neighbor<dim>(o_mesh.vertices()[0], n_sites, o_mesh.vertices()[o_mesh.vertices().n()-n_boundaries], n_boundaries, VNN, tree);
    timer.stop();
    LOG << fmt::format("All neighbors calculated in {} seconds", timer.seconds());

    timer.start();
    // site to polygon to ensure correct polygon indices
    index_t S2P[n_sites];
    for(index_t s = 0; s < n_sites; s++){
        S2P[o_mesh.polygons().group(s)] = s;
    }

    // nearest polygon with parallel-sorted polygons
    //std::vector<index_t> VNP(n_boundaries);
    // adjacent vertices of boundary edges (use vector reserve [upper bound] and shrink to fit after loop) : future work
    std::vector<std::pair<index_t, index_t>> ABV(n_boundaries);
    // P2BV: polygon to boundary vertices
    std::vector<std::vector<coord_t>> P2BV(o_mesh.polygons().n());
    for (auto& vect : P2BV) vect.reserve(10);
    for(index_t vn = 0; vn < n_boundaries; vn++){
        //VNP[vn] = S2P[VNN[vn]];
        index_t curr = o_mesh.vertices().n() - n_boundaries + vn;
        if(vn == 0){
            ABV[vn] = std::make_pair(o_mesh.vertices().n()-1, curr+1);
        }else if(vn==n_boundaries-1){
            ABV[vn] = std::make_pair(curr-1,o_mesh.vertices().n()-n_boundaries);
        }else{
            ABV[vn] = std::make_pair(curr-1,curr+1);
        }
        P2BV[S2P[VNN[vn]]].push_back(o_mesh.vertices().n() - n_boundaries + vn);
    }

    // I was able to either remove these data structures or place them in already-existing loops, but I'll keep the code here
    /*
    // ABV: adjacent boundary vertices
    for(index_t b = 0; b < n_boundaries; b++){
        index_t curr = o_mesh.vertices().n() - n_boundaries + b;
        if(b == 0){
            ABV.push_back(std::make_pair(o_mesh.vertices().n()-1, curr+1));
        }else if(b==n_boundaries-1){
            ABV.push_back(std::make_pair(curr-1,o_mesh.vertices().n()-n_boundaries));
        }else{
            ABV.push_back(std::make_pair(curr-1,curr+1));
        }
    }
    */
   /*
    // P2BV: polygon to boundary vertices
    std::vector<std::vector<coord_t>> P2BV(o_mesh.polygons().n());
    for(index_t n = 0; n < n_boundaries; n++){
        P2BV[VNP[n]].push_back(o_mesh.vertices().n() - n_boundaries + n);
    }
    */

    // polygon to boundary edges
    std::vector<std::vector<std::pair<index_t, index_t>>> P2BE(o_mesh.polygons().n());
    for (auto& vect : P2BE) vect.reserve(10);
    for (index_t m = 0; m < P2BV.size(); m++) {
        if (!P2BV[m].empty()) {
            int count = 0;
            for (const auto& elem : P2BV[m]) {
                index_t curr = elem - o_mesh.vertices().n() + n_boundaries;
                if(count < 1){
                    P2BE[m].push_back(std::make_pair(ABV[curr].first, elem));
                }
                P2BE[m].push_back(std::make_pair(elem, ABV[curr].second));
                count++;
            }
        } 
    }
    timer.stop();
    LOG << fmt::format("All data structures made in {} seconds", timer.seconds());

    timer.start();
    // Creation of new restricted polygons
    for(index_t poly = 0; poly < o_mesh.polygons().n(); poly++){
        // remove inside cells (texture-based) - not 100% accurate
        bool water = false;
        for (size_t i = 0; i < o_mesh.polygons().length(poly); i++){
            vec3d p, q;
            p[0] = o_mesh.vertices()[o_mesh.polygons()[poly][i]][0];
            p[1] = o_mesh.vertices()[o_mesh.polygons()[poly][i]][1];
            p[2] = o_mesh.vertices()[o_mesh.polygons()[poly][i]][2];
            sphere_params(p, q);
            double t;
            texture.sample(q[0], q[1], &t);
            if (t >= 11){ 
                water = true;
                break;
            }
        }
        // if at least one vertex is in the water
        if(water){   
            index_t vertices_count = 0;
            if(!P2BE[poly].empty()){
                index_t v1 = 0;
                // find first intersection
                while(v1 < o_mesh.polygons().length(poly)){
                    index_t v2 = (v1 + 1)%o_mesh.polygons().length(poly);
                    if(planeSide(o_mesh.vertices()[o_mesh.polygons()[poly][v1]],o_mesh.vertices()[o_mesh.polygons()[poly][v2]],o_mesh.vertices()[P2BE[poly].front().first])!=planeSide(o_mesh.vertices()[o_mesh.polygons()[poly][v1]],o_mesh.vertices()[o_mesh.polygons()[poly][v2]],o_mesh.vertices()[P2BE[poly].front().second])){
                        vec3d intersection;
                        if(findIntersection(o_mesh.vertices()[P2BE[poly].front().first], o_mesh.vertices()[P2BE[poly].front().second], o_mesh.vertices()[o_mesh.polygons()[poly][v1]], o_mesh.vertices()[o_mesh.polygons()[poly][v2]],intersection)){
                            v1++;
                            coord_t ver[dim] = {intersection[0],intersection[1],intersection[2]};
                            mesh.vertices().add(ver);
                            vertices_count++;
                            break;
                        }
                    }
                    v1++;
                }
                index_t temp = 0;
                // loop through vertices until second intersection is found
                while(temp < o_mesh.polygons().length(poly)){
                    v1 = (v1)%o_mesh.polygons().length(poly);
                    index_t v2 = (v1 + 1)%o_mesh.polygons().length(poly);
                    coord_t ver[dim] = {o_mesh.vertices()[o_mesh.polygons()[poly][v1]][0],o_mesh.vertices()[o_mesh.polygons()[poly][v1]][1],o_mesh.vertices()[o_mesh.polygons()[poly][v1]][2]};
                    mesh.vertices().add(ver);
                    vertices_count++;
                    if(planeSide(o_mesh.vertices()[o_mesh.polygons()[poly][v1]],o_mesh.vertices()[o_mesh.polygons()[poly][v2]],o_mesh.vertices()[P2BE[poly].back().first])!=planeSide(o_mesh.vertices()[o_mesh.polygons()[poly][v1]],o_mesh.vertices()[o_mesh.polygons()[poly][v2]],o_mesh.vertices()[P2BE[poly].back().second])){
                        vec3d intersection;
                        if(findIntersection(o_mesh.vertices()[P2BE[poly].back().first], o_mesh.vertices()[P2BE[poly].back().second], o_mesh.vertices()[o_mesh.polygons()[poly][v1]], o_mesh.vertices()[o_mesh.polygons()[poly][v2]],intersection)){
                            v1++;
                            temp++;
                            coord_t ver[dim] = {intersection[0],intersection[1],intersection[2]};
                            mesh.vertices().add(ver);
                            vertices_count++;
                            break;
                        }
                    }
                    v1++;
                    temp++;
                }
                // add any remaining (interior) boundary vertices
                for(index_t bE = 1; bE < P2BE[poly].size(); bE++){
                    coord_t ver[dim] = {o_mesh.vertices()[P2BE[poly][P2BE[poly].size() - bE].first][0],o_mesh.vertices()[P2BE[poly][P2BE[poly].size() - bE].first][1],o_mesh.vertices()[P2BE[poly][P2BE[poly].size() - bE].first][2]};
                    mesh.vertices().add(ver);
                    vertices_count++;
                }
            // if not intersected, just add the original polygon
            }else{
                for(index_t v = 0; v < o_mesh.polygons().length(poly); v++){
                    coord_t ver[dim] = {o_mesh.vertices()[o_mesh.polygons()[poly][v]][0],o_mesh.vertices()[o_mesh.polygons()[poly][v]][1],o_mesh.vertices()[o_mesh.polygons()[poly][v]][2]};
                    mesh.vertices().add(ver);
                    vertices_count++;
                }
            }
            // add the polygon to the new mesh
            index_t new_polygon[vertices_count];
            for(index_t vc = 0; vc < vertices_count; vc++){
                new_polygon[vc] = mesh.vertices().n() - vertices_count + vc;
            }
            mesh.polygons().add(new_polygon,vertices_count);
        }
    }
    timer.stop();
    LOG << fmt::format("All polygons altered and added in {} seconds", timer.seconds());

    // add lines to the mesh for visualization purposes
    mesh.lines().reserve(n_boundaries);
    current = mesh.vertices().n();
    l = 0;
    while(l < l_mesh.lines().n()){
        int g = l_mesh.lines().group(l); 
        index_t c = 0;
        int p = l_mesh.lines()[l][0];
        coord_t v[dim] = {l_mesh.vertices()[p][0],l_mesh.vertices()[p][1],l_mesh.vertices()[p][2]};
        mesh.vertices().add(v);
        while(l + 1 < l_mesh.lines().n() && g == l_mesh.lines().group(l+1)){
            p = l_mesh.lines()[l][1];
            v[0] = l_mesh.vertices()[p][0]; v[1] = l_mesh.vertices()[p][1]; v[2] = l_mesh.vertices()[p][2];
            mesh.vertices().add(v);
            int a = current+c; int b = current+c+1;
            int line[2] = {a,b};
            size_t id = mesh.lines().n();
            mesh.lines().add(line);
            mesh.lines().set_group(id, l_mesh.lines().group(l));
            l++; c++;
        }
        int a = current+c; int b = current;
        int line[2] = {a,b};
        size_t id = mesh.lines().n();
        mesh.lines().add(line);
        mesh.lines().set_group(id, l_mesh.lines().group(l));
        l++; current+=(c+1);
    }
    std::cout << mesh.polygons().n();
    meshb::write(mesh, "finaltest.meshb");
    //meshb::write(mesh, "A_10M.meshb");
}
UT_TEST_CASE_END(voronoi_cut)


UT_TEST_SUITE_END(cut_test_suite)