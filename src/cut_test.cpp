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
vec3d findIntersection(const coord_t* edge1Start, const coord_t* edge1End, const coord_t* edge2Start, const coord_t* edge2End) {
    vec3d e1S = edge1Start; vec3d e1E = edge1End;vec3d e2S = edge2Start; vec3d e2E = edge2End;
    // Create plane of intersecting edge
    vec3d i = e1E - e1S;
    vec3d mid = (e1S + e1E)/2.0;
    vec3d ns = normalize(mid);
    vec3d nplane = normalize(cross(i,ns));

    vec3d e = e2E - e2S;
    double t = dot(nplane,(mid-e2S)) / dot(nplane,e);

    vec3d r = e2S + t*e;

    return normalize(r);
}

UT_TEST_CASE(voronoi_cut) {
    // read in mesh to be edited
    static const int dim = 3;
    Mesh o_mesh(dim);
    meshb::read("../build/release/example7.meshb", o_mesh);
    size_t n_sites = o_mesh.polygons().n();

    // create boundary object (l_mesh) and place lines on original mesh (o_mesh)
    Mesh l_mesh(dim);
    meshb::read("../build/release/water.meshb", l_mesh);
    index_t current = o_mesh.vertices().n();
    index_t l = 0;
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

    // KD tree
    size_t n_boundaries = o_mesh.lines().n();
    std::vector<index_t> VNN(o_mesh.lines().n());
    std::shared_ptr<trees::KdTreeNd<coord_t, index_t>> tree{nullptr};
    get_nearest_neighbor<dim>(o_mesh.vertices()[0], n_sites, o_mesh.vertices()[o_mesh.vertices().n()-n_boundaries], n_boundaries, VNN, tree);

    // site to polygon to ensure correct polygon indices
    index_t S2P[n_sites];
    for(index_t s = 0; s < n_sites; s++){
        S2P[o_mesh.polygons().group(s)] = s;
    }

    // nearest polygon with parallel-sorted polygons
    std::vector<index_t> VNP(n_boundaries);
    for(index_t vn = 0; vn < n_boundaries; vn++){
        VNP[vn] = S2P[VNN[vn]];
    }

    // adjacent vertices of boundary edges
    std::vector<std::pair<index_t, index_t>> ABV;
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

    // polygon to boundary vertices
    std::vector<std::vector<coord_t>> P2BV(o_mesh.polygons().n());
    for(index_t n = 0; n < n_boundaries; n++){
        P2BV[VNP[n]].push_back(o_mesh.vertices().n() - n_boundaries + n);
    }

    // polygon to boundary edges
    std::vector<std::vector<std::pair<index_t, index_t>>> P2BE(o_mesh.polygons().n());
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

    // Creation of new restricted polygons
    for(index_t poly = 0; poly < o_mesh.polygons().n(); poly++){
        index_t vertices_count = 0;
        if(!P2BE[poly].empty()){
            index_t v1 = 0;
            while(v1 < o_mesh.polygons().length(poly)){
                index_t v2 = (v1 + 1)%o_mesh.polygons().length(poly);
                if(planeSide(o_mesh.vertices()[o_mesh.polygons()[poly][v1]],o_mesh.vertices()[o_mesh.polygons()[poly][v2]],o_mesh.vertices()[P2BE[poly].front().first])!=planeSide(o_mesh.vertices()[o_mesh.polygons()[poly][v1]],o_mesh.vertices()[o_mesh.polygons()[poly][v2]],o_mesh.vertices()[P2BE[poly].front().second])){
                //if(planeSide(o_mesh.vertices()[P2BE[poly].front().first],o_mesh.vertices()[P2BE[poly].front().second],o_mesh.vertices()[o_mesh.polygons()[poly][v1]])!=planeSide(o_mesh.vertices()[P2BE[poly].front().first],o_mesh.vertices()[P2BE[poly].front().second],o_mesh.vertices()[o_mesh.polygons()[poly][v2]])){
                    vec3d intersection = findIntersection(o_mesh.vertices()[P2BE[poly].front().first], o_mesh.vertices()[P2BE[poly].front().second], o_mesh.vertices()[o_mesh.polygons()[poly][v1]], o_mesh.vertices()[o_mesh.polygons()[poly][v2]]);
                    v1++;
                    coord_t ver[dim] = {intersection[0],intersection[1],intersection[2]};
                    mesh.vertices().add(ver);
                    vertices_count++;
                    break;
                }
                v1++;
            }
            index_t temp = 0;
            while(temp < o_mesh.polygons().length(poly)){
                v1 = (v1)%o_mesh.polygons().length(poly);
                index_t v2 = (v1 + 1)%o_mesh.polygons().length(poly);
                coord_t ver[dim] = {o_mesh.vertices()[o_mesh.polygons()[poly][v1]][0],o_mesh.vertices()[o_mesh.polygons()[poly][v1]][1],o_mesh.vertices()[o_mesh.polygons()[poly][v1]][2]};
                mesh.vertices().add(ver);
                vertices_count++;
                if(planeSide(o_mesh.vertices()[o_mesh.polygons()[poly][v1]],o_mesh.vertices()[o_mesh.polygons()[poly][v2]],o_mesh.vertices()[P2BE[poly].back().first])!=planeSide(o_mesh.vertices()[o_mesh.polygons()[poly][v1]],o_mesh.vertices()[o_mesh.polygons()[poly][v2]],o_mesh.vertices()[P2BE[poly].back().second])){
                //if(planeSide(o_mesh.vertices()[P2BE[poly].back().first],o_mesh.vertices()[P2BE[poly].back().second],o_mesh.vertices()[o_mesh.polygons()[poly][v1]])!=planeSide(o_mesh.vertices()[P2BE[poly].back().first],o_mesh.vertices()[P2BE[poly].back().second],o_mesh.vertices()[o_mesh.polygons()[poly][v2]])){
                    vec3d intersection = findIntersection(o_mesh.vertices()[P2BE[poly].back().first], o_mesh.vertices()[P2BE[poly].back().second], o_mesh.vertices()[o_mesh.polygons()[poly][v1]], o_mesh.vertices()[o_mesh.polygons()[poly][v2]]);
                    v1++;
                    temp++;
                    coord_t ver[dim] = {intersection[0],intersection[1],intersection[2]};
                    mesh.vertices().add(ver);
                    vertices_count++;
                    break;
                }
                v1++;
                temp++;
            }
            for(index_t bE = 1; bE < P2BE[poly].size(); bE++){
                coord_t ver[dim] = {o_mesh.vertices()[P2BE[poly][P2BE[poly].size() - bE].first][0],o_mesh.vertices()[P2BE[poly][P2BE[poly].size() - bE].first][1],o_mesh.vertices()[P2BE[poly][P2BE[poly].size() - bE].first][2]};
                mesh.vertices().add(ver);
                vertices_count++;
            }
        }else{
            for(index_t v = 0; v < o_mesh.polygons().length(poly); v++){
                coord_t ver[dim] = {o_mesh.vertices()[o_mesh.polygons()[poly][v]][0],o_mesh.vertices()[o_mesh.polygons()[poly][v]][1],o_mesh.vertices()[o_mesh.polygons()[poly][v]][2]};
                mesh.vertices().add(ver);
                vertices_count++;
            }
        }
        index_t new_polygon[vertices_count];
        for(index_t vc = 0; vc < vertices_count; vc++){
            new_polygon[vc] = mesh.vertices().n() - vertices_count + vc;
        }
        mesh.polygons().add(new_polygon,vertices_count);
    }

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
    //meshb::write(mesh, "finaltest.meshb");

}
UT_TEST_CASE_END(voronoi_cut)


UT_TEST_SUITE_END(cut_test_suite)