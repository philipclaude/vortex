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

UT_TEST_SUITE(adjacent_test_suite)

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

UT_TEST_CASE(voronoi_cut_edge) {
    // read in mesh to be edited
    static const int dim = 3;
    Mesh o_mesh(dim);
    //meshb::read("../build/release/example7.meshb", o_mesh);
    meshb::read("../build/release/analytic_10M.meshb", o_mesh);
    size_t n_sites = o_mesh.polygons().n();

    // create boundary object (l_mesh) and place lines on original mesh (o_mesh)
    Mesh l_mesh(dim);
    meshb::read("../build/release/water_0025.meshb", l_mesh);
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
    mesh.vertices().reserve(o_mesh.vertices().n());
    mesh.polygons().reserve(n_sites);

    // texture options
    std::string tex_file =
        std::string(VORTEX_SOURCE_DIR) + "/../data/oceans_2048.png";
    TextureOptions tex_opts;
    tex_opts.format = TextureFormat::kGrayscale;
    Texture texture(tex_file, tex_opts);
    texture.make_binary(10, 10, 255);

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

    // adjacent vertices of boundary edges (use vector reserve [upper bound] and shrink to fit after loop)
    std::vector<std::pair<index_t, index_t>> ABV(n_boundaries);
    // polygon to boundary vertices (possible to reserve inner vector?)
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
    
    Timer timer;
    timer.start();
    // Create map of adjacent polygons
    std::unordered_map<std::pair<index_t, index_t>, std::pair<index_t, index_t>> adjMap;
    adjMap.reserve(o_mesh.polygons().n()*10);
    for(index_t poly = 0; poly < o_mesh.polygons().n(); poly++){
        for(index_t e = 0; e < o_mesh.polygons().length(poly); e++){
            index_t e0 = o_mesh.polygons()[poly][e];
            index_t e1 = o_mesh.polygons()[poly][(e+1)%o_mesh.polygons().length(poly)];
            std::pair<index_t, index_t> check_key = std::make_pair(e1, e0);
            auto it = adjMap.find(check_key);
            if (it != adjMap.end()) {
                adjMap[check_key].second = poly;
            } else {
                adjMap[std::make_pair(e0, e1)] = std::make_pair(poly, -1);
            }
        }
    }
    timer.stop();
    LOG << fmt::format("Adjacencies found in {} seconds", timer.seconds());

    // map of old polygon indices to new indices
    std::unordered_map<index_t, index_t> indices;
    indices.reserve(o_mesh.polygons().n());

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
        if(water){   
            index_t vertices_count = 0;
            if(!P2BE[poly].empty()){
                index_t v1 = 0;
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
            indices[poly] = mesh.polygons().n();
            mesh.polygons().add(new_polygon,vertices_count);
        }
    }

    timer.start();
    // adjacent test
    for(index_t poly = 0; poly < o_mesh.polygons().n(); poly++){
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
        if(water){  
            if(!P2BE[poly].empty()){
                // std::cout << "first polygon: " << poly << std::endl;
                // std::cout << "P2BE: " << P2BE[poly].back().second << std::endl;
                std::pair<index_t, index_t> edge;
                coord_t firstIntersection[dim];
                for(int e = 0; e < o_mesh.polygons().length(poly); e++){
                    index_t v1 = e; index_t v2 = (e + 1) % o_mesh.polygons().length(poly);
                    if(planeSide(o_mesh.vertices()[o_mesh.polygons()[poly][v1]],o_mesh.vertices()[o_mesh.polygons()[poly][v2]],o_mesh.vertices()[P2BE[poly].back().first])!=planeSide(o_mesh.vertices()[o_mesh.polygons()[poly][v1]],o_mesh.vertices()[o_mesh.polygons()[poly][v2]],o_mesh.vertices()[P2BE[poly].back().second])){
                        vec3d intersection;
                        if(findIntersection(o_mesh.vertices()[P2BE[poly].back().first], o_mesh.vertices()[P2BE[poly].back().second], o_mesh.vertices()[o_mesh.polygons()[poly][v1]], o_mesh.vertices()[o_mesh.polygons()[poly][v2]],intersection)){
                            edge.first = o_mesh.polygons()[poly][v1]; edge.second = o_mesh.polygons()[poly][v2];
                            firstIntersection[0] = intersection[0]; firstIntersection[1] = intersection[1]; firstIntersection[2] = intersection[2];
                            break;
                        }
                    }
                }
                index_t first_endpoint = P2BE[poly].back().first;
                index_t second_endpoint = P2BE[poly].back().second;
                //find the next polygon
                index_t next_poly;
                if (adjMap.find(edge) != adjMap.end()) {
                    next_poly = adjMap[edge].second;
                } else {
                    std::pair<index_t, index_t> key = {edge.second, edge.first};
                    next_poly = adjMap[key].first;
                }
                // std::cout << "next polygon: " << next_poly << std::endl;
                // std::cout << "P2BE: " << P2BE[next_poly].front().second << std::endl;
                // while(P2BE[next_poly].front().second != second_endpoint){
                while(P2BE[next_poly].empty()){
                    mesh.vertices().add(firstIntersection);
                    size_t vert_count = 1;
                    size_t n = o_mesh.polygons().length(next_poly);
                    index_t eIndex;
                    // find index of the edge
                    for(index_t e = 0; e < n; e++){
                        if(o_mesh.polygons()[next_poly][e] == edge.second && o_mesh.polygons()[next_poly][(e+1)%n] == edge.first){
                            eIndex = (e+1)%n;
                            break;
                        }
                    }
                    // loop through polygon edge until second intersection point of edge is found, adding vertices along the way
                    for(index_t e = 0; e < n; e++){
                        index_t nv1 = (eIndex + e)%n;
                        index_t nv2 = (nv1 + 1)%n;
                        coord_t ver[dim] = {o_mesh.vertices()[o_mesh.polygons()[next_poly][nv1]][0],o_mesh.vertices()[o_mesh.polygons()[next_poly][nv1]][1],o_mesh.vertices()[o_mesh.polygons()[next_poly][nv1]][2]};
                        mesh.vertices().add(ver);
                        vert_count++;
                        if(planeSide(o_mesh.vertices()[o_mesh.polygons()[next_poly][nv1]],o_mesh.vertices()[o_mesh.polygons()[next_poly][nv2]],o_mesh.vertices()[first_endpoint])!=planeSide(o_mesh.vertices()[o_mesh.polygons()[next_poly][nv1]],o_mesh.vertices()[o_mesh.polygons()[next_poly][nv2]],o_mesh.vertices()[second_endpoint])){
                            vec3d intersection;
                            if(findIntersection(o_mesh.vertices()[first_endpoint], o_mesh.vertices()[second_endpoint], o_mesh.vertices()[o_mesh.polygons()[next_poly][nv1]], o_mesh.vertices()[o_mesh.polygons()[next_poly][nv2]],intersection)){
                                firstIntersection[0] = intersection[0]; firstIntersection[1] = intersection[1]; firstIntersection[2] = intersection[2];
                                mesh.vertices().add(firstIntersection);
                                vert_count++;
                                edge.first = o_mesh.polygons()[next_poly][nv1]; edge.second = o_mesh.polygons()[next_poly][nv2];
                                break;
                            }
                        }
                    }
                    // create polygon
                    index_t new_polygon[vert_count];
                    for(index_t vc = 0; vc < vert_count; vc++){
                        new_polygon[vc] = mesh.vertices().n() - vert_count + vc;
                    }
                    mesh.polygons().add(new_polygon,vert_count);
                    // remove old version
                    index_t index = indices[next_poly];
                    for(int v = 1; v < mesh.polygons().length(index); v++){
                        mesh.vertices()[mesh.polygons()[index][v]][0] = mesh.vertices()[mesh.polygons()[index][0]][0];
                        mesh.vertices()[mesh.polygons()[index][v]][1] = mesh.vertices()[mesh.polygons()[index][0]][1];
                        mesh.vertices()[mesh.polygons()[index][v]][2] = mesh.vertices()[mesh.polygons()[index][0]][2];
                    }
                    // set the next polygon
                    if (adjMap.find(edge) != adjMap.end()) {
                        next_poly = adjMap[edge].second;
                    } else {
                        std::pair<index_t, index_t> key = {edge.second, edge.first};
                        next_poly = adjMap[key].first;
                    }
                    // std::cout << "last polygon: " << next_poly << std::endl;
                    // std::cout << "P2BE: " << P2BE[next_poly].front().second << std::endl;
                }
            }
        }
    }

    timer.stop();
    LOG << fmt::format("Adjacency test completed in {} seconds", timer.seconds());

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

    //meshb::write(mesh, "edgetest.meshb");
}
UT_TEST_CASE_END(voronoi_cut_edge)


UT_TEST_SUITE_END(adjacent_test_suite)