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
#include <iostream>
#include <cmath>
#include <vector>

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

UT_TEST_SUITE(boundary_test_suite)

template <int dim>
std::shared_ptr<trees::KdTreeNd<coord_t, index_t>> get_nearest_neighbor(
    const coord_t* p, uint64_t np, const coord_t* q, uint64_t nq,
    std::vector<index_t>& nn, const VoronoiDiagramOptions& options,
    std::shared_ptr<trees::KdTreeNd<coord_t, index_t>> ptree = nullptr) {
  Timer timer;
  trees::KdTreeOptions kdtree_opts;
  kdtree_opts.max_dim = options.max_kdtree_axis_dim;
  if (kdtree_opts.max_dim < 0) kdtree_opts.max_dim = dim;
  using kdtree_t = trees::KdTree<dim, coord_t, index_t>;
  if (!ptree) {
    timer.start();
    ptree = std::make_shared<kdtree_t>(p, np, kdtree_opts);
    timer.stop();
    if (options.verbose)
      LOG << "kdtree created in " << timer.seconds() << " s.";
  }
  if (options.interleave_neighbors) return ptree;

  auto* tree = static_cast<kdtree_t*>(ptree.get());
  timer.start();
  std::parafor_i(0, nq,
                 [&](int tid, int k) { nn[k] = tree->nearest(&q[k * dim]); });
  timer.stop();
  if (options.verbose)
    LOG << "nearest neighbors computed in " << timer.seconds() << " s.";
  return ptree;
}

// Function to find the intersection point of two line segments
std::vector<coord_t> findIntersection(const coord_t* edge1Start,
                                    const coord_t* edge1End,
                                    const coord_t* edge2Start,
                                    const coord_t* edge2End) {
    double x1 = edge1Start[0], y1 = edge1Start[1];
    double x2 = edge1End[0], y2 = edge1End[1];
    double x3 = edge2Start[0], y3 = edge2Start[1];
    double x4 = edge2End[0], y4 = edge2End[1];

    // std::cout << "Edge 1: (" << x1 << ", " << y1 << ")" << std::endl;
    // std::cout << "Edge 1: (" << x2 << ", " << y2 << ")" << std::endl;
    // std::cout << "Edge 2: (" << x3 << ", " << y3 << ")" << std::endl;
    // std::cout << "Edge 2: (" << x4 << ", " << y4 << ")" << std::endl;

    float det = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

    const double epsilon = 1e-9;

    if (fabs(det) < epsilon) {
        // Line segments are parallel
        return {-1.0, -1.0};
    }

    // std::cout << "det: " << det << std::endl;

    float t1 = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / det;
    float t2 = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / det;

    // std::cout << "t1: " << t1 << " & t2: " << t2 << "\n";

    if (t1 >= 0.0 && t1 <= 1.0 && t2 >= 0.0 && t2 <= 1.0) {
        double intersectX = x1 + t1 * (x2 - x1);
        double intersectY = y1 + t1 * (y2 - y1);

        return {intersectX, intersectY};
    }

    return {-1.0, -1.0};
}

UT_TEST_CASE(test_2D_voronoi) {

    //const double tol = 1e-12;
    static const int dim = 3;
    size_t n_sites = 20;
    std::vector<coord_t> sites(n_sites * dim, 0.0);

    /*
    sites[0] = 0.1;
    sites[1] = 0.1;
    sites[3] = 0.6;
    sites[4] = 0.2;
    sites[6] = 0.8;
    sites[7] = 0.1;
    sites[9] = 0.2;
    sites[10] = 0.4;
    sites[12] = 0.4;
    sites[13] = 0.5;
    sites[15] = 0.7;
    sites[16] = 0.6;
    sites[18] = 0.9;
    sites[19] = 0.3;
    sites[21] = 0.1;
    sites[22] = 0.8;
    sites[24] = 0.5;
    sites[25] = 0.9;
    sites[27] = 0.8;
    sites[28] = 0.9;
    */
    sites[0] = 0.1;
    sites[1] = 0.1;
    sites[3] = 0.6;
    sites[4] = 0.2;
    sites[6] = 0.8;
    sites[7] = 0.1;
    sites[9] = 0.2;
    sites[10] = 0.4;
    sites[12] = 0.4;
    sites[13] = 0.5;
    sites[15] = 0.7;
    sites[16] = 0.6;
    sites[18] = 0.9;
    sites[19] = 0.3;
    sites[21] = 0.1;
    sites[22] = 0.8;
    sites[24] = 0.5;
    sites[25] = 0.9;
    sites[27] = 0.8;
    sites[28] = 0.9;
    sites[30] = 0.15;
    sites[31] = 0.35;
    sites[33] = 0.43;
    sites[34] = 0.29;
    sites[36] = 0.88;
    sites[37] = 0.05;
    sites[39] = 0.77;
    sites[40] = 0.33;
    sites[42] = 0.57;
    sites[43] = 0.52;
    sites[45] = 0.79;
    sites[46] = 0.65;
    sites[48] = 0.98;
    sites[49] = 0.44;
    sites[51] = 0.12;
    sites[52] = 0.67;
    sites[54] = 0.54;
    sites[55] = 0.86;
    sites[57] = 0.77;
    sites[58] = 0.96;
    
    std::vector<index_t> order(n_sites);
    sort_points_on_zcurve(sites.data(), n_sites, dim, order);
    std::vector<coord_t> ordered_sites(n_sites*dim, 0.0);

    Vertices vertices(dim);
    vertices.reserve(n_sites);
    coord_t x[dim];
    for (size_t i = 0; i < n_sites; i++) {
        for (int d = 0; d < dim; d++){
            x[d] = sites[dim * order[i] + d];
            ordered_sites[dim*i + d] = x[d];
        }
        vertices.add(x);
    }

    SquareDomain domain;
    VoronoiDiagram voronoi(dim, vertices[0], n_sites);
    VoronoiDiagramOptions options;
    options.n_neighbors = 10;
    options.allow_reattempt = false;
    options.parallel = true;
    options.store_mesh = true;
    options.verbose = true;
    voronoi.vertices().clear();
    voronoi.vertices().set_dim(3);
    voronoi.polygons().clear();
    voronoi.triangles().clear();
    voronoi.compute(domain, options);

    //create test boundary object
    size_t n_boundaries = 6;
    std::vector<coord_t> eVertices(n_boundaries * dim, 0.0);
    index_t curr_n_vertices = voronoi.vertices().n();
    std::cout << "verts: " << curr_n_vertices << std::endl;

    /*
    eVertices[0] = 0.55;
    eVertices[1] = 0.2;
    eVertices[3] = 0.6;
    eVertices[4] = 0.25;
    eVertices[6] = 0.7;
    eVertices[7] = 0.6;
    eVertices[9] = 0.5;
    eVertices[10] = 0.85;
    eVertices[12] = 0.15;
    eVertices[13] = 0.75;
    eVertices[15] = 0.2;
    eVertices[16] = 0.3;
    */

    eVertices[0] = 0.55;
    eVertices[1] = 0.2;
    eVertices[3] = 0.625;
    eVertices[4] = 0.225;
    eVertices[6] = 0.7;
    eVertices[7] = 0.375;
    eVertices[9] = 0.5;
    eVertices[10] = 0.65;
    // test point on vertex
    // eVertices[9] = 0.46476;
    // eVertices[10] = 0.682038;
    eVertices[12] = 0.4;
    eVertices[13] = 0.5;
    eVertices[15] = 0.45;
    eVertices[16] = 0.3;

    for (size_t i = 0; i < n_boundaries; ++i) {
        x[0] = eVertices[dim * i];
        x[1] = eVertices[dim * i + 1];
        x[2] = eVertices[dim * i + 2];
        voronoi.vertices().add(x);
    }

    index_t line[2] = {curr_n_vertices, curr_n_vertices+1};
    size_t id = voronoi.lines().n();
    voronoi.lines().add(line);
    voronoi.lines().set_group(id, 0);
    
    line[0] = curr_n_vertices+1;
    line[1] = curr_n_vertices+2;
    id = voronoi.lines().n();
    voronoi.lines().add(line);
    voronoi.lines().set_group(id, 0);

    line[0] = curr_n_vertices+2;
    line[1] = curr_n_vertices+3;
    id = voronoi.lines().n();
    voronoi.lines().add(line);
    voronoi.lines().set_group(id, 0);

    line[0] = curr_n_vertices+3;
    line[1] = curr_n_vertices+4;
    id = voronoi.lines().n();
    voronoi.lines().add(line);
    voronoi.lines().set_group(id, 0);

    line[0] = curr_n_vertices+4;
    line[1] = curr_n_vertices+5;
    id = voronoi.lines().n();
    voronoi.lines().add(line);
    voronoi.lines().set_group(id, 0);

    line[0] = curr_n_vertices+5;
    line[1] = curr_n_vertices;
    id = voronoi.lines().n();
    voronoi.lines().add(line);
    voronoi.lines().set_group(id, 0);

    // new mesh object
    Mesh mesh(dim);

    //KD tree
    std::vector<index_t> VNN(n_boundaries);
    std::shared_ptr<trees::KdTreeNd<coord_t, index_t>> tree{nullptr};

    // switch to sites if too much time/memory for ordered sites (or if not needed)
    get_nearest_neighbor<dim>(ordered_sites.data(), n_sites, eVertices.data(), n_boundaries, VNN, options, tree);

    // site to polygon to ensure correct polygon indices
    index_t S2P[n_sites];
    for(index_t s = 0; s < n_sites; s++){
        S2P[voronoi.polygons().group(s)] = s;
    }

    // nearest polygon with parallel-sorted polygons
    std::vector<index_t> VNP(n_boundaries);
    for(index_t vn = 0; vn < n_boundaries; vn++){
        VNP[vn] = S2P[VNN[vn]];
    }

    // adjacent vertices of boundary edges
    std::vector<std::pair<index_t, index_t>> ABV;
    for(index_t b = 0; b < n_boundaries; b++){
        index_t curr = voronoi.vertices().n() - n_boundaries + b;
        if(b == 0){
            ABV.push_back(std::make_pair(voronoi.vertices().n()-1, curr+1));
        }else if(b==n_boundaries-1){
            ABV.push_back(std::make_pair(curr-1,voronoi.vertices().n()-n_boundaries));
        }else{
            ABV.push_back(std::make_pair(curr-1,curr+1));
        }
    }

    // polygon to boundary vertices
    std::vector<std::vector<coord_t>> P2BV(voronoi.polygons().n());
    for(index_t n = 0; n < n_boundaries; n++){
        P2BV[VNP[n]].push_back(voronoi.vertices().n() - n_boundaries + n);
    }

    // polygon to boundary edges
    std::vector<std::vector<std::pair<index_t, index_t>>> P2BE(voronoi.polygons().n());
    for (index_t m = 0; m < P2BV.size(); m++) {
        if (!P2BV[m].empty()) {
            int count = 0;
            for (const auto& elem : P2BV[m]) {
                index_t curr = elem - voronoi.vertices().n() + n_boundaries;
                if(count < 1){
                    P2BE[m].push_back(std::make_pair(ABV[curr].first, elem));
                }
                P2BE[m].push_back(std::make_pair(elem, ABV[curr].second));
                count++;
            }
        } 
    }

    // Creation of new restricted polygons
    for(index_t poly = 0; poly < voronoi.polygons().n(); poly++){
        index_t vertices_count = 0;
        //std::cout << "polygon: " << poly << std::endl;
        //std::cout << "len: " << voronoi.polygons().length(poly) << std::endl;
        if(!P2BE[poly].empty()){
            //std::cout << "first b vertex: " << voronoi.vertices()[P2BE[poly].front().first][0] << ", " << voronoi.vertices()[P2BE[poly].front().first][1] << std::endl;
            //std::cout << "second b vertex: " << voronoi.vertices()[P2BE[poly].front().second][0] << ", " << voronoi.vertices()[P2BE[poly].front().second][1] << std::endl;
            index_t v1 = 0;
            while(v1 < voronoi.polygons().length(poly)){
                index_t v2 = (v1 + 1)%voronoi.polygons().length(poly);
                //std::cout << "first vertex: " << voronoi.vertices()[voronoi.polygons()[poly][v1]][0] << ", " << voronoi.vertices()[voronoi.polygons()[poly][v1]][1] << std::endl;
                //std::cout << "second vertex: " << voronoi.vertices()[voronoi.polygons()[poly][v2]][0] << ", " << voronoi.vertices()[voronoi.polygons()[poly][v2]][1] << std::endl;
                std::vector<coord_t> intersection = findIntersection(voronoi.vertices()[P2BE[poly].front().first], voronoi.vertices()[P2BE[poly].front().second], voronoi.vertices()[voronoi.polygons()[poly][v1]], voronoi.vertices()[voronoi.polygons()[poly][v2]]);
                v1++;
                if(intersection[0] >= 0.0){
                    coord_t ver[dim] = {intersection[0],intersection[1],0.0};
                    mesh.vertices().add(ver);
                    vertices_count++;
                    //std::cout << "intersection: " << intersection[0] << ", " << intersection[1] << std::endl;
                    break;
                }
            }
            index_t temp = 0;
            while(temp < voronoi.polygons().length(poly)){
                v1 = (v1)%voronoi.polygons().length(poly);
                index_t v2 = (v1 + 1)%voronoi.polygons().length(poly);
                coord_t ver[dim] = {voronoi.vertices()[voronoi.polygons()[poly][v1]][0],voronoi.vertices()[voronoi.polygons()[poly][v1]][1],0.0};
                mesh.vertices().add(ver);
                vertices_count++;
                //std::cout << "vertex: " << voronoi.vertices()[voronoi.polygons()[poly][v1]][0] << ", " << voronoi.vertices()[voronoi.polygons()[poly][v1]][1] << std::endl;
                //std::cout << "second vertex: " << voronoi.vertices()[voronoi.polygons()[poly][v2]][0] << ", " << voronoi.vertices()[voronoi.polygons()[poly][v2]][1] << std::endl;
                std::vector<coord_t> intersection = findIntersection(voronoi.vertices()[P2BE[poly].back().first], voronoi.vertices()[P2BE[poly].back().second], voronoi.vertices()[voronoi.polygons()[poly][v1]], voronoi.vertices()[voronoi.polygons()[poly][v2]]);
                v1++;
                temp++;
                if(intersection[0] >= 0.0){
                    coord_t ver[dim] = {intersection[0],intersection[1],0.0};
                    mesh.vertices().add(ver);
                    vertices_count++;
                    //std::cout << "intersection: " << intersection[0] << ", " << intersection[1] << std::endl;
                    break;
                }
            }
            for(index_t bE = 1; bE < P2BE[poly].size(); bE++){
                coord_t ver[dim] = {voronoi.vertices()[P2BE[poly][P2BE[poly].size() - bE].first][0],voronoi.vertices()[P2BE[poly][P2BE[poly].size() - bE].first][1],0.0};
                mesh.vertices().add(ver);
                vertices_count++;
                //std::cout << "bvertex: " <<  voronoi.vertices()[P2BE[poly][bE].first][0] << ", " << voronoi.vertices()[P2BE[poly][bE].first][1] << std::endl;
            }
        }else{
            for(index_t v = 0; v < voronoi.polygons().length(poly); v++){
                coord_t ver[dim] = {voronoi.vertices()[voronoi.polygons()[poly][v]][0],voronoi.vertices()[voronoi.polygons()[poly][v]][1],0.0};
                mesh.vertices().add(ver);
                vertices_count++;
                //std::cout << "vertex: " << voronoi.vertices()[voronoi.polygons()[poly][v]][0] << ", " << voronoi.vertices()[voronoi.polygons()[poly][v]][1] << std::endl;
            }
        }
        index_t new_polygon[vertices_count];
        for(index_t vc = 0; vc < vertices_count; vc++){
            new_polygon[vc] = mesh.vertices().n() - vertices_count + vc;
        }
        mesh.polygons().add(new_polygon,vertices_count);
    }

    for(int i = 0; i < mesh.polygons().n(); i++){
        std::cout << "Polygon " << i << ": ";
        for(int j = 0; j < mesh.polygons().length(i); j++){
            std::cout << "(" << mesh.vertices()[mesh.polygons()[i][j]][0] << ", " << mesh.vertices()[mesh.polygons()[i][j]][1] << ") ";
        }
        std::cout << std::endl;
    }

    meshb::write(voronoi, "2dVtest.meshb");
    meshb::write(mesh, "2dVtestUpdated.meshb");

}
UT_TEST_CASE_END(test_2D_voronoi)

/*

UT_TEST_CASE(test_2D_voronoi) {

    //const double tol = 1e-12;
    static const int dim = 4;
    size_t n_sites = 10;
    std::vector<coord_t> sites(n_sites * dim, 0.0);

    sites[0] = 0.1;
    sites[1] = 0.1;
    sites[4] = 0.6;
    sites[5] = 0.2;
    sites[8] = 0.8;
    sites[9] = 0.1;
    sites[12] = 0.2;
    sites[13] = 0.4;
    sites[16] = 0.4;
    sites[17] = 0.5;
    sites[20] = 0.7;
    sites[21] = 0.6;
    sites[24] = 0.9;
    sites[25] = 0.3;
    sites[28] = 0.1;
    sites[29] = 0.8;
    sites[32] = 0.5;
    sites[33] = 0.9;
    sites[36] = 0.8;
    sites[37] = 0.9;

    for (size_t k = 0; k < n_sites; k++) {
        sites[k * dim + 2] = 0.0;
        if (dim > 3) sites[k * dim + 3] = 0.0;
    }

    std::vector<index_t> order(n_sites);
    sort_points_on_zcurve(sites.data(), n_sites, dim, order);
    std::vector<coord_t> ordered_sites(n_sites*dim, 0.0);

    Vertices vertices(dim);
    vertices.reserve(n_sites);
    coord_t x[dim];
    for (size_t i = 0; i < n_sites; i++) {
        for (int d = 0; d < dim; d++){
            x[d] = sites[dim * order[i] + d];
            ordered_sites[dim*i + d] = x[d];
        }
        vertices.add(x);
        //std::cout << "Site " << i << ": " << x[0] << ", " << x[1] << std::endl;
    }

    SquareDomain domain;
    VoronoiDiagram voronoi(dim, vertices[0], n_sites);
    VoronoiDiagramOptions options;
    options.n_neighbors = 10;
    options.allow_reattempt = false;
    options.parallel = true;
    options.store_mesh = true;
    options.verbose = true;
    voronoi.vertices().clear();
    voronoi.vertices().set_dim(3);
    voronoi.polygons().clear();
    voronoi.triangles().clear();
    voronoi.compute(domain, options);

    //create test boundary object
    size_t n_eVertices = 6;
    std::vector<coord_t> eVertices(n_eVertices * dim, 0.0);

    eVertices[0] = 0.55;
    eVertices[1] = 0.2;
    eVertices[4] = 0.6;
    eVertices[5] = 0.25;
    eVertices[8] = 0.7;
    eVertices[9] = 0.6;
    eVertices[12] = 0.5;
    eVertices[13] = 0.85;
    eVertices[16] = 0.15;
    eVertices[17] = 0.75;
    eVertices[20] = 0.2;
    eVertices[21] = 0.3;

    for (size_t i = 0; i < n_eVertices; ++i) {
        x[0] = eVertices[dim * i];
        x[1] = eVertices[dim * i + 1];
        x[2] = eVertices[dim * i + 2];
        if(dim > 3){
            x[3] = eVertices[dim * i + 3];
        }
        voronoi.vertices().add(x);
    }

    int line[2] = {59, 60};
    size_t id = voronoi.lines().n();
    voronoi.lines().add(line);
    voronoi.lines().set_group(id, 0);
    
    line[0] = 60;
    line[1] = 61;
    id = voronoi.lines().n();
    voronoi.lines().add(line);
    voronoi.lines().set_group(id, 0);

    line[0] = 61;
    line[1] = 62;
    id = voronoi.lines().n();
    voronoi.lines().add(line);
    voronoi.lines().set_group(id, 0);

    line[0] = 62;
    line[1] = 63;
    id = voronoi.lines().n();
    voronoi.lines().add(line);
    voronoi.lines().set_group(id, 0);

    line[0] = 63;
    line[1] = 64;
    id = voronoi.lines().n();
    voronoi.lines().add(line);
    voronoi.lines().set_group(id, 0);

    line[0] = 64;
    line[1] = 59;
    id = voronoi.lines().n();
    voronoi.lines().add(line);
    voronoi.lines().set_group(id, 0);

    //KD tree
    std::vector<index_t> vnn(n_eVertices);
    std::shared_ptr<trees::KdTreeNd<coord_t, index_t>> tree{nullptr};

    // switch to sites if too much time/memory for ordered sites (or if not needed)
    get_nearest_neighbor<4>(ordered_sites.data(), n_sites, eVertices.data(), n_eVertices, vnn, options, tree);

    for(int c = 0; c < n_eVertices; ++c){
    //for(int c = 1; c < 2; ++c){
        int endpoint1 = vnn[c];
        int endpoint2 = vnn[(c + 1)%n_eVertices];

        //std::cout << "end1: " << endpoint1 << std::endl;
        //std::cout << "end2: " << endpoint2 << std::endl;

        if(endpoint1==endpoint2){
            //in same polygon
        }else{
            // more efficient way to find polygons?
            int matched = 0;
            int polygon1 = -1;
            int polygon2 = -1;
            for(int p = 0; p < voronoi.polygons().n(); p++){
                if(voronoi.polygons().group(p) == endpoint1){
                    polygon1 = p;
                    matched++;
                }else if(voronoi.polygons().group(p) == endpoint2){
                    polygon2 = p;
                    matched++;
                }
                if(matched > 1){
                    break;
                }
            }
            coord_t e2S[dim];
            coord_t e2E[dim];
            matched = 0;
            const double epsilon = 1e-9;
            for(int v = 0; v < voronoi.polygons().length(polygon1); v++){
                for(int w = 0; w < voronoi.polygons().length(polygon2); w++){
                    if(abs(voronoi.vertices()[voronoi.polygons()[polygon1][v]][0]-voronoi.vertices()[voronoi.polygons()[polygon2][w]][0]) < epsilon
                    && abs(voronoi.vertices()[voronoi.polygons()[polygon1][v]][1]-voronoi.vertices()[voronoi.polygons()[polygon2][w]][1]) < epsilon){
                        if(matched == 0){
                            e2S[0] = voronoi.vertices()[voronoi.polygons()[polygon1][v]][0];
                            e2S[1] = voronoi.vertices()[voronoi.polygons()[polygon1][v]][1];
                            e2S[2] = e2S[3] = 0.0;
                        }else{
                            e2E[0] = voronoi.vertices()[voronoi.polygons()[polygon1][v]][0];
                            e2E[1] = voronoi.vertices()[voronoi.polygons()[polygon1][v]][1];
                            e2E[2] = e2E[3] = 0.0;
                        }
                        matched++;
                        break;
                    }
                }
                if(matched > 1){
                    break;
                }
            }
            int a = voronoi.vertices().n() - n_eVertices + c;
            int b = voronoi.vertices().n() - n_eVertices + ((c + 1)%n_eVertices);
            std::vector<coord_t> intersection = findIntersection(e2S, e2E, voronoi.vertices()[a], voronoi.vertices()[b]);
            std::cout << "(" << intersection[0] << ", " << intersection[1] << ")" << std::endl;
        }
    }
    //meshb::write(voronoi, "2dVtest.meshb");

}
UT_TEST_CASE_END(test_2D_voronoi)

// Function to determine if a point is on the left side of a line segment
bool isOutside(const coord_t* edgeStart, const coord_t* edgeEnd, const coord_t* point) {
    // init values
    double p1x = edgeStart[0];
    double p1y = edgeStart[1];
    double p2x = edgeEnd[0];
    double p2y = edgeEnd[1];
    double qx = point[0];
    double qy = point[1];

    // Calculate vectors from p1 to p2 and from p1 to the point
    double vec1_x = p2x - p1x;
    double vec1_y = p2y - p1y;
    double vec2_x = qx - p1x;
    double vec2_y = qy - p1y;

    // Calculate the cross product
    double crossProduct = vec1_x * vec2_y - vec1_y * vec2_x;

    // Check if the cross product is positive (point is on the left side) or negative (point is on the right side)
    return crossProduct < 0;
}

// Function to check if a point is on the line segment defined by two points
bool isPointOnSegment(const Point& p1, const Point& p2, const Point& q) {
    // Check if the point q lies on the line segment p1p2
    float crossProduct = (q.y - p1.y) * (p2.x - p1.x) - (q.x - p1.x) * (p2.y - p1.y);
    if (std::abs(crossProduct) > 1e-6) // Adjust epsilon as needed for your precision requirements
        return false; // Point q is not collinear with p1 and p2

    // Check if the point q lies within the bounding box formed by the endpoints p1 and p2
    return q.x >= std::min(p1.x, p2.x) && q.x <= std::max(p1.x, p2.x) &&
           q.y >= std::min(p1.y, p2.y) && q.y <= std::max(p1.y, p2.y);
}

// Function to determine if a point is inside, on the boundary, or outside the polygon
int pointInPolygon(const std::vector<std::pair<int, int>>& polygon, const std::vector<Point>& vertices, const Point& q) {
    int n = polygon.size();
    if (n < 3) return -1; // Invalid polygon

    int crossings = 0;
    for (int i = 0; i < n; ++i) {
        int v1_index = polygon[i].first;
        int v2_index = polygon[i].second;
        Point p1 = vertices[v1_index];
        Point p2 = vertices[v2_index];

        //std::cout << "Edge: (" << p1.x << "," << p1.y << ") & (" << p2.x << "," << p2.y << ")" << std::endl;

        // Check if the point is on the line segment
        if (isPointOnSegment(p1, p2, q))
            return 1; // Point is on the boundary

        // Check for edge intersections
        if (((p1.y > q.y) != (p2.y > q.y)) &&
            (q.x < (p2.x - p1.x) * (q.y - p1.y) / (p2.y - p1.y) + p1.x))
            crossings++;
    }

    return (crossings % 2 == 1) ? 0 : 2; // Inside or outside
}
*/

/*
UT_TEST_CASE(test_grid_line) {
    int n = 10;
    Grid<Polygon> mesh({n, n});
    mesh.fields().set_defaults(mesh);

    LOG << fmt::format("writing {} vertices", mesh.vertices().n());

    //setting boundary lines
    int line[2] = {0, 13};
    size_t id = mesh.lines().n();
    mesh.lines().add(line);
    mesh.lines().set_group(id, 0);

    line[0] = 13;
    line[1] = 37;
    id = mesh.lines().n();
    mesh.lines().add(line);
    mesh.lines().set_group(id, 0);

    line[0] = 37;
    line[1] = 80;
    id = mesh.lines().n();
    mesh.lines().add(line);
    mesh.lines().set_group(id, 0);

    line[0] = 80;
    line[1] = 115;
    id = mesh.lines().n();
    mesh.lines().add(line);
    mesh.lines().set_group(id, 0);

    line[0] = 115;
    line[1] = 110;
    id = mesh.lines().n();
    mesh.lines().add(line);
    mesh.lines().set_group(id, 0);

    line[0] = 110;
    line[1] = 0;
    id = mesh.lines().n();
    mesh.lines().add(line);
    mesh.lines().set_group(id, 0);

    LOG << fmt::format("writing {} lines", mesh.lines().n());
    LOG << fmt::format("writing {} polygons\n", mesh.polygons().n());

    std::vector<std::vector<std::pair<int, int>>> continents;

    int l = 0;
    while(l < mesh.lines().n()){
        std::vector<std::pair<int, int>> continent;
        int start = mesh.lines()[l][0];
        int end = mesh.lines()[l][1];
        do {
            int p1 = mesh.lines()[l][0];
            int q1 = mesh.lines()[l][1];
            continent.push_back({p1, q1});
            end = q1;
            l++;
        }while(start != end);
        continents.push_back(continent);
    }

    
    // Print continents
    for (const auto& inner_vec : continents) {
        for (const auto& pair : inner_vec) {
            std::cout << "(" << pair.first << ", " << pair.second << ") ";
        }
        std::cout << std::endl;
    }
    
    std::vector<Point> temp_vertices;
    for(int i = 0; i < mesh.vertices().n(); i++){
        float tX = mesh.vertices()[i][0];
        float tY = mesh.vertices()[i][1];
        Point temp = {tX,tY};
        temp_vertices.push_back(temp);
    }

    Point point = {0.5, 0.9};

    for(int c = 0; c < continents.size(); c++){

        // Determine if the point is inside, on the boundary, or outside the polygon
        int result = pointInPolygon(continents[c], temp_vertices, point);

        // Output the result
        if (result == 0) {
            std::cout << "Point is inside the polygon" << std::endl;
        } else if (result == 1) {
            std::cout << "Point is on the boundary of the polygon" << std::endl;
        } else if (result == 2) {
            std::cout << "Point is outside the polygon" << std::endl;
        } else {
            std::cout << "Invalid polygon" << std::endl;
        }
    }

    //meshb::write(mesh, "gridtest.meshb");
}
UT_TEST_CASE_END(test_grid_line)

UT_TEST_CASE(test_grid_line) {
    int n = 10;
    Grid<Polygon> mesh({n, n});
    mesh.fields().set_defaults(mesh);

    LOG << fmt::format("writing {} vertices", mesh.vertices().n());

    //setting boundary lines
    int line[2] = {0, 13};
    size_t id = mesh.lines().n();
    mesh.lines().add(line);
    mesh.lines().set_group(id, 0);

    line[0] = 13;
    line[1] = 37;
    id = mesh.lines().n();
    mesh.lines().add(line);
    mesh.lines().set_group(id, 0);

    line[0] = 37;
    line[1] = 80;
    id = mesh.lines().n();
    mesh.lines().add(line);
    mesh.lines().set_group(id, 0);

    line[0] = 80;
    line[1] = 115;
    id = mesh.lines().n();
    mesh.lines().add(line);
    mesh.lines().set_group(id, 0);

    line[0] = 115;
    line[1] = 110;
    id = mesh.lines().n();
    mesh.lines().add(line);
    mesh.lines().set_group(id, 0);

    line[0] = 110;
    line[1] = 0;
    id = mesh.lines().n();
    mesh.lines().add(line);
    mesh.lines().set_group(id, 0);

    LOG << fmt::format("writing {} lines", mesh.lines().n());
    LOG << fmt::format("writing {} polygons\n", mesh.polygons().n());

    std::vector<std::pair<float, float>> new_vertices;
    std::vector<std::string> new_polygons;

    //for(int i = 0; i < 400; i += 4){
    for(int i = 380; i < 384; i += 4){
        //loop through polygon edges
        for(int j = 0; j < 4; j++){
            //endpoints of current polygon edge
            int p2 = mesh.polygons().data()[i+j];
            int q2 = mesh.polygons().data()[i + (j + 1) % 4];
            //coordinates of endpoints (p2 and q2)
            float p2X = mesh.vertices()[p2][0];
            float p2Y = mesh.vertices()[p2][1];
            float q2X = mesh.vertices()[q2][0];
            float q2Y = mesh.vertices()[q2][1];
            LOG << "edge: (" << p2X << "," << p2Y << ") & (" << q2X << "," << q2Y << ")";

            std::vector<float> edge2Start = {p2X, p2Y};
            std::vector<float> edge2End = {q2X, q2Y};

            //loop through boundary lines
            for(int l = 0; l < mesh.lines().n(); l++){
            //for(int l = 5; l < 6; l++){
                //endpoints of current boundary edge
                int p1 = mesh.lines()[l][0];
                int q1 = mesh.lines()[l][1];
                //coordinates of endpoints (p1 and q1)
                float p1X = mesh.vertices()[p1][0];
                float p1Y = mesh.vertices()[p1][1];
                float q1X = mesh.vertices()[q1][0];
                float q1Y = mesh.vertices()[q1][1];
                LOG << "boundary: (" << p1X << "," << p1Y << ") & (" << q1X << "," << q1Y << ")";

                //intersecting boundary edge
                std::vector<float> edge1Start = {p1X, p1Y};
                std::vector<float> edge1End = {q1X, q1Y};

                std::vector<float> intersection = findIntersection(edge1Start, edge1End, edge2Start, edge2End);

                if (intersection[0] != -1.0 && intersection[1] != -1.0) {
                    new_vertices.push_back(std::make_pair(intersection[0], intersection[1]));
                    LOG << "Intersection Point: (" << intersection[0] << ", " << intersection[1] << ")";
                } else {
                    LOG << "Line segments do not intersect.";
                }
                bool isOnOutside = isOutside(p1X, p1Y, q1X, q1Y, q2X, q2Y);
                if(isOnOutside){
                    float lowerX = (p1X < q1X) ? p1X : q1X;
                    float higherX = (p1X < q1X) ? q1X : p1X;
                    float lowerY = (p1Y < q1Y) ? p1Y : q1Y;
                    float higherY = (p1Y < q1Y) ? q1Y : p1Y;

                    if(q2X >= lowerX && q2X <= higherX && q2Y >= lowerY && q2Y <= higherY){
                        new_vertices.push_back(std::make_pair(q2X, q2Y));
                    }
                    LOG << "Endpoint is outside.\n";
                }else{
                    LOG << "Endpoint is inside.\n";
                }
            }    
        }
    }

    
    for (const auto& pair : new_vertices) {
        std::cout << "(" << pair.first << ", " << pair.second << ")" << std::endl;
    }
    

    //meshb::write(mesh, "gridtest.meshb");
}
UT_TEST_CASE_END(test_grid_line)
*/

UT_TEST_SUITE_END(boundary_test_suite)