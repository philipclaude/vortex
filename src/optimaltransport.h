#pragma once

#include "mesh.h"
#include "voronoi.h"
#include "math/spmat.h"

namespace vortex
{
    double calc_triangle_area(vec3d &p0, vec3d &p1, vec3d &p2);

    double calc_spherical_triangle_area(vec3d &p0, vec3d &p1, vec3d &p2);

    double calc_spherical_distance(vec3d p1, vec3d p2);

    double calc_rsme_error(VoronoiDiagram &voronoi, const std::vector<double> &cell_sizes);

    double calc_gradient_norm(std::vector<double> &de_dw);

    double geogram_CVT(vec3d &site, vec3d &p1, vec3d &p2, vec3d &p3);

    double calculate_energy(VoronoiDiagram &voronoi, std::vector<double> &ideal_size, std::vector<double> &de_dw, std::set<int> &sites_visited);

    template <typename T>
    double objective_func(const std::vector<double> &x, std::vector<double> &grad, void *data0);

    void create_edge_map(VoronoiDiagram &voronoi, std::unordered_map<std::pair<int, int>, int> &edgeSiteMap);

    void build_hessian(VoronoiDiagram &voronoi, spmat<double> &hessian, std::vector<double> &de_dw, std::vector<double> &cell_area);

    template <typename T>
    struct nlopt_data
    {
        VoronoiDiagram &voronoi;
        T &domain;
        Vertices &vertices;
        VoronoiDiagramOptions &options;
        std::vector<double> &cell_sizes;
        int iter;
        double error;
        double energy;
        bool output_converge;
        std::ofstream &outputFile;
        double &energy_time;
    };
}