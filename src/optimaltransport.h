#pragma once

#include "mesh.h"
#include "voronoi.h"

namespace vortex
{
    double calc_triangle_area(vec3d &p0, vec3d &p1, vec3d &p2);

    double calc_spherical_triangle_area(vec3d &p0, vec3d &p1, vec3d &p2);

    double calc_rsme_error(VoronoiDiagram &voronoi, const std::vector<double> &cell_sizes);

    template <typename T>
    double calc_energy(unsigned n, const double *x, double *de_dw, void *data0);

    double calc_energy_2D(unsigned n, const double *x, double *de_dw, void *data0);

    template <typename T>
    struct nlopt_data
    {
        VoronoiDiagram &voronoi;
        T &domain;
        Vertices &vertices;
        VoronoiDiagramOptions &options;
        const std::vector<double> &cell_sizes;
        int iter;
    };
}