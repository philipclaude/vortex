#include "optimaltransport.h"
#include "math/mat.hpp"
#include "math/vec.hpp"
#include "util.h"

namespace vortex
{

    double calc_triangle_area(vec3d &p0, vec3d &p1, vec3d &p2)
    {
        vec3d edge1 = p1 - p0;
        vec3d edge2 = p2 - p0;

        vec3d crossProduct = cross(edge1, edge2);
        return 0.5 * length(crossProduct);
    };

    double calc_spherical_triangle_area(vec3d &p0, vec3d &p1, vec3d &p2)
    {
        coord_t num = std::fabs(dot(p0, cross(p1, p2)));

        coord_t den = 1.0 + dot(p0, p1) + dot(p1, p2) + dot(p0, p2);
        coord_t ak = 2.0 * std::atan2(num, den);
        return ak;
    };

    double calc_rsme_error(VoronoiDiagram &voronoi, const std::vector<double> &cell_sizes)
    {
        const Topology<Polygon> &polygons = voronoi.polygons();
        double error = 0.0;
        double sum = 0.0;
        for (int k = 0; k < polygons.n(); k++)
        {
            int group_index = polygons.group(k);
            auto props = voronoi.analyze();
            sum += pow((voronoi.properties()[group_index].mass - cell_sizes[group_index]), 2);
            LOG << fmt::format("cell_sizes[{}] = {}, mass = {}", group_index, cell_sizes[group_index], voronoi.properties()[group_index].mass);
        }
        error = pow((sum / polygons.n()), 0.5);

        return error;
    }

    template <typename T>
    double calc_energy(unsigned n, const double *x, double *de_dw, void *data0)
    {
        nlopt_data<T> &data = *static_cast<nlopt_data<T> *>(data0);
        VoronoiDiagram &voronoi = data.voronoi;
        const Topology<Polygon> &polygons = voronoi.polygons();
        auto &weights = voronoi.weights();
        weights.resize(n, 0.0);
        for (int i = 0; i < n; i++)
        {
            weights[i] = x[i];
        }
        lift_sites(data.vertices, voronoi.weights());

        data.options.store_mesh = true;
        data.domain.set_initialization_fraction(0.7);
        voronoi.vertices().clear();
        voronoi.polygons().clear();
        voronoi.triangles().clear();
        voronoi.compute(data.domain, data.options);

        data.iter++;

        double energy = 0.0;
        for (int k = 0; k < polygons.n(); k++)
        {
            double cell_area = 0.0;
            // get the site of the polygon
            int group_index = polygons.group(k);
            vec3d site(voronoi.vertices()[group_index]);
            double curr_weight = weights[group_index];

            int p1_index = polygons(k, 0);
            vec3d p1(voronoi.vertices()[p1_index]);

            for (int j = 1; j < voronoi.polygons().length(k) - 1; j++)
            {
                // get two vertices of the polygon
                int p2_index = polygons(k, j);
                int p3_index = polygons(k, j + 1);
                vec3d p2(voronoi.vertices()[p2_index]);
                vec3d p3(voronoi.vertices()[p3_index]);

                // calculate area of the triangle
                double t_area = calc_spherical_triangle_area(p1, p2, p3);

                // taken directly from Levy Geogram
                double cur_f = 0.0;
                for (index_t c = 0; c < 2; c++)
                {
                    double u0 = site[c] - p1[c];
                    double u1 = site[c] - p2[c];
                    double u2 = site[c] - p3[c];
                    cur_f += u0 * u0;
                    cur_f += u1 * (u0 + u1);
                    cur_f += u2 * (u0 + u1 + u2);
                }

                energy += t_area * cur_f / 6.0;
                cell_area += t_area;
            }

            if (de_dw)
            {
                de_dw[group_index] = cell_area - data.cell_sizes[group_index];
            }
            energy += (-curr_weight * cell_area);
            energy += data.cell_sizes[k] * curr_weight;
        }

        return -energy;
    };
}
