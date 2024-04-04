#include "optimaltransport.h"
#include "math/mat.hpp"
#include "math/vec.hpp"
#include "quadrature.hpp"
#include "util.h"

namespace vortex
{

    double calc_triangle_area(vec3d p0, vec3d p1, vec3d p2)
    {
        vec3d edge1 = p1 - p0;
        vec3d edge2 = p2 - p0;

        vec3d crossProduct = cross(edge1, edge2);
        return 0.5 * length(crossProduct);
    };

    double calc_spherical_triangle_area(vec3d p0, vec3d p1, vec3d p2)
    {
        coord_t num = std::fabs(dot(p0, cross(p1, p2)));

        coord_t den = 1.0 + dot(p0, p1) + dot(p1, p2) + dot(p0, p2);
        coord_t area = 2.0 * std::atan2(num, den);
        return area;
    };

    double calc_rsme_error(VoronoiDiagram &voronoi, const std::vector<double> &cell_sizes)
    {
        const Topology<Polygon> &polygons = voronoi.polygons();
        double error = 0.0;
        double sum = 0.0;
        for (int k = 0; k < polygons.n(); k++)
        {
            int group_index = polygons.group(k);
            sum += pow((voronoi.properties()[group_index].mass - cell_sizes[group_index]), 2);
            // auto props = voronoi.analyze();
            // LOG << fmt::format("cell_sizes[{}] = {}, mass = {}", group_index, cell_sizes[group_index], voronoi.properties()[group_index].mass);
        }
        error = pow((sum / polygons.n()), 0.5);

        return error;
    }

    double geogram_CVT(vec3d &site, vec3d &p1, vec3d &p2, vec3d &p3)
    {
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
        return cur_f;
    }

    void create_edge_map(VoronoiDiagram &voronoi, void *edgeMap)
    {
        const Topology<Polygon> &polygons = voronoi.polygons();

        std::unordered_map<std::pair<int, int>, int> &edgeSiteMap = *static_cast<std::unordered_map<std::pair<int, int>, int> *>(edgeMap);

        for (int k = 0; k < polygons.n(); k++)
        {
            int group_index = polygons.group(k);

            for (int j = 0; j < voronoi.polygons().length(k); j++)
            {
                // two points that make up edge
                int p1_index = polygons(k, j);
                int p2_index = polygons(k, (j + 1) % voronoi.polygons().length(k));

                edgeSiteMap.insert({std::pair(p1_index, p2_index), group_index});
            }
        }
    }

    double build_hessian(VoronoiDiagram &voronoi)
    {
        const Topology<Polygon> &polygons = voronoi.polygons();
        std::unordered_map<std::pair<int, int>, int> edgeSiteMap;
        // mapping from edge v1 - v2 to site left of v1 - v2
        create_edge_map(voronoi, static_cast<void *>(&edgeSiteMap));

        double hessian[polygons.n()][polygons.n()];

        for (int k = 0; k < polygons.n(); k++)
        {
            int group_index = polygons.group(k);
            vec3d site(voronoi.vertices()[group_index]);

            double Hii = 0.0;

            for (int j = 0; j < voronoi.polygons().length(k); j++)
            {
                // two points that make up edge
                int p1_index = polygons(k, j);
                int p2_index = polygons(k, (j + 1) % voronoi.polygons().length(k));

                vec3d p1(voronoi.vertices()[p1_index]);
                vec3d p2(voronoi.vertices()[p2_index]);

                double Aij = length(p1 - p2);

                int site2_index = edgeSiteMap[std::pair(p2_index, p1_index)];

                vec3d site2(voronoi.vertices()[site2_index]);

                double l = length(site - site2);

                hessian[group_index][site2_index] = 0.5 * Aij / l;

                Hii -= hessian[group_index][site2_index];
            }
            hessian[group_index][group_index] = Hii;
        }
    }

    template <>
    double calc_energy<SphereDomain>(unsigned int n, const double *x, double *de_dw, void *data0)
    {
        nlopt_data<SphereDomain> &data = *static_cast<nlopt_data<SphereDomain> *>(data0);
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

        TriangleQuadrature<SphericalTriangle> quad(10);

        double energy = 0.0;
        for (int k = 0; k < polygons.n(); k++)
        {
            double cell_area = 0.0;
            // get the site of the polygon
            int group_index = polygons.group(k);
            vec3d site(voronoi.vertices()[group_index]);
            double curr_weight = weights[group_index];

            const auto *p1 = voronoi.vertices()[polygons(k, 0)];

            for (int j = 1; j < voronoi.polygons().length(k) - 1; j++)
            {
                // get two vertices of the polygon
                const auto *p2 = voronoi.vertices()[polygons(k, j)];
                const auto *p3 = voronoi.vertices()[polygons(k, j + 1)];

                // calculate area of the triangle
                double t_area = calc_spherical_triangle_area(p1, p2, p3);

                // quadrature
                auto fn = [&site](const vec3d &x) -> double
                {
                    return std::pow(length(x - site), 2);
                };

                energy -= quad.integrate(fn, p1, p2, p3);
                cell_area += t_area;
            }

            if (de_dw)
            {
                de_dw[group_index] = cell_area - data.cell_sizes[group_index];
            }
            energy += (curr_weight * cell_area);
            energy -= data.cell_sizes[k] * curr_weight;
        }
        LOG << fmt::format("iter={}, energy = {}", data.iter, energy);

        return energy;
    };

    template <>
    double calc_energy<SquareDomain>(unsigned int n, const double *x, double *de_dw, void *data0)
    {
        nlopt_data<SquareDomain> &data = *static_cast<nlopt_data<SquareDomain> *>(data0);
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
            const auto *site = voronoi.vertices()[group_index];
            double curr_weight = weights[group_index];

            const auto *p1 = voronoi.vertices()[polygons(k, 0)];

            for (int j = 1; j < voronoi.polygons().length(k) - 1; j++)
            {
                // get two vertices of the polygon
                const auto *p2 = voronoi.vertices()[polygons(k, j)];
                const auto *p3 = voronoi.vertices()[polygons(k, j + 1)];

                // calculate area of the triangle
                double t_area = calc_triangle_area(p1, p2, p3);

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
        LOG << fmt::format("iter={}, energy = {}", data.iter, energy);
        return -energy;
    };
}
