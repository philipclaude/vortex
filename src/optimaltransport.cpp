#include "optimaltransport.h"
#include "math/mat.hpp"
#include "math/vec.hpp"
#include "quadrature.hpp"
#include <iostream>
#include <fstream>

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

    // https : // www.johndcook.com/blog/2021/11/29/area-of-spherical-triangle/
    double
    calc_spherical_triangle_area(vec3d p0, vec3d p1, vec3d p2)
    {
        double num = std::fabs(dot(p0, cross(p1, p2)));

        double den = 1.0 + dot(p0, p1) + dot(p1, p2) + dot(p0, p2);
        double area = 2.0 * std::atan2(num, den);
        return area;
    };

    double calc_spherical_distance(vec3d p1, vec3d p2)
    {
        double num = dot(p1, p2);
        double distance = std::acos(num);
        return distance;
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
        }
        error = pow((sum / polygons.n()), 0.5);

        return error;
    }

    double calc_gradient_norm(std::vector<double> &de_dw)
    {
        double sum = 0.0;
        for (int i = 0; i < de_dw.size(); i++)
        {
            sum += pow(de_dw[i], 2);
        }
        return pow(sum, 0.5);
    }

    double calc_gradient(std::vector<double> &de_dw)
    {
        double sum = 0.0;
        for (int i = 0; i < de_dw.size(); i++)
        {
            sum += de_dw[i];
        }
        return sum;
    }

    double
    geogram_CVT(vec3d site, vec3d p1, vec3d p2, vec3d p3)
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

    void create_edge_map(VoronoiDiagram &voronoi, std::unordered_map<std::pair<int, int>, int> &edgeSiteMap)
    {
        const Topology<Polygon> &polygons = voronoi.polygons();

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

    void build_hessian(VoronoiDiagram &voronoi, spmat<double> &hessian, std::vector<double> &de_dw, std::vector<double> &cell_area)
    {
        const Topology<Polygon> &polygons = voronoi.polygons();
        std::unordered_map<std::pair<int, int>, int> edgeSiteMap;

        // mapping from edge v1 - v2 to site left of v1 - v2
        create_edge_map(voronoi, edgeSiteMap);

        for (int k = 0; k < polygons.n(); k++)
        {
            int group_index = polygons.group(k);
            vec3d site(voronoi.vertices()[group_index]);

            double Hii = 0.0;
            double area = 0.0;
            const auto *p0 = voronoi.vertices()[polygons(k, 0)];
            for (int j = 0; j < voronoi.polygons().length(k); j++)
            {
                // two points that make up edge
                int p1_index = polygons(k, j);
                int p2_index = polygons(k, (j + 1) % voronoi.polygons().length(k));

                vec3d p1(voronoi.vertices()[p1_index]);
                vec3d p2(voronoi.vertices()[p2_index]);

                double Aij = calc_spherical_distance(p1, p2);

                int site2_index = edgeSiteMap[std::pair(p2_index, p1_index)];

                vec3d site2(voronoi.vertices()[site2_index]);

                double l = length(site - site2);

                double entry = 0.5 * Aij / l;
                hessian(group_index, site2_index) = entry;

                Hii -= entry;
                if (j != 0 && j != voronoi.polygons().length(k) - 1)
                {
                    double t_area = calc_spherical_triangle_area(p0, p1, p2);
                    area += t_area;
                }
            }
            de_dw[group_index] = cell_area[group_index] - area;

            hessian(group_index, group_index) = Hii;
        }
    }

    double calculate_energy(VoronoiDiagram &voronoi, std::vector<double> &ideal_size, std::vector<double> &de_dw, std::set<int> &sites_visited)
    {
        TriangleQuadrature<SphericalTriangle> quad(4);
        const Topology<Polygon> &polygons = voronoi.polygons();
        auto &weights = voronoi.weights();

        double energy = 0.0;
        for (int k = 0; k < polygons.n(); k++)
        {
            double cell_area = 0.0;
            double cell_energy = 0.0;

            // get the site of the polygon
            int group_index = polygons.group(k);
            vec3d site(voronoi.vertices()[group_index]);
            sites_visited.erase(group_index);
            double curr_weight = weights[group_index];

            const auto *p1 = voronoi.vertices()[polygons(k, 0)];

            for (int j = 1; j < voronoi.polygons().length(k) - 1; j++)
            {
                // get two vertices of the polygon
                const auto *p2 = voronoi.vertices()[polygons(k, j)];
                const auto *p3 = voronoi.vertices()[polygons(k, j + 1)];

                // calculate area of the triangle
                double t_area = calc_spherical_triangle_area(p1, p2, p3);

                auto fn = [&site](const vec3d &q) -> double
                {
                    return std::pow(length(q - site), 2);
                };

                cell_energy += quad.integrate(fn, p1, p2, p3);

                cell_area += t_area;
            }
            if (!de_dw.empty())
            {
                de_dw[group_index] = ideal_size[group_index] - cell_area;
            }
            energy += cell_energy + curr_weight * (-cell_area + ideal_size[group_index]);
        }

        // calculate for any missed vertices
        for (int index : sites_visited)
        {
            if (!de_dw.empty())
            {
                de_dw[index] = ideal_size[index];
            }
            energy += weights[index] * ideal_size[index];
        }

        for (int i = 0; i < de_dw.size(); i++)
        {
            de_dw[i] = -1 * de_dw[i];
        }

        return -energy;
    }

    template <>
    double objective_func<SphereDomain>(const std::vector<double> &x, std::vector<double> &de_dw, void *data0)
    {
        Timer time;
        nlopt_data<SphereDomain> &data = *reinterpret_cast<nlopt_data<SphereDomain> *>(data0);
        VoronoiDiagram &voronoi = data.voronoi;

        time.start();
        auto &weights = voronoi.weights();
        std::set<int> sites_visited;
        weights.resize(x.size(), 0.0);
        for (int i = 0; i < x.size(); i++)
        {
            sites_visited.insert(i);
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

        double energy = calculate_energy(voronoi, data.cell_sizes, de_dw, sites_visited);

        time.stop();

        data.energy_time += time.seconds();

        // double error = calc_gradient_norm(de_dw);
        // data.error = error;

        if (data.output_converge)
        {
            if (data.outputFile.is_open())
            {
                // data.outputFile << "iter: " << data.iter << " energy change: " << abs(energy - data.energy) << std::endl;
                data.outputFile << "iter: " << data.iter << " error: " << abs(energy - data.energy) << std::endl;
            }
            else
            {
                std::cout << "Error opening file" << std::endl;
            }
        }
        data.energy = energy;
        return energy;
    };

    template <>
    double objective_func<SquareDomain>(const std::vector<double> &x, std::vector<double> &de_dw, void *data0)
    {
        Timer time;
        nlopt_data<SquareDomain> &data = *reinterpret_cast<nlopt_data<SquareDomain> *>(data0);
        VoronoiDiagram &voronoi = data.voronoi;
        const Topology<Polygon> &polygons = voronoi.polygons();
        LOG << polygons.n();
        time.start();
        auto &weights = voronoi.weights();
        std::set<int> sites_visited;
        weights.resize(x.size(), 0.0);
        for (int i = 0; i < x.size(); i++)
        {
            sites_visited.insert(i);
            weights[i] = x[i];
        }
        lift_sites(data.vertices, voronoi.weights());

        data.options.store_mesh = true;

        voronoi.vertices().clear();
        voronoi.vertices().set_dim(3);
        voronoi.polygons().clear();
        voronoi.triangles().clear();
        voronoi.compute(data.domain, data.options);

        data.iter++;

        double energy = 0.0;
        for (int k = 0; k < polygons.n(); k++)
        {
            double cell_area = 0.0;
            double cell_energy = 0.0;

            // get the site of the polygon
            int group_index = polygons.group(k);
            vec3d site(voronoi.vertices()[group_index]);
            sites_visited.erase(group_index);
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
                cell_energy += t_area * cur_f / 6.0;

                cell_area += t_area;
            }
            if (!de_dw.empty())
            {
                de_dw[group_index] = data.cell_sizes[group_index] - cell_area;
            }
            energy += cell_energy + curr_weight * (-cell_area + data.cell_sizes[group_index]);
        }

        // calculate for any missed vertices
        for (int index : sites_visited)
        {
            if (!de_dw.empty())
            {
                de_dw[index] = data.cell_sizes[index];
            }
            energy += weights[index] * data.cell_sizes[index];
        }

        for (int i = 0; i < de_dw.size(); i++)
        {
            de_dw[i] = -1 * de_dw[i];
        }

        time.stop();

        data.energy_time += time.seconds();

        double error = calc_gradient_norm(de_dw);
        data.error = error;

        if (data.output_converge)
        {
            if (data.outputFile.is_open())
            {
                data.outputFile << "iter: " << data.iter << " energy change: " << abs(-energy - data.energy) << std::endl;
            }
            else
            {
                std::cout << "Error opening file" << std::endl;
            }
        }
        data.energy = -energy;

        return -energy;
    };
}
