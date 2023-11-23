#pragma once

#include <vector>

#include "defs.h"

namespace vortex {
class Mesh;
class Vertices;

template <typename coord_t>
void get_bounding_box(const coord_t* points, int64_t n_points, int8_t dim,
                      coord_t* xmin, coord_t* xmax);

template <typename coord_t, typename index_t>
void sort_points_on_zcurve(const coord_t* points, uint64_t n_points, int8_t dim,
                           std::vector<index_t>& order);

void sample_surface(const Mesh& mesh, Vertices& sites, index_t n);
void sample_volume(const Mesh& mesh, Vertices& sites, index_t n);

}  // namespace vortex