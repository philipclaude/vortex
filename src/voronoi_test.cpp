#include "voronoi.h"

#include <Predicates_psm.h>

#include "graphics.h"
#include "io.h"
#include "library.h"
#include "tester.h"
#include "util.h"

using namespace vortex;

UT_TEST_SUITE(voronoi_test_suite)

UT_TEST_CASE(test_pool) {
  // basic tests
  {
    std::vector<int> x_pool(300);

    pool<int> x(x_pool.data(), x_pool.capacity());
    UT_ASSERT_EQUALS(x.size(), 0);
    UT_ASSERT(x.empty());
    x.push_back(1);
    UT_ASSERT_EQUALS(x.capacity(), 300);
    UT_ASSERT_EQUALS(x.size(), 1);

    for (int i = 0; i < 200; i++) {
      x.push_back(i + 2);
    }
    UT_ASSERT_EQUALS(x.size(), 201);
    UT_ASSERT_EQUALS(x.capacity(), 300);

    for (size_t i = 0; i < x.size(); i++) {
      UT_ASSERT_EQUALS(size_t(x[i]), i + 1);
    }

    UT_ASSERT_EQUALS(x.back(), 201);

    x.clear();
    UT_ASSERT(x.empty());
    UT_ASSERT_EQUALS(x.capacity(), 300);
  }

  // swap test
  {
    std::vector<int> x_pool(3), y_pool(4);
    pool<int> x(x_pool.data(), 4);
    pool<int> y(y_pool.data(), 4);
    x.set_size(3);
    y.set_size(4);
    x[0] = 1;
    x[1] = 4;
    x[2] = 9;
    y[0] = 2;
    y[1] = 10;
    y[2] = 0;
    y[3] = 5;
    x.swap(y);

    UT_ASSERT_EQUALS(x.size(), 4);
    UT_ASSERT_EQUALS(y.size(), 3);
    UT_ASSERT_EQUALS(x[0], 2);
    UT_ASSERT_EQUALS(x[1], 10);
    UT_ASSERT_EQUALS(x[2], 0);
    UT_ASSERT_EQUALS(x[3], 5);
    UT_ASSERT_EQUALS(y[0], 1);
    UT_ASSERT_EQUALS(y[1], 4);
    UT_ASSERT_EQUALS(y[2], 9);
  }
}
UT_TEST_CASE_END(test_pool)

UT_TEST_CASE(test_sphere) {
  auto irand = [](int min, int max) {
    return min + double(rand()) / (double(RAND_MAX) + 1.0) * (max - min);
  };
  static const int dim = 3;
  size_t n_sites = 1e4;
  std::vector<coord_t> sites(n_sites * dim, 0.0);
  for (size_t k = 0; k < n_sites; k++) {
    coord_t theta = 2.0 * M_PI * irand(0, 1);
    coord_t phi = acos(2.0 * irand(0, 1) - 1.0);
    sites[k * dim + 0] = cos(theta) * sin(phi);
    sites[k * dim + 1] = sin(theta) * sin(phi);
    sites[k * dim + 2] = cos(phi);
  }

  std::vector<index_t> order(n_sites);
  sort_points_on_zcurve(sites.data(), n_sites, dim, order);

  Vertices vertices(dim);
  vertices.reserve(n_sites);
  coord_t x[dim];
  for (size_t i = 0; i < n_sites; i++) {
    for (int d = 0; d < dim; d++) x[d] = sites[dim * order[i] + d];
    vertices.add(x);
  }

  SphereDomain domain(1.0);
  VoronoiDiagram voronoi(dim, vertices[0], n_sites);
  VoronoiDiagramOptions options;
  options.n_neighbors = 75;
  options.allow_reattempt = false;
  options.parallel = true;
  int n_iter = 10;
  for (int iter = 1; iter <= n_iter; ++iter) {
    options.store_mesh = iter == n_iter;
    options.verbose = (iter == 1 || iter == n_iter - 1);
    voronoi.vertices().clear();
    voronoi.vertices().set_dim(3);
    voronoi.polygons().clear();
    voronoi.compute(domain, options);

    // move each site to the centroid of the corresponding cell
    const auto& properties = voronoi.properties();
    ASSERT(properties.size() == n_sites);
    vec3 x;
    double area = 0.0;
    for (size_t k = 0; k < n_sites; k++) {
      x = static_cast<float>(1.0 / properties[k].mass) * properties[k].moment;
      x = unit_vector(x);
      for (int d = 0; d < 3; d++) vertices[k][d] = x[d];
      area += properties[k].mass;
    }
    LOG << fmt::format("iter = {}, area = {}", iter, area);
  }

  voronoi.fields().set_defaults(voronoi);
  // Viewer viewer(voronoi, 7681);
  if (voronoi.polygons().n() > 0) meshb::write(voronoi, "sphere.meshb");
}
UT_TEST_CASE_END(test_sphere)

UT_TEST_CASE(test_sphere_triangulation) {
  auto irand = [](int min, int max) {
    return min + double(rand()) / (double(RAND_MAX) + 1.0) * (max - min);
  };
  Sphere sphere(7);
  static const int dim = 3;
  size_t n_sites = 1e6;
#if 0
  std::vector<coord_t> sites(n_sites * dim, 0.0);
  for (size_t k = 0; k < n_sites; k++) {
    coord_t theta = 2.0 * M_PI * irand(0, 1);
    coord_t phi = acos(2.0 * irand(0, 1) - 1.0);
    sites[k * dim + 0] = cos(theta) * sin(phi);
    sites[k * dim + 1] = sin(theta) * sin(phi);
    sites[k * dim + 2] = cos(phi);
  }
#else
  Vertices data(3);
  sample_surface(sphere, data, n_sites);
  auto& sites = data.data();
#endif

  std::vector<index_t> order(n_sites);
  sort_points_on_zcurve(sites.data(), n_sites, dim, order);

  Vertices vertices(dim);
  vertices.reserve(n_sites);
  coord_t x[dim];
  for (size_t i = 0; i < n_sites; i++) {
    for (int d = 0; d < dim; d++) x[d] = sites[dim * order[i] + d];
    vertices.add(x);
  }

  TriangulationDomain domain(sphere.vertices()[0], sphere.vertices().n(),
                             sphere.triangles()[0], sphere.triangles().n());
  VoronoiDiagram voronoi(dim, vertices[0], n_sites);
  VoronoiDiagramOptions options;
  options.n_neighbors = 75;
  options.allow_reattempt = false;
  options.parallel = true;
  int n_iter = 1;
  for (int iter = 1; iter <= n_iter; ++iter) {
    options.store_mesh = iter == n_iter;
    options.verbose = (iter == 1 || iter == n_iter - 1);
    voronoi.vertices().clear();
    voronoi.vertices().set_dim(3);
    voronoi.polygons().clear();
    voronoi.compute(domain, options);

    // move each site to the centroid of the corresponding cell
    const auto& properties = voronoi.properties();
    ASSERT(properties.size() == n_sites);
    vec3 x;
    double area = 0.0;
    for (size_t k = 0; k < n_sites; k++) {
      x = static_cast<float>(1.0 / properties[k].mass) * properties[k].moment;
      x = unit_vector(x);
      for (int d = 0; d < 3; d++) vertices[k][d] = x[d];
      area += properties[k].mass;
    }
    LOG << fmt::format("iter = {}, area = {}", iter, area);
  }

  // randomize the colors a bit, otherwise neighboring cells
  // will have similar colors and won't visually stand out

  size_t n_colors = 20;
  std::vector<int> site2color(n_sites);
  for (size_t k = 0; k < n_sites; k++)
    site2color[k] = int(n_colors * double(rand()) / double(RAND_MAX));
  for (size_t k = 0; k < voronoi.polygons().n(); k++) {
    int group = voronoi.polygons().group(k);  // the group is the site
    voronoi.polygons().set_group(k, site2color[group]);
  }

  LOG << fmt::format("writing {} polygons", voronoi.polygons().n());
  // voronoi.fields().set_defaults(voronoi);
  //  Viewer viewer(voronoi, 7681);
  if (voronoi.polygons().n() > 0)
    meshb::write(voronoi, "sphere_triangulation.meshb");
}
UT_TEST_CASE_END(test_sphere_triangulation)

UT_TEST_SUITE_END(voronoi_test_suite)