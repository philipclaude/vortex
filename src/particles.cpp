#include "particles.h"

#include "io.h"
#include "math/mat.hpp"
#include "math/vec.hpp"

namespace vortex {

template <>
bool has_boundary<SquareDomain>() {
  return true;
}

template <>
bool has_boundary<SphereDomain>() {
  return false;
}

template <>
void project_point<SquareDomain>(double* x) {
  x[2] = 0;
}

template <>
void project_point<SphereDomain>(double* x) {
  vec3d p(x);
  p = normalize(p);
  for (int d = 0; d < 3; d++) x[d] = p[d];
}

template <>
void project_velocity<SquareDomain>(double* v) {}

template <>
void project_velocity<SphereDomain>(double* v) {
  vec3d p(v);
  vec3d u(v);
  p = normalize(p);
  vec3d ur = dot(u, p) * p;  // radial component of velocity
  vec3d ut = u - ur;         // tangential component of velocity
  for (int d = 0; d < 3; d++) v[d] = ut[d];
}

void Particles::create(size_t np, const coord_t* xp, int dim) {
  Vertices::reserve(np);
  vec3d uvw;
  centroids_.reserve(np);
  velocity_.reserve(np);
  for (size_t k = 0; k < np; k++) {
    vec4d xk(xp + dim * k, dim);
    Vertices::add(&xk[0]);
    centroids_.add(&xk[0]);
    velocity_.add(&uvw[0]);
  }
  mass_.resize(np, 1);
  volume_.resize(np, 1);
}

void Particles::save(const std::string& filename) const {
  size_t n_data = n() * 3;
  std::vector<float> data(n_data);
  size_t m = 0;
  for (size_t k = 0; k < n(); k++)
    for (int d = 0; d < 3; d++) data[m++] = (*this)[k][d];

  std::string extension = filename.substr(filename.find_last_of(".") + 1);

  if (extension == "vtk") {
    FILE* fid = fopen(filename.c_str(), "wb");
    fprintf(fid, "# vtk DataFile Version 2.0\nvortex vertices\n");
    fprintf(fid, "BINARY\nDATASET UNSTRUCTURED_GRID\nPOINTS %zu float\n", n());
    std::parafor_i(0, n_data,
                   [&data](int tid, size_t k) { io::swap_end(data[k]); });
    fwrite(&data[0], 4, n_data, fid);
    fprintf(fid, "\nCELLS 0 0\nCELL_TYPES 0\n");
    fprintf(fid, "POINT_DATA %zu\n", n());
    fprintf(fid, "FIELD FieldData 1\n");
    fprintf(fid, "density 1 %zu float\n", n());
    std::vector<float> density_data(n());
    for (size_t k = 0; k < n(); k++)
      density_data[k] = density_[k] > 1 ? 1000 : 0;
    std::parafor_i(0, n(), [&density_data](int tid, size_t k) {
      io::swap_end(density_data[k]);
    });
    fwrite(&density_data[0], 4, n(), fid);
    fprintf(fid, "\n");
    fclose(fid);
  } else if (extension == "meshb") {
    Mesh mesh(3);
    mesh.vertices().reserve(n());
    for (size_t k = 0; k < n(); k++) {
      mesh.vertices().add((*this)[k]);
    }
    meshb::write(mesh, filename);
  } else if (extension == "solb") {
    Mesh mesh(3);
    mesh.vertices().reserve(n());
    for (size_t k = 0; k < n(); k++) {
      mesh.vertices().add((*this)[k]);
    }

    // Collect density information
    std::vector<float> density_data(n());
    for (size_t k = 0; k < n(); k++)
      density_data[k] = density_[k] > 1 ? 1000 : 0;

    meshb::write(mesh, filename, false,
                 density_data);  // Pass densities to write function
  } else {
    std::cerr << "Unsupported file extension: " << extension << std::endl;
  }
}

void ParticleSimulation::compute_search_direction(
    const std::vector<double>& target_volumes) {
  hessian_.clear();

  // add regularization?
  // for (size_t k = 0; k < particles_.n(); k++)
  //  hessian_(k, k) = 1e-3 * target_volumes[k];

  // set up the sparse matrix
  for (const auto& [b, volume] : voronoi_.facets()) {
    size_t site_i = b.first;
    size_t site_j = b.second;
    vec3d pi(particles_[site_i]);
    vec3d pj(particles_[site_j]);
    double delta_ij = 0.5 * volume / length(pi - pj);
    hessian_(site_i, site_j) = delta_ij;
    hessian_(site_j, site_i) = delta_ij;
    hessian_(site_i, site_i) -= delta_ij;
    hessian_(site_j, site_j) -= delta_ij;
  }

  // solve for the search direction
  // a tolerance of 1e-3 should give a good enough direction
  hessian_.solve_nl(gradient_, dw_, 1e-3, true);
}

void ParticleSimulation::calculate_properties() {
  for (size_t k = 0; k < particles_.n(); k++) {
    ASSERT(voronoi_.properties()[k].site == k);
    particles_.volume()[k] = voronoi_.properties()[k].volume;
    ASSERT(particles_.volume()[k] > 0);
    for (int d = 0; d < 3; d++)
      particles_.centroids()(k, d) =
          voronoi_.properties()[k].moment[d] / particles_.volume()[k];
    max_displacement_[k] = std::numeric_limits<double>::max();
  }
  // voronoi_.weights().resize(particles_.n(), 0);

  // determine the maximum displacement as the min (bisector distance)/2
  const auto& facets = voronoi_.facets();
  for (const auto& [b, _] : facets) {
    size_t site_i = b.first;
    size_t site_j = b.second;
    vec3d pi(particles_[site_i]);
    vec3d pj(particles_[site_j]);
    double d = 0.5 * length(pi - pj);
    if (d < max_displacement_[site_i]) max_displacement_[site_i] = d;
    if (d < max_displacement_[site_j]) max_displacement_[site_j] = d;
  }
}

}  // namespace vortex