### About

`vortex` is a Voronoi mesher, visualizer and fluid simulator for geophysical fluids on the sphere. Currently, the main capability consists of solving the shallow water equations with a Lagrangian method to simulate the air in the Earth's atmosphere. Future work consists of applying this method to simulate the oceans.

**Note:** `vortex` currently contains new unpublished research which is why it is a private repository - please do not share the code with anyone. It could be released publicly after the ideas are published.

### Installation

To get started, `vortex` only requires `git`, `CMake` and a `C++` compiler that supports `C++17`. The `CMake` configuration will automatically download all required dependencies to the `extern` directory (which will be ignored by `git`). If you need to update a dependency (e.g. `wings`), delete the `extern/wings` directory and re-run the `CMake` configuration in your build directory.

Since `vortex` includes a visualizer built upon `wings`, it also requires `OpenGL` libraries (and `EGL` for Linux or `CGL` for macOS). `CGL` (Core `OpenGL`) should already be included in `macOS` - for Linux, the `dev/setup.sh` script will download `OpenGL`, including `EGL`. See the steps below for getting started.

#### Quickstart

To get started building `vortex`, please use the following steps:

```sh
git clone git@github.com:middleburygcl/vortex.git
cd vortex
dev/setup.sh
mkdir build
mkdir build/release
cd build/release
cmake ../../
make
```

Note that `dev/setup.sh` only needs to be run **once** on your system. Pass `-j` to the `make` command to build with multiple threads.

The steps above will build the `vortex` library in `build/release/lib/`, unit tests (in `build/release/test`), as well as the main vortex executable (`build/release/bin/vortex`).

You can also compile `vortex` in "debug" mode (which will be slower). If you configure `vortex` from `build/debug`, the `CMake` configuration will detect that you wish to build in debug mode, and will add debugging symbols when compiling. The address sanitizer can be useful for catching memory leaks.

#### Running tests

Before running the tests, you'll need to create some sample files. From the root of the repository, run the `dev/examples.sh` script by providing the location of the main `vortex` executable, as well as a mesh size for some of the sample files (here 0.02):

```sh
dev/examples.sh build/release/bin/vortex 0.02
```

Then you can navigate back to the build directory and run all the unit tests:

```sh
cd build/release
make unit
```

All unit test output (including details about possible failures) will be written to `unit_tests_output.txt`. The steps above are the same steps used in the workflow that runs on every pull request to ensure unit tests always pass.

To run a single suite of unit tests (for example the suite of Voronoi diagram tests):

```sh
make voronoi_test
```

### Examples

Many of the examples below are run by the `dev/examples.sh` script and elaborated upon in this section. All of these examples assume you are currently in the `build/release` directory.

All `vortex` functionality is accessed through the `vortex` executable (in the `bin` directory). It provides "subprograms" to run specific functions. The subprogram should be the first argument passed to `vortex`. Current subprograms include `mesh`, `extract`, `voronoi`, `merge`, `viz`, `swe` and `gm`. Each of these subprograms will have it's own set of required and optional arguments. For more information, type:

```sh
bin/vortex --help
```

Details about a specific subprogram can be printed using (e.g. for the `voronoi` subprogram):

```sh
bin/vortex voronoi --help
```

#### 1. Creating an adapted mesh of the Earth.

Some of the algorithms below rely on being able to distinguish between land and water, as well as resolve and extract the coastlines. This can be done by creating a mesh which adapts an initial icosahedron mesh to the size defined by [this image](data/oceans_2048.png). Triangles in darker regions (land) will have a target size of `hmin` and lighter regions (water) will have a target size of `hmax` (the default is `hmax = 2 * hmin`). To create the adapted mesh using the `data/oceans_2048.png` image for sizing information:

```sh
bin/vortex mesh ../../data/oceans_2048.png --hmin 0.005 --hmax 0.01 --output earth.meshb
```

You can also pass `--n_iter X` to use `X` adaptation iterations (the default is 5). Then visualize the mesh using:

```sh
bin/vortex viz earth.meshb
```

and open `vortex.html`.

At the end of the `mesh` subprogram, the background image (to determine triangle sizes) is used to assign "group" numbers to each triangle. There are two groups in the final mesh: one for water and one for land. Press the `f` key in the user interface to cycle between the fields and note that the message box will display which field is currently active.

#### 2. Extracting the coast lines and separating the mesh from Example 1 into land and water meshes.

We can directly use the group numbers in the output mesh from Example 1 to find the edges along the coast (i.e. edges whose left and right triangle groups differ). These edges can then be chained together to form closed loops which represent the coasts. This can be done with the `extract` subprogram, which also accepts arguments for where the continents and oceans mesh files should be written after the extraction:

```sh
bin/vortex extract earth.meshb --oceans water.meshb --continents land.meshb
```

#### 3. Voronoi diagram in a square.

All Voronoi diagram functionality can be accessed through the `voronoi` subprogram. The main inputs are `--domain` and `--points`. The `--domain` specifies how the domain on which the Voronoi diagram will be restricted and can be an analytic domain ("square" or "sphere" or a triangulation). The `--points` input defines how points should be initialized. For example,

```sh
bin/vortex voronoi --domain square --points random --n_points 100000 --output example3.meshb
```

This will create a Voronoi diagram in a square with 100k points and write the final mesh to `example3.meshb`.

You can also perform `X` iterations of Lloyd relaxation using the `--n_smooth` option (here, 10 iterations):

```sh
bin/vortex voronoi --domain square --points random --n_points 100000 --n_smooth 10 --output example3.meshb
```

#### 4. Voronoi diagram on a sphere.

The Voronoi diagram of a sphere is calculated in a similar way to Example 3:

```sh
bin/vortex voronoi --domain sphere --points random --n_points 100000 --n_smooth 10 --output example4.meshb
```

#### 5. Voronoi diagram on the triangulation of a sphere.

Instead of representing a surface analytically (e.g. for a sphere), we can represent it as a triangulation and calculate the Voronoi diagram restricted to this triangulation.

**Using the vertices of the triangulation for the sites:**

In this example we will use a subdivided icosahedron as the triangulation and will pass `--n_subdiv` to specify the number of subdivisions. It's also possible to pass a mesh file, and the triangles in this mesh will be used for the domain (see Example 6).

The `--points` can be initialized to the vertices of the triangulation:

```sh
bin/vortex voronoi --domain icosahedron --n_subdiv 3 --points vertices --output example5a.meshb
```

In this case, we don't specify `--n_points` since it is controlled by the vertices of the triangulation.

**Using a sampling of the triangulation for the sites:**

Alternatively, we can randomly sample points on the triangulation:

```sh
bin/vortex voronoi --domain icosahedron --n_subdiv 3 --points random --n_points 100000 --output example5b.meshb
```

When visualizing `example5b.meshb`, press `f` until the the message box shows that the "group" is plotted. Voronoi polygons that are associated with a particular site are assigned the same group.

#### 6. Voronoi diagram of the oceans using a triangulation.

The output `water.meshb` mesh of Example 2 can be used to calculate the Voronoi diagram of the oceans. Voronoi sites will be sampled along the triangulation:

```sh
bin/vortex voronoi --domain water.meshb --points random --n_points 100000 --output example6.meshb
```

#### 7. Voronoi diagram of the oceans without a triangulation.

Another idea for creating a Voronoi diagram of the oceans is to initially sample the oceans (using a texture) and then smooth the vertices. This will move some of the points into the continents and eventually converges to covering the entire sphere (this is a recent idea, so I have no theortical guarantees of this). The number of smoothing iterations needed to cover the sphere seems to depend on how many points are used:

```sh
bin/vortex voronoi --domain sphere --points random_oceans --n_points 1000000 --n_smooth 10 --output example7.meshb
```

#### 8. Merging vertices geometrically.

One limitation in `vortex` is that exact predicates are not used when calculating Voronoi diagrams. This is possible for the square and triangulation domains, but I haven't figured out a way yet for spherical domains. As such, it's possible that two Voronoi vertices have the exact same coordinates but different symbolic information (i.e. represent different Delaunay triangles). To overcome this, we can merge nearby Voronoi vertices. This is done somewhat efficiently by constructing a kd-tree and checking the nearest neighbors to the vertices - if the vertices are too close, then they are considered duplicates and element references to the duplicate will be updated.

This is all achieved using the `merge` subprogram, which can also combine polygons based on their group number. The latter (using the `--combine` option) is useful when postprocessing Voronoi diagrams of triangulations since a single Voronoi cell might be split across several triangles and, hence, produces several polygons for a single Voronoi site.

For example, using the output mesh of Example 7:

```sh
bin/vortex merge example6.meshb --combine --output example8.meshb
```

In contrast to visualizing the output meshes from Examples 5 or 6, now the "cell" field should display a single color for the Voronoi cell.

#### 9. Shallow water equation fluid simulations

`vortex` can simulate the shallow water equations using Voronoi diagrams with the `swe` subprogram. Each Voronoi cell represents a particle that moves with the fluid, and a semi-discrete optimal transport (SDOT) problem is solved at every time step to conserve the total mass of the fluid. Some of the benchmarks proposed by [Williamson et al.](https://www.sciencedirect.com/science/article/abs/pii/S0021999105800166) are supported (cases 1, 2, 5, 6). For example:

```sh
bin/vortex swe --case williamson5 --output results --particles icosahedron6 --save_every 24 --step 60
```

will use the shallow water equation solver on test case 5 (zonal flow over an isolated mountain) from the Williamson paper. The depth field will be saved to `.vtk` files every 24 hours of simulated time in the `results` directory. The number of particles can be controlled by changing the number at the end of the `icosahedron` keyword (here, the particles are initialized as the vertices obtained from subdividing an icosahedron `6` times). The time step used in the simulation here (`--step`) is 60 seconds.

The output will contain the following information:

```sh
-----------------------------------------------------------------------
| Step   | day:hr:mn | dt (s)  | Rw  | Ra  | Rm | Rp | Re | SDPD | Eh |
-----------------------------------------------------------------------
```

which refers to the time step counter, the time (in the simulation) at that time step, the time step used, and information about the simulation residuals (lower is better). Generally, `dt` will match what was prescribed by `--step` unless the time step is too big.

`Rw` indicates how well the semi-discrete optimal transport problem succeeded at finding the weights and `Ra` refers to the error in representing the true sphere area with the Voronoi diagram at that step.

`Rm`, `Rp` and `Re` refer to the residuals of mass, momentum and energy, respectively. `SDPD` refers to "simulated days per day", i.e. how many days of real-world time is the simulation achieving per day of computation time. `Eh` is the norm of the error in the height of the fluid, as compared with the analytical solution (if available for that test case).


#### 10. Gallouët-Mérigot fluid simulations

The `gm` subprogram can be used to perform fluid simulations using a technique similar to that of Gallouët-Mérigot (2017) and Lévy (2018, 2021), which also uses semi-discrete optimal transport. The output will contain the following information:

```sh
-------------------------------------------------------------------
| step | time | dt_max | dt | |rm| | nm | max(dx) | |dx| | max(v) | @X.Y steps / sec.
-------------------------------------------------------------------
```

`dt_max` is an estimate for the largest time step that can be used so that particles are always displaced somewhere within their current Voronoi cell and `dt` is the effective step size, which is the minimum of `dt_max` and the user-defined step size (controlled by `--time_step_scale`). The value of `|rm|` is the residual for the mass conservation equation (solving the SDOT problem) and `nm` is the number of Newton-iterations it took to solve this SDOT problem. The remaining values correspond to the norm of the total particle displacement from one time step to the next (and the `max`) as well as the maximum particle velocity.

Currently, the particle simulation results are saved as a series of `.vtk` files containing the particle positions and densities. Not every time step is stored - the number of time steps after which the particles are saved is controlled by the `--save_every` command-line parameter (default is every 50 time steps). The particle data can then be imported into ParaView and animated as the simulation runs.

#### A. Rayleigh-Taylor instability.

The following command will simulate a Rayleigh-Taylor instability within a rectangle in which a high-density fluid is placed on top of a lower-density one:

```sh
bin/vortex gm --domain rectangle --corners -1 -3 1 3  --n_particles 50000 --output_directory rt50k --total_time_steps 10000 --density_ratio 3 --epsilon_scale 10 --time_step_scale 0.1
```

The resulting fluid motion should look like [this](https://drive.google.com/file/d/1aGFSnAb4-hU3VrhjRZgdv3j2yqLOpPuz/view?usp=drive_link).

#### B. Rotating sphere.

This example simulates a band of high density fluid around the equator with 100k particles on a sphere. The particles are given an initial velocity of $\vec{\Omega} \times \vec{p}$ where $\vec{\Omega} = (0, \Omega_y, 0)$ and $\vec{p}$ is the initial position of a particle.
The simulation is done in a frame of reference that rotates with the sphere, so the Coriolis force is included on the right-hand-side of the momentum equation (no other external forces are applied).

```
bin/vortex gm --domain sphere  --n_particles 100000 --output_directory sphere100k --total_time_steps 10000 --density_ratio 10 --epsilon_scale 10 --time_step_scale 0.15 --omega 0.01
```

This takes approximately 1 second per time step (on a 2022 M1 MacBook Pro) so the full simulation will take about 3  hours to run. The resulting fluid motion should look like [this](https://drive.google.com/file/d/1OnLbWHb4ANTVTb80AoijXW1l_eiDheE0/view?usp=drive_link).

### Developing

The `main` branch of `vortex` is protected, so every contribution must be made in the form of a pull request (PR). Please create a branch (`git checkout -b [branchname]`) when developing a new feature.

Upon every PR, `vortex` will check (using a GitHub workflow) that the unit tests pass, and will also check that the format is in the correct style - this project uses the `Google` style for `C++`. The `clang-format` tool (downloaded by `dev/setup.sh` mentioned above) is used to check the style, which can also be used to update files to conform to the style. I strongly recommend using VS Code and in the Settings, set the **C_Cpp: Clang_format_fallback Style** to **Google** and set the **C_Cpp: Formatting** to **clangFormat**. In the **Text Editor -> Formatting** section, I would also recommend enabling the **Format On Save** option.

#### Programming with `vortex`

The first thing to note is the use of `index_t` and `coord_t` (defined in `src/defs.h`) for defining the integer type to use for element indices, and the floating-point type to use for vertex coordinates, respectively.

The main container for elements and vertices is the `Mesh` structure, which stores all the possible mesh topologies. Each topology is of type `Topology<T>` where the template parameter `T` can be `Line`, `Triangle`, `Quad` or `Polygon`. These topologies inherit from `array2d<index_t>`, thus representing how vertices are connected to form elements. The number of entities in a topology is accessed using `.n()` and the `[]` operator will return a pointer to the first item stored in an entity. For example,

```c++
Mesh mesh(3); // initialize a mesh with 3d vertices

// create 20 random vertices
size_t n_vertices = 20;
coord_t x[3];
for (size_t k = 0; k < n_vertices; k++) {
   for (int d = 0; d < 3; d++)
      x[d] = coord_t(rand()) / coord_t(RAND_MAX);
   mesh.vertices().add(x);
}
ASSERT(mesh.vertices().n() == n_vertices); // use of ASSERT
LOG << fmt::format("mesh has {} vertices", mesh.vertices().n()); // use of LOG and fmt::format

// retrieve y-coordinate of vertex 5
LOG << fmt::format("y[5] = {}", mesh.vertices()[5][1]);

auto& polygons = mesh.polygons(); // retrieve a reference to the polygons
ASSERT(polygons.n() == 0);

index_t polygon[5] = {0, 4, 12, 7, 2}; // some polygon
polygons.add(polygon, 5);

ASSERT(polygons.n() == 1);
ASSERT(polygons.length(0) == 5); // there are 5 items in the first element (polygon)
ASSERT(polygons[0][3] == 7); // see the polygon indices above
```

A half-edge representation of the mesh can be created using:

```c++
HalfMesh hmesh(mesh);

for (auto& node : hmesh.nodes()) {} // loop through the HalfNode objects
for (auto& edge : hmesh.edges()) {} // loop through the HalfEdge objects
for (auto& face : hmesh.faces()) {} // loop through the HalfFace objects
```

Note that the half-edge mesh will be created using **all** polygons, triangles and quads. If you are using an output mesh from the Voronoi diagram calculation, you should either (1) clear the triangles with `mesh.triangles().clear()` (if you wish to work with the polygons) or (2) copy the polygons to a new mesh (using `mesh.polygons().copy(mesh2.polygons())`). Remember to copy the vertices as well if you choose to do the latter (`mesh.vertices().copy(mesh2.vertices())`).

Please see the member functions and variables for `HalfNode`, `HalfEdge` and `HalfFace` in `src/halfedges.h`. One way to become familiar with this data structure is to study the `HalfMesh::flip`, `HalfMesh::split` or `HalfMesh::collapse` functions, in addition to how they are used in `mesher.cpp`. It may also be useful to look at the `run_extract` function in `src/vortex.cpp`.

#### Checking code coverage.

Code coverage reports let us know how many lines of the `vortex` source code are "covered" by the testing suite. Generating these reports requires a few extra dependencies (`gcov`, `lcov`, and `genhtml` which requires the Perl `GD.pm` module), which are installed by the `dev/setup.sh` script.

Assuming you are starting from the root of the `vortex` repository, create a directory that will automatically be detected for code coverage:

```sh
mkdir build/coverage
```

Navigate to this directory, configure and build:

```sh
cd build/coverage
cmake ../../
make unit_coverage
```

This time, we are running the `unit_coverage` target which will run all of the unit tests and generate code coverage information. Now, display the code coverage results:

```sh
make coverage_show
```

and navigate through the `HTML` report to see which lines are being executed.

### Acknowledgements

Many thanks to the following projects which `vortex` depends on:

- `argparse`: https://github.com/p-ranav/argparse
- `fmtlib`: https://github.com/fmtlib/fmt
- `libMeshb`: https://github.com/LoicMarechal/libMeshb
- `morton-nd`: https://github.com/morton-nd/morton-nd
- `OpenNL`: https://github.com/BrunoLevy/geogram.psm.OpenNL
- `PCK`: https://github.com/BrunoLevy/geogram.psm.Predicates
- `stb` (via `wings`): https://github.com/nothings/stb
- `tinyobjloader`: https://github.com/tinyobjloader/tinyobjloader

### License

All `vortex` source code (`C++`, `HTML`, `JavaScript` and `GLSL`) is distributed under the Apache-2.0 License.

Copyright 2023 - 2025 Philip Claude Caplan

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

### Contributors

- Philip Caplan ([philipclaude](https://github.com/philipclaude))
- Otis Milliken ([Hokalaka2](https://github.com/Hokalaka2))
- Col McDermott ([Col-McDermott](https://github.com/Col-McDermott))
- Zeyi Tong ([zeyiiiii](https://github.com/zeyiiiii))
- Toby Pouler ([tpouler](https://github.com/tpouler))