Octree-mg
==

The octree-mg library implements parallel geometric multigrid methods on
quadtree/octree grids, which can be used to solve elliptic PDEs such as
Poissons's equation. The provided solvers can be used in existing
adaptive-mesh-refinement (AMR) frameworks that employ quadtree/octree grids.

## Usage

Type `make` in the top folder, and run the programs in the `tests` folder.

## Features

* MPI parallelization
* The same code can be used in 2D/3D (compiled into two different libraries)
* Support for adaptively refined grids, with a consistent discretization near refinement boundaries
* Rectangular grids can easily be created (e.g. 512 x 256 x 256 cells).
* Operators with sparse stencils (5/7-point in 2D/3D) are supported
* Support for periodic, Dirichlet, Neumann, and continuous boundary conditions
* Coarse grids are automatically created

## Restrictions

* One layer of ghost cells is used, and diagonal ghost cells are currently not set. This
  means stencils with diagonal elements are not possible.
* Point-wise smoothers are employed, so currently a requirement is that dx, dy and dz are similar

## TODO

* Add test with refinement boundaries (is tested elsewhere already)
* Provide better load balancing for stand-alone usage
