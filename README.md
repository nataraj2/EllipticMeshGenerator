# Elliptic mesh generator for 3D geometries

## Compilation and running

mpicxx -std=c++14 EllipticMeshGenerator.cpp -o out
./out

## Inputs


## Outputs
1. The code outputs a mesh `file_mesh.vtk` which is a ascii VTK file of the structured mesh 
2. A `bc.dat` file which contains the boundary conditions for the bottom wall
3. Files `file_bc_wall.vtk` and `file_bc_farfield.vtk` which are vtk files to visualize 
the wall bc regions for the hump region vicinity


## The smoothing algorithm
The elliptic governing equation is solved for the y coordinate to smooth  

$\alpha_{11}\cfrac{\partial^2 y}{\partial \xi^2} + 2\alpha_{12}\cfrac{\partial^2 y}{\partial\xi\partial\eta} + 2\alpha_{13}\cfrac{\partial^2y}{\partial\xi\partial\zeta}+\alpha_{22}\cfrac{\partial^2y}{\partial\eta^2}+2\alpha_{23}\cfrac{\partial^2y}{\partial\eta\partial\zeta}+\alpha_{33}\cfrac{\partial^2y}{\partial\zeta^2} = 0$ 




