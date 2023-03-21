# Elliptic mesh generator for 3D geometries

## Compilation and running

mpicxx -std=c++14 EllipticMeshGenerator.cpp -o out
./out

## Inputs

## The smoothing algorithm
The elliptic governing equation is solved for the y coordinate to smooth  

$\alpha_{11}\cfrac{\partial^2 y}{\partial \xi^2} + 2\alpha_{12}\cfrac{\partial^2 y}{\partial\xi\partial\eta} + 2\alpha_{13}\cfrac{\partial^2y}{\partial\xi\partial\zeta}+\alpha_{22}\cfrac{\partial^2y}{\partial\eta^2}+2\alpha_{23}\cfrac{\partial^2y}{\partial\eta\partial\zeta}+\alpha_{33}\cfrac{\partial^2y}{\partial\zeta^2} = 0$ 




