# Elliptic mesh generator for 3D geometries

## Compilation and running

mpicxx -std=c++14 EllipticMeshGenerator.cpp -o out
./out

## Inputs


## The smoothing algorithm
The elliptic governing equation is solved for the y coordinate to smooth    
$\alpha_{11}\frac{\partial^2 y}{\partial \xi^2} + 2\alpha_{12}\frac{\partial^2 y{\partial\xi\partial\eta} + 2\alpha_{13}\frac{\partial^2y}{\partial\xi\partial\zeta}$
$+\alpha_{22}\frac{\partial^2y}{\partial\eta^2}+2\alpha_{23}\frac{\partial^2y}{\partial\eta\partial\zeta}+\alpha_{33}\frac{\partial^2y}{\partial\zeta^2}$ 




