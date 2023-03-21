# Elliptic mesh generator for 3D geometries

## Compilation and running

mpicxx -std=c++14 EllipticMeshGenerator.cpp -o out
./out

## Inputs


## The smoothing algorithm
The elliptic governing equation is solved for the y coordinate to smooth 
$\alpha\frac{\partial y}{\partial \xi}^2$




