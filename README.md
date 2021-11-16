Written by Hyunmin Yang, HANUL, Korea University.

This program simulates tripple GEMs using Garfield++, ElmerFem and Gmsh.

# Prerequisites
ROOT version 6  
Garfield++  
Geant4  
Gmsh  
Elmer

# Executables
tgem               : Plot equipotential lines and drift lines of tripple GEMs.

# Garfield++
## Meshing
First, you have to generate a 3D mesh file describing an unit cell of tripple GEM.
```
mkdir build
cd build
cmake ..
gmsh tgemcell.geo -3 -order 2
```
**tgemcell.msh** will be created on your directory.
## Generating field map data.
```
ElmerGrid 14 2 tgemcell.msh -autoclean
ElmerSolver tgemcell.sif
```
Running ElmerGrid converts .msh file to **mesh.boundary**, **mesh.elements**, **mesh.header** and **mesh.nodes** in the tgemcell diretory. The second command solves the field equation with **tgemcell.sif**(start information file).
**tgemcell.ep** and **tgemcell.result** will be created the tgemcell directory.

## Importing field map data to Garfield++
Following all steps above, required files are ready in the tgemcell directory.
1. **tgemcell/dielectrics.dat** : First line tells the number of material(different volumes with same material are considered to have different material ID.) and the list of material IDs and its relative permitivity.
2. **mesh.boundary**, **mesh.elements**, **mesh.header** and **mesh.nodes** : geometry information.
3. **tgemcell/tgemcell.result** : field map solution.
