#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <stdio.h>

#include "UniformGrid.h"
#include "tools.h"
#include "Array.h"
#include "Cell.h"

struct OutputFiles;

void RunSimulation(int N_cells, Cell* old_cells, Cell* new_cells, int** NeighbourList, int maxNeighbours, UniformGrid& Grid, OutputFiles Files, bool append, DoubleArray2D& Height, DoubleArray3D& Density, DoubleArray3D& Density1, DoubleArray3D& Density2, DoubleArray2D& WallDensity,DoubleArray2D& WallDensity1,DoubleArray2D& WallDensity2,EnvArray3D& Environment, EnvArray3D& oldEnvironment, AgaArray3D** FieldAgar, AgaArray3D** oldFieldAgar, AgaArray2D** FieldWall, AgaArray2D** oldFieldWall, CoordArray2D& Normal);

#endif // _SIMULATION_H_
