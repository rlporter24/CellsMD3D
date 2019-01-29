#ifndef COMPUTE_H_
#define COMPUTE_H_

#include "UniformGrid.h"
#include "tools.h"
#include "Array.h"
#include "Cell.h"

void mean_stress(const Cell& cell, const Cell* cell_array, const int* neighbours, UniformGrid& Grid, DoubleArray2D& Wall, DoubleArray2D& Height, CoordArray2D& Normal, Tensor& stressTensor, DoubleCoord& Fnet);
void GetDensity(DoubleArray3D& DensityVec, DoubleArray3D& Density1Vec, DoubleArray3D& Density2Vec, DoubleArray3D& insideColonyDen,DoubleArray2D& Height, UniformGrid& Grid, Cell* cells, int& minx, int& maxx, int& miny, int& maxy, int& maxz);
void GetHeight(DoubleArray3D& DensityVec, DoubleArray2D& Height, UniformGrid& Grid, Cell* cells, int& minx, int& maxx, int& miny, int& maxy, int& maxz);
void ShiftDensity(DoubleArray3D& DensityVec, DoubleArray3D& Density1Vec, DoubleArray3D& Density2Vec,DoubleArray2D& DensityWallVec, DoubleArray2D& DensityWallVec1, DoubleArray2D& DensityWallVec2, DoubleArray3D& DensityVecShiftP, DoubleArray3D& Density1VecShiftP, DoubleArray3D& Density2VecShiftP,DoubleArray2D& DensityWallVecShiftP, DoubleArray2D& DensityWallVec1ShiftP, DoubleArray2D& DensityWallVec2ShiftP,int BoxX, int BoxY, int BoxZ);
void GetSurfaceNormal(CoordArray2D& Normal, DoubleArray2D& Height, int minX, int maxX, int minY, int maxY);
void BoxAverage(DoubleArray3D& RoughDensity, DoubleArray3D& Density, DoubleArray2D& RoughWallDensity, DoubleArray2D& WallDensity, DoubleArray3D& Filter, int FilterDim, int minx, int maxx, int miny, int maxy, int maxz, int BoxX, int BoxY, int BoxZ);
//void GetHeight(DoubleArray2D& Height, UniformGrid& Grid, Cell* cells, int& minx, int& maxx, int& miny, int& maxy, int& maxz);
void Height_Average(DoubleArray2D& InHeight, DoubleArray2D& OutHeight, int minx, int maxx, int miny, int maxy);
void Smooth(DoubleArray2D& InArray, DoubleArray2D& OutArray, DoubleArray2D& Filter, int FilterDim, int minx, int maxx, int miny, int maxy);

#endif /* COMPUTE_H_ */
