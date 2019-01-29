#ifndef _CELLS_H_
#define _CELLS_H_

#include "tools.h"
#include "Array.h"

class UniformGrid;
struct Cell;


void MoveCell(int cellID, UniformGrid& Grid, const Cell* old_cells, Cell* new_cells, const int* NeighbourList, double dt, DoubleArray2D& Height, CoordArray2D& Normal, DoubleArray2D& Wall);
void GrowCell(Cell& cell, int cellID, double dt, int* dividingCells, int& numDivide, EnvArray3D& Env, AgaArray2D** Wal, UniformGrid& Grid);
void DivideCell(int parentID, int daughterID, Cell* cells, UniformGrid& Grid, const int* neighbours, DoubleArray2D& Wall, DoubleArray2D& Height, CoordArray2D& Normal, double t);

#endif // _CELLS_H_
