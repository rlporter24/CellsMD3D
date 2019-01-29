#ifndef NEIGHBOURS_H_
#define NEIGHBOURS_H_

#include "tools.h"
struct Cell;
class UniformGrid;

int** InitializeNeighbourList(int maxCells, int maxNeighbours);
void getNeighbours(const Cell* cell_array, int N_cells, UniformGrid& Grid, int** NeighbourList, int maxNeighbours);
void getNeighbours(const Cell* cell_array, UniformGrid& Grid, int** NeighbourList, int maxNeighbours, int* ID, int ID_length);
int reduceNeighbours(const Cell& cell, const Cell* cell_array, int* neighbours);
double clamp(double n, double minn, double maxn);
void min_distance(const Cell& cell1, const Cell& cell2, double& d, DoubleCoord& c1, DoubleCoord& c2);


#endif /* NEIGHBOURS_H_ */
