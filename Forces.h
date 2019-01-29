#ifndef FORCES_H_
#define FORCES_H_
#include "tools.h"
#include "UniformGrid.h"
#include "Array.h"
#include "Cell.h"

void F_cc(const Cell& cell1, const Cell& cell2, DoubleCoord& F, DoubleCoord& r, double& d);
void F_cw(const Cell& cell, double Wall, DoubleCoord& F1, DoubleCoord& F2, DoubleCoord& r1, DoubleCoord& r2, DoubleCoord& cwStaFric, DoubleCoord& cwDynFric);
void F_surf_tension(const Cell& cell, UniformGrid& Grid, const IntCoord& XYAddress, DoubleArray2D& Height, CoordArray2D& Normal, DoubleCoord& F, DoubleCoord& T);
void F_v(const Cell& cell, DoubleCoord& F, DoubleCoord& T);
double clamp(double n, double minn, double maxn);
void min_distance(const Cell& cell1, const Cell& cell2, double& t, DoubleCoord& c1, DoubleCoord& c2); // double& d, double& s,
void sum_forces(const Cell& cell, const Cell* cell_array, const int* neighbours, DoubleCoord& Fnet, DoubleCoord& Tnet, DoubleArray2D& Height, CoordArray2D& Normal, UniformGrid& Grid, const IntCoord& XYAddress, DoubleArray2D& Wall, DoubleCoord& cwStaFric, DoubleCoord& cwDynFric);
void reduce_neighbours(const Cell& cell, const Cell* cell_array, int* neighbours, int& numNeighbours);
//void F_int(double k, const Cell& cell, Segment& F);

#endif /* FORCES_H_ */

