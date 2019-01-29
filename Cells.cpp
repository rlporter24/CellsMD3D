#include "Cells.h"

#include "Compute.h"
#include "Forces.h"
#include "grow.h"
#include "Integrate.h"
#include "Constants.h"
#include "UniformGrid.h"
#include "Cell.h"
#include <omp.h>


// Main function to move cell
void MoveCell(int cellID, UniformGrid& Grid, const Cell* old_cells, Cell* new_cells, const int* neighbours, double dt, DoubleArray2D& Height, CoordArray2D& Normal, DoubleArray2D& Wall)
{
	UniformGrid::Address oldAddress, newAddress;	// for storing addresses of cells in the Grid

	DoubleCoord v, Fnet;
	IntCoord XYAddress;

	const Cell& oldCell = old_cells[cellID];
	Cell& newCell = new_cells[cellID];

	// gives the current neighbours of the cell
	oldAddress = Grid.GetAddress(average(oldCell.Position));
	XYAddress = Grid.GetXY(oldAddress);

	// integrates one step, updates positions from old to new
	integrate(dt, cellID, old_cells, new_cells, neighbours, Height, Normal, Grid, XYAddress, Wall);

	// check if the cell has moved out of its box
	newAddress = Grid.GetAddress(average(newCell.Position));
    if (newAddress.a!=oldAddress.a) {
#pragma omp critical
        {
		Grid.Move(cellID, oldAddress, newAddress);
        }
    }
}

// Main function to grow cell
void GrowCell(Cell& cell, int cellID, double dt, int* dividingCells, int& numDivide, EnvArray3D& Env, AgaArray2D** Wal, UniformGrid& Grid)
{
	// grow new cells
	grow(dt, cell, Env, Wal, Grid);
    int index;
	// check if cell will divide
	if (cell.Length>L_divide)
	{
#pragma omp critical
        {
        index = numDivide++;
        }
		dividingCells[index] = cellID;
	}

}

// Main function to divide cell
void DivideCell(int parentID, int daughterID, Cell* cells, UniformGrid& Grid, const int* neighbours, DoubleArray2D& Wall, DoubleArray2D& Height, CoordArray2D& Normal, double t)
{
	UniformGrid::Address oldAddress, newAddress;	// for storing addresses of cells in the Grid
	//double F_centre;
	DoubleCoord Fnet;
	Tensor stressTensor;
	Cell& parentCell = cells[parentID];
	Cell& daughterCell = cells[daughterID];

	// find stress on the mother cell
	oldAddress = Grid.GetAddress(average(parentCell.Position));
	
	// remove the ID from the grid
	Grid.Remove(parentID, oldAddress);

	// divide and create a new cell with ID N_cells
	divide(parentCell, daughterCell, t);

	// add mother and daughter to grid
	Grid.Add(parentID, Grid.GetAddress(average(parentCell.Position)));
	Grid.Add(daughterID, Grid.GetAddress(average(daughterCell.Position)));
}
