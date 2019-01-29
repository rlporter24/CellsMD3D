#include <omp.h>
#include <iostream>
#include <math.h>
#include <time.h>
//#include <clapack.h>
#include "tools.h"
#include "UniformGrid.h"
#include "InputOutput.h"
#include "Constants.h"
#include "Simulation.h"
#include "Array.h"
#include "Cell.h"
#include "Neighbours.h"

using namespace std;

int main( int argc, char* argv[] ) 
{
	srand((unsigned)time(NULL));
	bool append;
    
    strcpy(DirName, argv[3]);
    printf("Output saved to directory: %s\n", DirName);
    
    int threadCount = atoi(argv[2]);
    
	if(threadCount > omp_get_max_threads()){threadCount = omp_get_max_threads();}
	if(threadCount <= 0)                   {threadCount = omp_get_max_threads();}
    omp_set_num_threads(threadCount);
	
    // read in parameters, or else use all defaults

	Inputs IniConditions;

	IniConditions = ReadParameters(argv[1]);
	if (argc==5 || argc==8)
		append = true;	// if a restart file is specified append new cell data to end of file

	else append = false;

	// Predefine arrays of cells for storing cell information
	Cell* old_cells = new Cell[maxCells];
	Cell* new_cells = new Cell[maxCells];
	
	printf("Created cell arrays \n");

	// Define a grid to store the cell positions
	BoxLength = (L_divide+2*cellRadius);	// physical size of each box
	//const int maxCellsPerBox = int(ceil(BoxLength*BoxLength/(PI*cellRadius*cellRadius)*2*1.5));  // max cells per box is based on the maximum density of cells
	const int maxCellsPerBox = int(ceil(BoxLength*BoxLength*BoxLength/(PI*cellRadius*cellRadius*L_divide*density_threshold/2)));
	UniformGrid Grid(BoxX, BoxY, BoxZ, maxCellsPerBox, BoxLength);

	//	printf("Created grid \n");
	printf("Created grid; max num of cells: %d \n", maxCellsPerBox);

	// initialize a configuration with N_cells cells in it
	int N_cells;
	if (argc==5 || argc==8)
	{
		N_cells = LoadCells(argv[4], old_cells, Grid, t0, initial_dt);
	}
	else N_cells = AddFirstCells(old_cells, L_divide, cellRadius, Grid, IniConditions);
	printf("Added initial cells \n");

	// files for output
	OutputFiles Files;
    CreateOutputFileLineage(0, Files, append);

	// Define an array to store the density of cells
	DoubleArray2D Height(BoxX*refinementGridHeight, BoxY*refinementGridHeight);	//	for storing the density of cells in each box
	DoubleArray3D Density(BoxX, BoxY, BoxZ);
	DoubleArray3D Density1(BoxX, BoxY, BoxZ);
	DoubleArray3D Density2(BoxX, BoxY, BoxZ);
	DoubleArray2D WallDensity(BoxX, BoxY);
	DoubleArray2D WallDensity1(BoxX, BoxY);
	DoubleArray2D WallDensity2(BoxX, BoxY);


	// Define double buffered environment arrays. Env array is initialized to zeros.
	EnvArray3D Environment(BoxX, BoxY, BoxZ);
	EnvArray3D oldEnvironment(BoxX, BoxY, BoxZ);
    AgaArray3D* FieldAgar[maxLevels];
    AgaArray3D* oldFieldAgar[maxLevels];
    AgaArray2D* FieldWall[maxLevels];
    AgaArray2D* oldFieldWall[maxLevels];
    for (int i=0; i<maxLevels; i++)
    {
        FieldAgar[i]= new AgaArray3D(BoxX, BoxY, BoxZAgar);
        oldFieldAgar[i]=new AgaArray3D(BoxX, BoxY, BoxZAgar);
        FieldWall[i]=new AgaArray2D(BoxX, BoxY);
        oldFieldWall[i]=new AgaArray2D(BoxX, BoxY);
    }
    
	printf("Created fields \n");

	CoordArray2D Normal(BoxX*refinementGridHeight, BoxY*refinementGridHeight);

	// Array for storing neighbours so that they don't need to be recalculated every time step
	int maxNeighbours = Grid.MaxCellsPerBox()*27+1;
	int** NeighbourList = InitializeNeighbourList(maxCells, maxNeighbours);

	// Start timing simulation
	time_t t_start, t_end;
	t_start = time(NULL);

	// Main simulation function
	RunSimulation(N_cells, old_cells, new_cells, NeighbourList, maxNeighbours, Grid, Files, append, Height, Density, Density1, Density2, WallDensity, WallDensity1, WallDensity2, Environment, oldEnvironment, FieldAgar, oldFieldAgar, FieldWall, oldFieldWall, Normal);
	
	// save end time
	t_end = time(NULL);

	// cleanup
    CloseOutputFileLineage(Files);
    
	delete[] old_cells;
	delete[] new_cells;
    
    for (int i=0; i<maxLevels; i++)
    {
        delete FieldAgar[i];
        delete oldFieldAgar[i];
        delete FieldWall[i];
        delete oldFieldWall[i];
    }
	
	printf("Total time for simulation %ld\n",t_end-t_start);
	return 0;
}
