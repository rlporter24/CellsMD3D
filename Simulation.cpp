#include "Simulation.h"
#include <stdio.h>
#include <math.h>

#include "Array.h"
#include "Cell.h"
#include "Cells.h"
#include "Compute.h"
#include "InputOutput.h"
#include "Constants.h"
#include "tools.h"
#include "UniformGrid.h"
#include "Nutrients.h"
#include "Neighbours.h"
#include "Forces.h"
#include "ClockIt.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

void SimSingleThread(int N_cells, Cell* old_cells, Cell* new_cells, int** NeighbourList, int maxNeighbours, UniformGrid& Grid, OutputFiles Files, bool append, DoubleArray2D& Height, DoubleArray3D& Density, DoubleArray3D& Density1, DoubleArray3D& Density2, DoubleArray2D& WallDensity,DoubleArray2D& WallDensity1,DoubleArray2D& WallDensity2,EnvArray3D& Environment, EnvArray3D& oldEnvironment, AgaArray3D** FieldAgar, AgaArray3D** oldFieldAgar, AgaArray2D** FieldWall, AgaArray2D** oldFieldWall, CoordArray2D& Normal);

int GetProcessorCount()
{
#ifdef _WIN32
    SYSTEM_INFO info;
    GetSystemInfo(&info);
    return info.dwNumberOfProcessors;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);    
#endif
}

void RunSimulation(int N_cells, Cell* old_cells, Cell* new_cells, int** NeighbourList, int maxNeighbours, UniformGrid& Grid, OutputFiles Files, bool append, DoubleArray2D& Height, DoubleArray3D& Density, DoubleArray3D& Density1, DoubleArray3D& Density2,DoubleArray2D& WallDensity,DoubleArray2D& WallDensity1,DoubleArray2D& WallDensity2, EnvArray3D& Environment, EnvArray3D& oldEnvironment, AgaArray3D** FieldAgar, AgaArray3D** oldFieldAgar, AgaArray2D** FieldWall, AgaArray2D** oldFieldWall, CoordArray2D& Normal)
{
    SimSingleThread(N_cells, old_cells, new_cells, NeighbourList, maxNeighbours, Grid, Files, append, Height, Density, Density1, Density2, WallDensity, WallDensity1, WallDensity2, Environment, oldEnvironment, FieldAgar, oldFieldAgar, FieldWall, oldFieldWall, Normal);
}

void SimSingleThread(int N_cells, Cell* old_cells, Cell* new_cells, int** NeighbourList, int maxNeighbours, UniformGrid& Grid, OutputFiles Files, bool append, DoubleArray2D& Height, DoubleArray3D& Density, DoubleArray3D& Density1, DoubleArray3D& Density2, DoubleArray2D& WallDensity,DoubleArray2D& WallDensity1,DoubleArray2D& WallDensity2,EnvArray3D& Environment, EnvArray3D& oldEnvironment, AgaArray3D** FieldAgar, AgaArray3D** oldFieldAgar, AgaArray2D** FieldWall, AgaArray2D** oldFieldWall, CoordArray2D& Normal)
{
    // initialize time (generations)
	double dt = initial_dt;	// time step	
	double t = t0;			// current time
	int minx, maxx, miny, maxy, maxz;

	// **********************Initialize*******************************

    // counter to determine when to output data
	double NextOutTime = OutputTime;
    double NextUpdateTime = UpdateTime;
    double HeightDensityNextUpdateTime = UpdateTime/10.0;
    double HeightDensityUpdateTime = UpdateTime/10.0;
	bool OutFlag = true;
    bool UpdateFlag = true;
    bool HeightDensityUpdateFlag = true;

	// filter for doing some smoothing and averaging of fields
	DoubleArray2D Filter(FilterLen,FilterLen);
	double value = 1.0/(double)(FilterLen*FilterLen);

	for (int ii = 0; ii<FilterLen; ii++)
	{
		for (int jj = 0; jj<FilterLen; jj++)
		{
			Filter.Set(ii,jj,value);
		}
	}
	int FilterDim = FilterLen/2;

	// make filter function for averaging the density
/*	DoubleArray3D Filter3D(3,3,3);
	double filter3D[27] = {	0.0010,    0.0081,    0.0010,
							0.0081,    0.0642,    0.0081,
							0.0010,    0.0081,    0.0010,

							0.0081,    0.0642,    0.0081,
							0.0642,    0.5102,    0.0642,
							0.0081,    0.0642,    0.0081,

							0.0010,    0.0081,    0.0010,
							0.0081,    0.0642,    0.0081,
							0.0010,    0.0081,    0.0010	};


	Filter3D.SetData(filter3D,27);*/
    
    DoubleArray3D Filter3D(FilterLen,FilterLen,FilterLen);
    double norm = 0;
    for (int ii=0; ii<FilterLen; ii++)
    {
        for (int jj=0; jj<FilterLen; jj++)
        {
            for (int kk=0; kk<FilterLen; kk++)
            {
     	        value = exp(-(ii-FilterLen/2)*(ii-FilterLen/2)-(jj-FilterLen/2)*(jj-FilterLen/2)-(kk-FilterLen/2)*(kk-FilterLen/2));
                norm = norm + value;
            }
        }
    }
    for (int ii=0; ii<FilterLen; ii++)
    {
        for (int jj=0; jj<FilterLen; jj++)
        {
            for (int kk=0; kk<FilterLen; kk++)
            {
 	        value = exp(-(ii-FilterLen/2)*(ii-FilterLen/2)-(jj-FilterLen/2)*(jj-FilterLen/2)-(kk-FilterLen/2)*(kk-FilterLen/2));
                Filter3D.Set(ii,jj,kk,value/norm);
            }
        }
    }


	// *******************Calculate important fields******************
	DoubleArray3D RoughDensity(Density.Size().x, Density.Size().y, Density.Size().z);
	DoubleArray3D RoughDensity1(Density1.Size().x, Density1.Size().y, Density1.Size().z);
	DoubleArray3D RoughDensity2(Density2.Size().x, Density2.Size().y, Density2.Size().z);
    DoubleArray3D insideColonyDen(Density.Size().x, Density.Size().y, Density.Size().z);
    
	DoubleArray2D RoughWallDensity(WallDensity.Size().x, WallDensity.Size().y);
	DoubleArray2D RoughWallDensity1(WallDensity1.Size().x, WallDensity1.Size().y);
    DoubleArray2D RoughWallDensity2(WallDensity2.Size().x, WallDensity2.Size().y);
	DoubleArray3D RoughDensityShiftP(Density.Size().x, Density.Size().y, Density.Size().z);
	DoubleArray3D RoughDensity1ShiftP(Density1.Size().x, Density1.Size().y, Density1.Size().z);
	DoubleArray3D RoughDensity2ShiftP(Density2.Size().x, Density2.Size().y, Density2.Size().z);
	DoubleArray2D RoughWallDensityShiftP(WallDensity.Size().x, WallDensity.Size().y);
	DoubleArray2D RoughWallDensity1ShiftP(WallDensity1.Size().x, WallDensity1.Size().y);
	DoubleArray2D RoughWallDensity2ShiftP(WallDensity2.Size().x, WallDensity2.Size().y);
	DoubleArray3D DensityShiftP(Density.Size().x, Density.Size().y, Density.Size().z);
	DoubleArray3D Density1ShiftP(Density1.Size().x, Density1.Size().y, Density1.Size().z);
	DoubleArray3D Density2ShiftP(Density2.Size().x, Density2.Size().y, Density2.Size().z);
	DoubleArray2D WallDensityShiftP(WallDensity.Size().x, WallDensity.Size().y);
	DoubleArray2D WallDensity1ShiftP(WallDensity1.Size().x, WallDensity1.Size().y);
	DoubleArray2D WallDensity2ShiftP(WallDensity2.Size().x, WallDensity2.Size().y);
	DoubleArray2D RoughHeight(Height.Size().x, Height.Size().y);
	RoughHeight.Initialize(cellRadius);
	RoughDensity.Initialize(0.0);
	RoughDensity1.Initialize(0.0);
	RoughDensity2.Initialize(0.0);
    insideColonyDen.Initialize(0.0);
    
	RoughWallDensity.Initialize(0.0);
	RoughWallDensity1.Initialize(0.0);
	RoughWallDensity2.Initialize(0.0);
	RoughDensityShiftP.Initialize(0.0);
	RoughDensity1ShiftP.Initialize(0.0);
	RoughDensity2ShiftP.Initialize(0.0);
	RoughWallDensityShiftP.Initialize(0.0);
	RoughWallDensity1ShiftP.Initialize(0.0);
	RoughWallDensity2ShiftP.Initialize(0.0);
	DensityShiftP.Initialize(0.0);
	Density1ShiftP.Initialize(0.0);
	Density2ShiftP.Initialize(0.0);
	WallDensityShiftP.Initialize(0.0);
	WallDensity1ShiftP.Initialize(0.0);
	WallDensity2ShiftP.Initialize(0.0);

	// Make rough wall
	DoubleArray2D Wall(BoxX,BoxY);
	for (int ii = 0; ii < BoxX; ii++)
	{
		for (int jj = 0; jj < BoxY; jj++)
		{
			Wall.Set(ii,jj,((float)rand()/RAND_MAX-0.5)*wall_rough);
		}
	}
	printf("Created wall \n");

	int* IDlist = new int[maxNeighbours];
	int IDlen = 0;

    // cell stress tensor
	Tensor stressTensor;

	// net force
	DoubleCoord Fnet;

    // dividing cells
	int* dividingCells = new int[maxCells];		// list of cells that need to divide at the end of a timestep
	int numDivide;

	
	// calculate fields, output to file
    CreateOutputFiles(0, Files, append);
	GetDensity(RoughDensity, RoughDensity1, RoughDensity2, insideColonyDen, Height, Grid, old_cells, minx, maxx, miny, maxy, maxz);
	ShiftDensity(RoughDensity, RoughDensity1, RoughDensity2, RoughWallDensity, RoughWallDensity1, RoughWallDensity2,RoughDensityShiftP, RoughDensity1ShiftP, RoughDensity2ShiftP,RoughWallDensityShiftP, RoughWallDensity1ShiftP, RoughWallDensity2ShiftP,BoxX, BoxY, BoxZ);
	RoughDensityShiftP.Output(Files.roughDensity,3);
	RoughDensity1ShiftP.Output(Files.roughDensity1,3);
	RoughDensity2ShiftP.Output(Files.roughDensity2,3);

	BoxAverage(RoughDensity, Density, RoughWallDensity, WallDensity, Filter3D, int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
	BoxAverage(RoughDensity1, Density1,RoughWallDensity1, WallDensity1, Filter3D, int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
	BoxAverage(RoughDensity2, Density2, RoughWallDensity2, WallDensity2,Filter3D, int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
	BoxAverage(RoughDensityShiftP, DensityShiftP, RoughWallDensityShiftP, WallDensityShiftP, Filter3D, int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
	BoxAverage(RoughDensity1ShiftP, Density1ShiftP,RoughWallDensity1ShiftP, WallDensity1ShiftP, Filter3D, int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
	BoxAverage(RoughDensity2ShiftP, Density2ShiftP, RoughWallDensity2ShiftP, WallDensity2ShiftP,Filter3D, int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
	DensityShiftP.Output(Files.density,3);
	Density1ShiftP.Output(Files.density1,3);
	Density2ShiftP.Output(Files.density2,3);
	WallDensityShiftP.Output(Files.walldensity);
	WallDensity1ShiftP.Output(Files.walldensity1);
    WallDensity2ShiftP.Output(Files.walldensity2);
	Height_Average(Height, RoughHeight, minx, maxx, miny, maxy);
	RoughHeight.Output(Files.roughheight);
	Smooth(RoughHeight, Height, Filter, FilterDim, minx, maxx, miny, maxy);
	Height.Output(Files.height);


	// find surface normal
	GetSurfaceNormal(Normal, Height, minx, maxx, miny, maxy);
	Normal.Output(Files.normal);


	int Nconv = 0;

	Nconv = UpdateEnvArray(&Environment, &oldEnvironment, FieldAgar, oldFieldAgar, FieldWall, oldFieldWall, DensityShiftP, Density1ShiftP, Density2ShiftP, WallDensityShiftP, WallDensity1ShiftP, WallDensity2ShiftP, minx, maxx, miny, maxy, maxz, Height, insideColonyDen);
    
	//ShiftGrowthrate(&Environment, &FieldWall,BoxX, BoxY, BoxZ);

	Environment.Output(Files.env,1);
    for (int level=0;level<maxLevels;level++)
    {
        FieldAgar[level]->Append(Files.aga,1);
        FieldWall[level]->Append(Files.wal);
    }
    

	for (int cellID=0;cellID<N_cells;cellID++)
	{
		// output cell data
		if (OutFlag)
		{

			IntCoord XYAddress = Grid.GetXY( Grid.GetAddress(average(old_cells[cellID].Position)) );
			DoubleCoord T(0,0,0);
			F_surf_tension(old_cells[cellID], Grid, XYAddress, Height, Normal, Fnet, T);
			Output(Files.cells, cellID, t, old_cells[cellID], Fnet);
			//printf("%6f %6d %6d\n", t, N_cells, cellID);
		}

	}
	fflush(Files.cells);
	OutFlag = false;
    CloseOutputFiles(Files);

	// **************************** Loop through time and evolve colony ****************************
    
    ClockIt     t0, t1, t2, t3, t4, t00;
    double      s0, s1, s2, s3, s4, st, s00;
    s0=0; s1=0; s2=0; s3=0; s4=0; st=0; s00=0;
    int koutput=0;
    
    rewind(Files.lineage);
    while (t<=t_max)
	{
		numDivide = 0;					// number of cells dividing

		// update fields only when UpdateFlag is true
        t0.start();
        
        if (UpdateFlag)
        {
            UpdateFlag = false;
            
            // find density
            GetDensity(RoughDensity, RoughDensity1, RoughDensity2, insideColonyDen, Height, Grid, old_cells, minx, maxx, miny, maxy, maxz);
            ShiftDensity(RoughDensity, RoughDensity1, RoughDensity2, RoughWallDensity, RoughWallDensity1, RoughWallDensity2,RoughDensityShiftP, RoughDensity1ShiftP, RoughDensity2ShiftP,RoughWallDensityShiftP, RoughWallDensity1ShiftP, RoughWallDensity2ShiftP,BoxX, BoxY, BoxZ);
            BoxAverage(RoughDensity, Density, RoughWallDensity, WallDensity, Filter3D, int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
            BoxAverage(RoughDensity1, Density1,RoughWallDensity1, WallDensity1, Filter3D, int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
            BoxAverage(RoughDensity2, Density2, RoughWallDensity2, WallDensity2,Filter3D, int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
            BoxAverage(RoughDensityShiftP, DensityShiftP, RoughWallDensityShiftP, WallDensityShiftP, Filter3D, int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
            BoxAverage(RoughDensity1ShiftP, Density1ShiftP,RoughWallDensity1ShiftP, WallDensity1ShiftP, Filter3D, int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
            BoxAverage(RoughDensity2ShiftP, Density2ShiftP, RoughWallDensity2ShiftP, WallDensity2ShiftP,Filter3D, int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
            Height_Average(Height, RoughHeight, minx, maxx, miny, maxy);
            //			RoughHeight.Output(Files.height);
            Smooth(RoughHeight, Height, Filter, FilterDim, minx, maxx, miny, maxy);
            //			Height.Output(Files.height);
            
            // find surface normal
            GetSurfaceNormal(Normal, Height, minx, maxx, miny, maxy);
            
            // find nutrient concentrations
            Nconv = UpdateEnvArray(&Environment, &oldEnvironment, FieldAgar, oldFieldAgar, FieldWall, oldFieldWall, DensityShiftP, Density1ShiftP, Density2ShiftP, WallDensityShiftP, WallDensity1ShiftP, WallDensity2ShiftP, minx, maxx, miny, maxy, maxz, Height, insideColonyDen);
            //ShiftGrowthrate(&Environment, &FieldWall,BoxX, BoxY, BoxZ);
            if (Nconv>=maxIter)
                printf("Environment has not converged");
//            MyAssert(Nconv<maxIter,"Environment has not converged");
            
            // update neighbour lists
            getNeighbours(old_cells, N_cells, Grid, NeighbourList, maxNeighbours);
            
        }
	s0+=t0.stop()/1000.0;
        t00.start();
        if (HeightDensityUpdateFlag)
        {
            HeightDensityUpdateFlag = false;
	    
            // find height
            GetHeight(RoughDensity, Height, Grid, old_cells, minx, maxx, miny, maxy, maxz);
            Height_Average(Height, RoughHeight, minx, maxx, miny, maxy);
            Smooth(RoughHeight, Height, Filter, FilterDim, minx, maxx, miny, maxy);
            // find surface normal
            GetSurfaceNormal(Normal, Height, minx, maxx, miny, maxy);
            
            // update neighbour lists
	    //            getNeighbours(old_cells, N_cells, Grid, NeighbourList, maxNeighbours);
        }
       	s00+=t00.stop()/1000.0;
        t1.start();
        
		// calculate forces and move cells
        int cellID;
        if (OutFlag)
        {
            //CreateOutputFiles(int(t*10), Files, append);
            CreateOutputFiles(koutput++, Files, append);
            for (cellID=0;cellID<N_cells;cellID++)
            {
            // output cell data
            
                mean_stress(old_cells[cellID], old_cells, NeighbourList[cellID], Grid, Wall, Height, Normal, stressTensor, Fnet);
                Output(Files.cells, cellID, t, old_cells[cellID], stressTensor);
            }
            
        }
        
#pragma omp parallel for default(shared) private(cellID) schedule(static)
        for (cellID=0;cellID<N_cells;cellID++)
        {
            MoveCell(cellID, Grid, old_cells, new_cells, NeighbourList[cellID], dt, Height, Normal, Wall);
            GrowCell(new_cells[cellID], cellID, dt, dividingCells, numDivide, Environment, FieldWall, Grid);
        }
        
	//	fflush(Files.cells);
        s1+=t1.stop()/1000.0;
        t2.start();
        
		// output fields to file
		if (OutFlag)
		{
			OutFlag = false;

			// save restart file with cell information
			SaveCells(Files.restart, new_cells, N_cells, t, dt);

			// save field data
            RoughHeight.Output(Files.roughheight);
			Height.Output(Files.height);
			RoughDensityShiftP.Output(Files.roughDensity,3);
			RoughDensity1ShiftP.Output(Files.roughDensity1,3);
			RoughDensity2ShiftP.Output(Files.roughDensity2,3);
			DensityShiftP.Output(Files.density,3);
			Density1ShiftP.Output(Files.density1,3);
			Density2ShiftP.Output(Files.density2,3);
			WallDensityShiftP.Output(Files.walldensity);
			WallDensity1ShiftP.Output(Files.walldensity1);
			WallDensity2ShiftP.Output(Files.walldensity2);
			Normal.Output(Files.normal);
			Environment.Output(Files.env,1);
            for (int level=0;level<maxLevels;level++)
            {
                FieldAgar[level]->Append(Files.aga,1);
                FieldWall[level]->Append(Files.wal);
            }
            printf("t = %6f, %5d cells;      ", t, N_cells);
            st += s0+s1+s2+s3+s4+s00;
            printf("sim time = %3d:%2d:%2d;      ", int(st/60.0/60.0), int(st/60.0)%60, int(st)%60);
            printf("time cost = %g %g %g %g %g %g. \n", s0, s00, s1, s2, s3, s4);
            s0=0; s1=0; s2=0; s3=0; s4=0; s00=0;
			fflush(stdout);
            CloseOutputFiles(Files);
		}
        s2+=t2.stop()/1000.0;
        t3.start();

		// Division (must be done after integration step because new neighbours are created)
		for (int cellCount = 0; cellCount < numDivide; cellCount++)
		{
			if (N_cells>(maxCells-1))
			{
				break;
				printf("Too many cells!");
			}
		
			// ID of cell that is undergoing division
			const int cellID = dividingCells[cellCount];

			DivideCell(cellID, N_cells, new_cells, Grid, NeighbourList[cellID], Wall, Height, Normal, t);
            fprintf(Files.lineage, "%d %d\n", cellID, N_cells);
			N_cells++;

			// make list of all neighbours, new cell, and old cell
			IDlen = NeighbourList[cellID][0];
			memcpy(IDlist, NeighbourList[cellID], (IDlen+1)*sizeof(int)); // copy neighbours
			IDlist[0] = cellID;
			IDlist[IDlen+1] = N_cells-1;
			IDlen+=2;

			// update the neighbours of all of these cells
			getNeighbours(new_cells, Grid, NeighbourList, maxNeighbours, IDlist, IDlen);
		}
        s3+=t3.stop()/1000.0;
        t4.start();

		// switch positions of old and new cells
		{
			Cell* temp = new_cells;
			new_cells = old_cells;
			old_cells = temp;
		}

		t += dt;
		
		// determine if we're writing output next time
		OutFlag = (NextOutTime<=0);
		NextOutTime = (OutFlag? OutputTime: NextOutTime) - dt;

		// determine if we're updating fields
		UpdateFlag = (NextUpdateTime<=0);
		NextUpdateTime = (UpdateFlag? UpdateTime: NextUpdateTime) - dt;
        
        HeightDensityUpdateFlag = (HeightDensityNextUpdateTime<=0);
        HeightDensityNextUpdateTime = (HeightDensityUpdateFlag? HeightDensityUpdateTime: HeightDensityNextUpdateTime) - dt;

        s4+=t4.stop()/1000.0;
        

	}
    fflush(Files.lineage);
    printf("Done\n");
}

