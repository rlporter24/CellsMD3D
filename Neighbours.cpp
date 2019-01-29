#include <math.h>
#include <float.h>

#include "Neighbours.h"
#include "Cell.h"
#include "Constants.h"
#include "tools.h"
#include "UniformGrid.h"

/**********************
 * The functions in this file create a list of neighboring cells to speed up calculation of cell-cell forces
 *********************/
int** InitializeNeighbourList(int maxCells, int maxNeighbours)
{
	int** NeighbourList = new int*[maxCells];
	for (int icell = 0; icell<maxCells; icell++)
	{
		NeighbourList[icell] = new int[maxNeighbours];
		NeighbourList[icell][0] = 0;
	}
	return NeighbourList;
}

void getNeighbours(const Cell* cell_array, int N_cells, UniformGrid& Grid, int** NeighbourList, int maxNeighbours)
{
	UniformGrid::Address Address;
	int numNeighbours;

	// loop through all cells
	for (int cellID=0;cellID<N_cells;cellID++)
	{
		// get coarse-grained neighbours from Grid
		Address = Grid.GetAddress(average(cell_array[cellID].Position));
		numNeighbours = Grid.GetNeighbours(cellID, Address, NeighbourList[cellID], maxNeighbours);
		numNeighbours = reduceNeighbours(cell_array[cellID], cell_array, NeighbourList[cellID]);
	}
}

void getNeighbours(const Cell* cell_array, UniformGrid& Grid, int** NeighbourList, int maxNeighbours, int* ID, int ID_length)
{
	UniformGrid::Address Address;
	int numNeighbours;
	int cellID;

	// loop through all cells
	for (int ii=0;ii<ID_length;ii++)
	{
		cellID = ID[ii];

		// get coarse-grained neighbours from Grid
		Address = Grid.GetAddress(average(cell_array[cellID].Position));
		numNeighbours = Grid.GetNeighbours(cellID, Address, NeighbourList[cellID], maxNeighbours);
		numNeighbours = reduceNeighbours(cell_array[cellID], cell_array, NeighbourList[cellID]);
	}
}

int reduceNeighbours(const Cell& cell, const Cell* cell_array, int* neighbours)
{
	int ID;
	double d;
	DoubleCoord c1, c2, v;
	int numNeighReduced = 1;
	double dp, dq;

	// position of cell center of mass
	DoubleCoord cm = average(cell.Position);
	int numNeighbours = neighbours[0];

	// loop through neighbours and find the forces
	for (int neighbourID = 1; neighbourID<numNeighbours+1; neighbourID++)
	{
		// quick check if cell is in right ballpark
		ID = neighbours[neighbourID];	// the ID of the current neighbour
		v = diff(cell_array[ID].Position.p,cm);
		dp = sqrt(dot(v,v));	// distance squared between p and cm
		v = diff(cell_array[ID].Position.q,cm);
		dq = sqrt(dot(v,v));	// distance squared between p and cm

		if ((dp<(cell.Length/2.0+Rc*cell.Radius))||(dq<(cell.Length/2.0+Rc*cell.Radius)))	// check to see if it is possible that the cells are overlapping
		{
			min_distance(cell,cell_array[ID],d,c1,c2);	// get actual minimum distance between cells
			if (d<Rc*cell.Radius)
			{
				// reassign neighbour list based on which neighbours are within a distance r_th of the current cell
				// will reduce the number of neighbours for future sum forces calls, which should speed things up considerably
				neighbours[numNeighReduced++] = ID;
			}
		}
	}
	neighbours[0] = numNeighReduced-1;
	return numNeighReduced-1;

}

double clamp(double n, double minn, double maxn)
{
	if (n<minn)
		return minn;
	else if (n>maxn)
		return maxn;
	else
		return n;

}

// find the minimum distance between any two cells
void min_distance(const Cell& cell1, const Cell& cell2, double& d, DoubleCoord& c1, DoubleCoord& c2)
{
	DoubleCoord p1,q1,p2,q2,v1,v2,r,cv;
	double l1,l2,f,c,b,denom,tnom,s,t;

    // S1 and S2 are line Segments defined by the points (p1,q1) and (p2,q2)
    // respectively
	p1 = cell1.Position.p;
	q1 = cell1.Position.q;
	p2 = cell2.Position.p;
	q2 = cell2.Position.q;

	// find the direction Segment of S1 and S2
	v1 = diff(q1,p1);
	v2 = diff(q2,p2);
	r = diff(p1,p2);	// vec between two start points

	l1 = dot(v1,v1);
	l2 = dot(v2,v2);
	f = dot(v2,r);
	c = dot(v1,r);
	b = dot(v1,v2);

	denom = l1*l2-b*b;

	if (denom!=0) // if not parallel Segments
		s = clamp((b*f-c*l2)/denom,0,1);	// closest point on cell1 to cell2 within the Segment cell1
	else
		s = 0;	// can choose arbitrary number if parallel

	// compute point on L2 closest to S1
	tnom = b*s+f;

    // check if this is within the Segment (t between 0 and 1)
	if (tnom<0)
	{
		t = 0;
		s = clamp(-c/l1,0,1);
	}
	else if (tnom>l2)
	{
	    t = 1;
	    s = clamp((b-c)/l1,0,1);
	}
	else
	{
	    t = tnom/l2;
	}
	c1 = sum(p1,scale(v1,s));
	c2 = sum(p2,scale(v2,t));
	cv = diff(c1,c2);
	d = sqrt(dot(cv,cv));		// length of segment connecting two cells

}



