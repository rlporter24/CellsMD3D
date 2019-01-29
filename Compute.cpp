#include <float.h>
#include <math.h>

#include "Array.h"
#include "Cell.h"
#include "Compute.h"
#include "Constants.h"
#include "Forces.h"
#include "tools.h"
#include "UniformGrid.h"

/*************************************************
 * Computes important fields from averaged cell quantities
 *************************************************/

// Computes the mean stress tensor on a cell, decomposed into principal axes
// remember that the forces are decomposed into forces acting on two spherical caps of cells
void mean_stress(const Cell& cell, const Cell* cell_array, const int* neighbours, UniformGrid& Grid, DoubleArray2D& Wall, DoubleArray2D& Height, CoordArray2D& Normal, Tensor& stressTensor, DoubleCoord& Fnet)
{
	int ID;
	DoubleCoord F, r, F2, r2, dynFric, staFric; // gives force, position of force

	double d;

	double V = cell.Length*PI*cell.Radius*cell.Radius + 4.0/3.0*PI*cell.Radius*cell.Radius*cell.Radius;

	IntCoord XYAddress = Grid.GetXY( Grid.GetAddress(average(cell.Position)));

	int numNeighbours = neighbours[0];

	// initialize stress tensor
	stressTensor = Tensor(0,0,0,0,0);

	// initialize net force
	Fnet = DoubleCoord(0,0,0);

	// find center of mass of the cell
	DoubleCoord CM = average(cell.Position);

	// loop through neighbours and find the forces
	for (int neighbourID = 1; neighbourID<numNeighbours+1; neighbourID++){
		ID = neighbours[neighbourID];	// the ID of the current neighbour
		F_cc(cell, cell_array[ID], F, r, d);
		Fnet = sum(Fnet, F);	// net force is the sum of all forces

		// position of contact force relative to center of mass
		r = diff(r,CM);

		// calc stress tensor
		stressTensor.xx -= F.x*r.x;
		stressTensor.yy -= F.y*r.y;
		stressTensor.zz -= F.z*r.z;

	}

	// calculate cell-wall forces
	if (min(cell.Position.q.z,cell.Position.p.z)<1.2*cell.Radius){
		double wall_y = Wall.Get(XYAddress.x, XYAddress.y);
		F_cw(cell, wall_y, F, F2, r, r2, dynFric, staFric);
		Fnet = sum(Fnet, F);	// net force is the sum of all forces
		// position of contact force relative to center of mass
		r = diff(r,CM);

		// calc stress tensor
		stressTensor.xx -= F.x*r.x;
		stressTensor.yy -= F.y*r.y;
		stressTensor.zz -= F.z*r.z;

		Fnet = sum(Fnet, F2);	// net force is the sum of all forces
		// position of contact force relative to center of mass
		r2 = diff(r2,CM);

		// calc stress tensor
		stressTensor.xx -= F2.x*r2.x;
		stressTensor.yy -= F2.y*r2.y;
		stressTensor.zz -= F2.z*r2.z;
	}


	// surface tension
	DoubleCoord T(0,0,0);
	F_surf_tension(cell, Grid, XYAddress, Height, Normal, F, T);
	Fnet = sum(Fnet, F);	// net force is the sum of all forces

	// position of forces are at nodes
	r = diff(cell.Position.p,CM);
	r.z += cell.Radius;

	r2 = diff(cell.Position.q,CM);
	r2.z += cell.Radius;


	// calc stress tensor
	stressTensor.xx -= F.x*r.x*0.5 - F.x*r2.x*0.5;
	stressTensor.yy -= F.y*r.y*0.5 - F.y*r2.y*0.5;
	stressTensor.zz -= F.z*r.z*0.5 - F.z*r2.z*0.5;

	// normalize by the volume
	stressTensor.xx /= V;
	stressTensor.yy /= V;
	stressTensor.zz /= V;

	// find radial and azimuthal components of stress tensor
	double Lcm = sqrt(CM.x*CM.x+CM.y*CM.y);
	if (Lcm>DBL_EPSILON)
	{
		DoubleCoord uv = scale(CM, 1.0/Lcm);	// unit vector in plane to center of mass of cell
		stressTensor.rr = stressTensor.xx*uv.x + stressTensor.yy*uv.y;
		stressTensor.tt = -stressTensor.xx*uv.y + stressTensor.yy*uv.x;
	}
	else
	{
		stressTensor.rr = (stressTensor.xx+stressTensor.yy);
		stressTensor.tt = 0.0;
	}

}

// compute height and cell density fields from positions of individual cells
void GetDensity(DoubleArray3D& DensityVec, DoubleArray3D& Density1Vec, DoubleArray3D& Density2Vec, DoubleArray3D& insideColonyDen, DoubleArray2D& Height, UniformGrid& Grid, Cell* cells, int& minx, int& maxx, int& miny, int& maxy, int& maxz)
{
	double density, density1, density2, height[refinementGridHeight][refinementGridHeight];
	int count;
	int* cellID;
	minx = DensityVec.Size().x;
	maxx = 0;
	miny = DensityVec.Size().y;
	maxy = 0;
	maxz = 0;
	Cell cell;
    DoubleCoord cm;
    int i_ht, j_ht;
	double top;
    double Vbox = BoxLength*BoxLength*BoxLength;

	for (int boxx = 0; boxx < DensityVec.Size().x; boxx++)
	{
		for (int boxy = 0; boxy < DensityVec.Size().y; boxy++)
		{
            for (i_ht = 0; i_ht < refinementGridHeight; i_ht++)
            {
                for (j_ht = 0; j_ht < refinementGridHeight; j_ht++)
                {
                    height[i_ht][j_ht] = cellRadius;
                }
            }
			for (int boxz = 0; boxz < DensityVec.Size().z; boxz++)
			{
				// number of cells per box
				count = Grid.GetCells(boxx, boxy, boxz, cellID);

				// find a bounding box where the cells exist
				if (count>0)
				{
					minx = min(minx,boxx);
					maxx = max(maxx,boxx);
					miny = min(miny,boxy);
					maxy = max(maxy,boxy);
					maxz = max(maxz,boxz);
				}

				density = 0;
				density1 = 0;
				density2 = 0;
				for (int celli = 0; celli < count; celli++)
				{
					// add up lengths of all cells to get an approximation of the density
					cell = cells[cellID[celli]];
                    cm = average(cell.Position);
                    i_ht = int((cm.x/BoxLength-(boxx-BoxX/2))*refinementGridHeight);
                    j_ht = int((cm.y/BoxLength-(boxy-BoxY/2))*refinementGridHeight);
                    if (i_ht<0)
                        i_ht=0;
                    if (i_ht>refinementGridHeight-1)
                        i_ht=refinementGridHeight-1;
                    if (j_ht<0)
                        j_ht=0;
                    if (j_ht>refinementGridHeight-1)
                        j_ht=refinementGridHeight-1;
                    top = max(cell.Position.p.z,cell.Position.q.z);
					height[i_ht][j_ht] = max(height[i_ht][j_ht],top);
                    //height = max(height,top);
					switch(cell.Type)
					{
						case 1:
						{
							//density1 += cell.Length + 2.0*cell.Radius;
							density1 += cell.Length*PI*cell.Radius*cell.Radius;
							//printf("type 1 added\n");
							break;
						}
						case 2:
						{
							//density2 += cell.Length + 2.0*cell.Radius;
							density2 += cell.Length*PI*cell.Radius*cell.Radius;
							//printf("type 2 added\n");
							break;
						}
						default:
							printf("CellType error!\n");
					}
					//density += cell.Length + 2.0*cell.Radius;
					density += cell.Length*PI*cell.Radius*cell.Radius;

				}
				if (density==0)
				{
					Density1Vec.Set(boxx,boxy,boxz,0);
					Density2Vec.Set(boxx,boxy,boxz,0);
					DensityVec.Set(boxx,boxy, boxz, 0);//min(1.0,density/density_threshold));	// density is normalized to the close packing fraction
                    insideColonyDen.Set(boxx,boxy, boxz, 0);

				}
				else
				{
					Density1Vec.Set(boxx,boxy,boxz,density1/density);
					Density2Vec.Set(boxx,boxy,boxz,density2/density);
					DensityVec.Set(boxx,boxy, boxz, density/Vbox);//min(1.0,density/density_threshold));	// density is normalized to the close packing fraction
                    insideColonyDen.Set(boxx,boxy, boxz, 1.0);
				}
                for (i_ht = 0; i_ht < refinementGridHeight; i_ht++)
                {
                    for (j_ht = 0; j_ht < refinementGridHeight; j_ht++)
                    {
                        Height.Set(boxx*refinementGridHeight+i_ht,boxy*refinementGridHeight+j_ht,height[i_ht][j_ht]);
                    }
                }
//				Height.Set(boxx,boxy,height);
			}
		}
	}
    for (int boxx = 0; boxx < DensityVec.Size().x; boxx++)
    {
        for (int boxy = 0; boxy < DensityVec.Size().y; boxy++)
        {
            for (int boxz = 0; boxz < DensityVec.Size().z; boxz++)
            {
                if (insideColonyDen.Get(boxx,boxy,boxz)<0.5)
                {
                    if (boxz>0 && boxx>0 && boxy>0 && boxx<DensityVec.Size().x-1 && boxy<DensityVec.Size().y-1 && boxz < DensityVec.Size().z-1)
                    {
                        if ((insideColonyDen.Get(boxx+1,boxy,boxz)>0.5) || (insideColonyDen.Get(boxx-1,boxy,boxz)>0.5) || (insideColonyDen.Get(boxx,boxy+1,boxz)>0.5) || (insideColonyDen.Get(boxx,boxy-1,boxz)>0.5) || (insideColonyDen.Get(boxx,boxy,boxz+1)>0.5) || (insideColonyDen.Get(boxx,boxy,boxz-1)>0.5))
                            insideColonyDen.Set(boxx,boxy, boxz, 0.4);
                    }
                    else if (boxz==0 && boxx>0 && boxy>0 && boxx<DensityVec.Size().x-1 && boxy<DensityVec.Size().y-1)
                    {
                        if ((insideColonyDen.Get(boxx+1,boxy,boxz)>0.5) || (insideColonyDen.Get(boxx-1,boxy,boxz)>0.5) || (insideColonyDen.Get(boxx,boxy+1,boxz)>0.5) || (insideColonyDen.Get(boxx,boxy-1,boxz)>0.5) || (insideColonyDen.Get(boxx,boxy,boxz+1)>0.5) )
                            insideColonyDen.Set(boxx,boxy, boxz, 0.4);
                    }
                }
            }
        }
    }
	minx = minx - 5;
	maxx = maxx + 5;
	miny = miny - 5;
	maxy = maxy + 5;
	maxz = maxz + 5;
}

// compute height field from positions of individual cells
void GetHeight(DoubleArray3D& DensityVec, DoubleArray2D& Height, UniformGrid& Grid, Cell* cells, int& minx, int& maxx, int& miny, int& maxy, int& maxz)
{
    double density, density1, density2;
    int count;
    int* cellID;
    minx = DensityVec.Size().x;
    maxx = 0;
    miny = DensityVec.Size().y;
    maxy = 0;
    maxz = 0;
    Cell cell;
    double top;
    DoubleCoord cm;
    double height[refinementGridHeight][refinementGridHeight];
    double Vbox = BoxLength*BoxLength*BoxLength;
    int i_ht, j_ht;
    int boxes_pos_ct;
    
    for (int boxx = 0; boxx < DensityVec.Size().x; boxx++)
    {
        for (int boxy = 0; boxy < DensityVec.Size().y; boxy++)
        {
            for (i_ht = 0; i_ht < refinementGridHeight; i_ht++)
            {
                for (j_ht = 0; j_ht < refinementGridHeight; j_ht++)
                {
                    height[i_ht][j_ht] = cellRadius;
                }
            }
            boxes_pos_ct = 0;
//            height = cellRadius;
            for (int boxz = DensityVec.Size().z-1; boxz >= 0; boxz--)
            {
                // number of cells per box
                count = Grid.GetCells(boxx, boxy, boxz, cellID);
                
                // find a bounding box where the cells exist
                if (count>0)
                {
                    minx = min(minx,boxx);
                    maxx = max(maxx,boxx);
                    miny = min(miny,boxy);
                    maxy = max(maxy,boxy);
                    maxz = max(maxz,boxz);
                    boxes_pos_ct++;
                }
                
                for (int celli = 0; celli < count; celli++)
                {
                    // add up lengths of all cells to get an approximation of the density
                    cell = cells[cellID[celli]];
                    cm = average(cell.Position);
                    i_ht = int((cm.x/BoxLength-(boxx-BoxX/2))*refinementGridHeight);
                    j_ht = int((cm.y/BoxLength-(boxy-BoxY/2))*refinementGridHeight);
                    if (i_ht<0)
                        i_ht=0;
                    if (i_ht>refinementGridHeight-1)
                        i_ht=refinementGridHeight-1;
                    if (j_ht<0)
                        j_ht=0;
                    if (j_ht>refinementGridHeight-1)
                        j_ht=refinementGridHeight-1;
                    top = max(cell.Position.p.z,cell.Position.q.z);
                    height[i_ht][j_ht] = max(height[i_ht][j_ht],top);
                    
                }

                for (i_ht = 0; i_ht < refinementGridHeight; i_ht++)
                {
                    for (j_ht = 0; j_ht < refinementGridHeight; j_ht++)
                    {
                        Height.Set(boxx*refinementGridHeight+i_ht,boxy*refinementGridHeight+j_ht,height[i_ht][j_ht]);
                    }
                }
//                Height.Set(boxx,boxy,height);

                if (boxes_pos_ct>3)
                {
                    break;
                }
                
            }
        }
    }
    minx = minx - 5;
    maxx = maxx + 5;
    miny = miny - 5;
    maxy = maxy + 5;
    maxz = maxz + 5;
}


// shift the density profile from box center to grid point
void ShiftDensity(DoubleArray3D& DensityVec, DoubleArray3D& Density1Vec, DoubleArray3D& Density2Vec,DoubleArray2D& DensityWallVec, DoubleArray2D& DensityWallVec1, DoubleArray2D& DensityWallVec2, DoubleArray3D& DensityVecShiftP, DoubleArray3D& Density1VecShiftP, DoubleArray3D& Density2VecShiftP,DoubleArray2D& DensityWallVecShiftP, DoubleArray2D& DensityWallVec1ShiftP, DoubleArray2D& DensityWallVec2ShiftP,int BoxX, int BoxY, int BoxZ)
{
	double densityvectemp, densityvectemp1, densityvectemp2;

	for (int boxx = 0; boxx < BoxX-1; boxx++)
	{
		for (int boxy = 0; boxy < BoxY-1; boxy++)
		{
			for (int boxz = 0; boxz < BoxZ-1; boxz++)
			{
				densityvectemp=(DensityVec.Get(boxx,boxy,boxz)+DensityVec.Get(boxx+1,boxy,boxz)+DensityVec.Get(boxx,boxy+1,boxz)+DensityVec.Get(boxx+1,boxy+1,boxz)+DensityVec.Get(boxx,boxy,boxz+1)+DensityVec.Get(boxx+1,boxy,boxz+1)+DensityVec.Get(boxx,boxy+1,boxz+1)+DensityVec.Get(boxx+1,boxy+1,boxz+1))/8.0;
				densityvectemp1=(Density1Vec.Get(boxx,boxy,boxz)+Density1Vec.Get(boxx+1,boxy,boxz)+Density1Vec.Get(boxx,boxy+1,boxz)+Density1Vec.Get(boxx+1,boxy+1,boxz)+Density1Vec.Get(boxx,boxy,boxz+1)+Density1Vec.Get(boxx+1,boxy,boxz+1)+Density1Vec.Get(boxx,boxy+1,boxz+1)+Density1Vec.Get(boxx+1,boxy+1,boxz+1))/8.0;
				densityvectemp2=(Density2Vec.Get(boxx,boxy,boxz)+Density2Vec.Get(boxx+1,boxy,boxz)+Density2Vec.Get(boxx,boxy+1,boxz)+Density2Vec.Get(boxx+1,boxy+1,boxz)+Density2Vec.Get(boxx,boxy,boxz+1)+Density2Vec.Get(boxx+1,boxy,boxz+1)+Density2Vec.Get(boxx,boxy+1,boxz+1)+Density2Vec.Get(boxx+1,boxy+1,boxz+1))/8.0;
				DensityVecShiftP.Set(boxx,boxy,boxz,densityvectemp);
				Density1VecShiftP.Set(boxx,boxy,boxz,densityvectemp1);
				Density2VecShiftP.Set(boxx,boxy,boxz,densityvectemp2);
			}

		}
	}
	for (int boxy = 0; boxy < BoxY-1; boxy++)
	{
		for (int boxz = 0; boxz < BoxZ-1; boxz++)
		{
			densityvectemp=(DensityVec.Get(BoxX-1,boxy,boxz)+DensityVec.Get(BoxX-1,boxy,boxz+1)+DensityVec.Get(BoxX-1,boxy+1,boxz)+DensityVec.Get(BoxX-1,boxy+1,boxz+1))/4.0;
			densityvectemp1=(Density1Vec.Get(BoxX-1,boxy,boxz)+Density1Vec.Get(BoxX-1,boxy,boxz+1)+Density1Vec.Get(BoxX-1,boxy+1,boxz)+Density1Vec.Get(BoxX-1,boxy+1,boxz+1))/4.0;
			densityvectemp2=(Density2Vec.Get(BoxX-1,boxy,boxz)+Density2Vec.Get(BoxX-1,boxy,boxz+1)+Density2Vec.Get(BoxX-1,boxy+1,boxz)+Density2Vec.Get(BoxX-1,boxy+1,boxz+1))/4.0;
			DensityVecShiftP.Set(BoxX-1,boxy,boxz,densityvectemp);
			Density1VecShiftP.Set(BoxX-1,boxy,boxz,densityvectemp1);
			Density2VecShiftP.Set(BoxX-1,boxy,boxz,densityvectemp2);
		}
	}
	for (int boxx = 0; boxx < BoxX-1; boxx++)
	{
		for (int boxz = 0; boxz < BoxZ-1; boxz++)
		{
			densityvectemp=(DensityVec.Get(boxx,BoxY-1,boxz)+DensityVec.Get(boxx,BoxY-1,boxz+1)+DensityVec.Get(boxx+1,BoxY-1,boxz)+DensityVec.Get(boxx+1,BoxY-1,boxz+1))/4.0;
			densityvectemp1=(Density1Vec.Get(boxx,BoxY-1,boxz)+Density1Vec.Get(boxx,BoxY-1,boxz+1)+Density1Vec.Get(boxx+1,BoxY-1,boxz)+Density1Vec.Get(boxx+1,BoxY-1,boxz+1))/4.0;
			densityvectemp2=(Density2Vec.Get(boxx,BoxY-1,boxz)+Density2Vec.Get(boxx,BoxY-1,boxz+1)+Density2Vec.Get(boxx+1,BoxY-1,boxz)+Density2Vec.Get(boxx+1,BoxY-1,boxz+1))/4.0;
			DensityVecShiftP.Set(boxx,BoxY-1,boxz,densityvectemp);
			Density1VecShiftP.Set(boxx,BoxY-1,boxz,densityvectemp1);
			Density2VecShiftP.Set(boxx,BoxY-1,boxz,densityvectemp2);
		}
	}
	for (int boxx = 0; boxx < BoxX-1; boxx++)
	{
		for (int boxy = 0; boxy < BoxY-1; boxy++)
		{
			densityvectemp=(DensityVec.Get(boxx,boxy,BoxZ-1)+DensityVec.Get(boxx,boxy+1,BoxZ-1)+DensityVec.Get(boxx+1,boxy,BoxZ-1)+DensityVec.Get(boxx+1,boxy+1,BoxZ-1))/4.0;
			densityvectemp1=(Density1Vec.Get(boxx,boxy,BoxZ-1)+Density1Vec.Get(boxx,boxy+1,BoxZ-1)+Density1Vec.Get(boxx+1,boxy,BoxZ-1)+Density1Vec.Get(boxx+1,boxy+1,BoxZ-1))/4.0;
			densityvectemp2=(Density2Vec.Get(boxx,boxy,BoxZ-1)+Density2Vec.Get(boxx,boxy+1,BoxZ-1)+Density2Vec.Get(boxx+1,boxy,BoxZ-1)+Density2Vec.Get(boxx+1,boxy+1,BoxZ-1))/4.0;
			DensityVecShiftP.Set(boxx,boxy,BoxZ-1,densityvectemp);
			Density1VecShiftP.Set(boxx,boxy,BoxZ-1,densityvectemp1);
			Density2VecShiftP.Set(boxx,boxy,BoxZ-1,densityvectemp2);
		}
	}
	for (int boxx = 0; boxx<BoxX-1 ; boxx++)
	{
		densityvectemp=(DensityVec.Get(boxx,BoxY-1,BoxZ-1)+DensityVec.Get(boxx+1,BoxY-1,BoxZ-1))/2.0;
		densityvectemp1=(Density1Vec.Get(boxx,BoxY-1,BoxZ-1)+Density1Vec.Get(boxx+1,BoxY-1,BoxZ-1))/2.0;
		densityvectemp2=(Density2Vec.Get(boxx,BoxY-1,BoxZ-1)+Density2Vec.Get(boxx+1,BoxY-1,BoxZ-1))/2.0;
		DensityVecShiftP.Set(boxx,BoxY-1,BoxZ-1,densityvectemp);
		Density1VecShiftP.Set(boxx,BoxY-1,BoxZ-1,densityvectemp1);
		Density2VecShiftP.Set(boxx,BoxY-1,BoxZ-1,densityvectemp2);
	}
	for (int boxy = 0; boxy<BoxY-1 ; boxy++)
	{
		densityvectemp=(DensityVec.Get(BoxX-1,boxy,BoxZ-1)+DensityVec.Get(BoxX-1,boxy+1,BoxZ-1))/2.0;
		densityvectemp1=(Density1Vec.Get(BoxX-1,boxy,BoxZ-1)+Density1Vec.Get(BoxX-1,boxy+1,BoxZ-1))/2.0;
		densityvectemp2=(Density2Vec.Get(BoxX-1,boxy,BoxZ-1)+Density2Vec.Get(BoxX-1,boxy+1,BoxZ-1))/2.0;
		DensityVecShiftP.Set(BoxX-1,boxy,BoxZ-1,densityvectemp);
		Density1VecShiftP.Set(BoxX-1,boxy,BoxZ-1,densityvectemp1);
		Density2VecShiftP.Set(BoxX-1,boxy,BoxZ-1,densityvectemp2);
	}
	for (int boxz = 0; boxz<BoxZ-1 ; boxz++)
	{
		densityvectemp=(DensityVec.Get(BoxX-1,BoxY-1,boxz)+DensityVec.Get(BoxX-1,BoxY-1,boxz))/2.0;
		densityvectemp1=(Density1Vec.Get(BoxX-1,BoxY-1,boxz)+Density1Vec.Get(BoxX-1,BoxY-1,boxz))/2.0;
		densityvectemp2=(Density2Vec.Get(BoxX-1,BoxY-1,boxz)+Density2Vec.Get(BoxX-1,BoxY-1,boxz))/2.0;
		DensityVecShiftP.Set(BoxX-1,BoxY-1,boxz,densityvectemp);
		Density1VecShiftP.Set(BoxX-1,BoxY-1,boxz,densityvectemp1);
		Density2VecShiftP.Set(BoxX-1,BoxY-1,boxz,densityvectemp2);
	}
	DensityVecShiftP.Set(BoxX-1,BoxY-1,BoxZ-1,DensityVec.Get(BoxX-1,BoxY-1,BoxZ-1));
	Density1VecShiftP.Set(BoxX-1,BoxY-1,BoxZ-1,Density1Vec.Get(BoxX-1,BoxY-1,BoxZ-1));
	Density2VecShiftP.Set(BoxX-1,BoxY-1,BoxZ-1,Density2Vec.Get(BoxX-1,BoxY-1,BoxZ-1));

	for (int boxx = 0; boxx<BoxX-1; boxx++)
	{
		for (int boxy = 0; boxy<BoxY-1; boxy++)
		{
			densityvectemp=(DensityVec.Get(boxx,boxy,0)+DensityVec.Get(boxx,boxy+1,0)+DensityVec.Get(boxx+1,boxy,0)+DensityVec.Get(boxx+1,boxy+1,0))/4.0;
			densityvectemp1=(Density1Vec.Get(boxx,boxy,0)+Density1Vec.Get(boxx,boxy+1,0)+Density1Vec.Get(boxx+1,boxy,0)+Density1Vec.Get(boxx+1,boxy+1,0))/4.0;
			densityvectemp2=(Density2Vec.Get(boxx,boxy,0)+Density2Vec.Get(boxx,boxy+1,0)+Density2Vec.Get(boxx+1,boxy,0)+Density2Vec.Get(boxx+1,boxy+1,0))/4.0;
			DensityWallVecShiftP.Set(boxx,boxy,densityvectemp);
			DensityWallVec1ShiftP.Set(boxx,boxy,densityvectemp1);
			DensityWallVec2ShiftP.Set(boxx,boxy,densityvectemp2);
		}
	}
	for (int boxy = 0; boxy<BoxY-1; boxy++)
	{
		densityvectemp=(DensityVec.Get(BoxX-1,boxy,0)+DensityVec.Get(BoxX-1,boxy+1,0))/2.0;
		densityvectemp1=(Density1Vec.Get(BoxX-1,boxy,0)+Density1Vec.Get(BoxX-1,boxy+1,0))/2.0;
		DensityWallVecShiftP.Set(BoxX-1,boxy,densityvectemp);
		DensityWallVec1ShiftP.Set(BoxX-1,boxy,densityvectemp1);
		DensityWallVec2ShiftP.Set(BoxX-1,boxy,densityvectemp2);
	}
	for (int boxx = 0; boxx<BoxX-1; boxx++)
	{
		densityvectemp=(DensityVec.Get(boxx,BoxY-1,0)+DensityVec.Get(boxx+1,BoxY-1,0))/2.0;
		densityvectemp1=(Density1Vec.Get(boxx,BoxY-1,0)+Density1Vec.Get(boxx+1,BoxY-1,0))/2.0;
		DensityWallVecShiftP.Set(boxx,BoxY-1,densityvectemp);
		DensityWallVec1ShiftP.Set(boxx,BoxY-1,densityvectemp1);
		DensityWallVec2ShiftP.Set(boxx,BoxY-1,densityvectemp2);
	}
	DensityWallVecShiftP.Set(BoxX-1,BoxY-1,DensityVec.Get(BoxX-1,BoxY-1,0));
	DensityWallVec1ShiftP.Set(BoxX-1,BoxY-1,Density1Vec.Get(BoxX-1,BoxY-1,0));
	DensityWallVec2ShiftP.Set(BoxX-1,BoxY-1,Density2Vec.Get(BoxX-1,BoxY-1,0));

}







// computes surface normal from height array
void GetSurfaceNormal(CoordArray2D& Normal, DoubleArray2D& Height, int minX, int maxX, int minY, int maxY)
{
	//IntCoord2D Size = Height.Size();
	double d_x, d_y, d_z;
	double l_n = 1.0;
	DoubleCoord N;

	// compute normal

    for (int yi = minY*refinementGridHeight; yi < (maxY+1)*refinementGridHeight; yi++)
    {
        for (int xi = minX*refinementGridHeight; xi < (maxX+1)*refinementGridHeight; xi++)
        {

			// gradient of density gives you the surface normal
            d_x = (Height.Get(xi+1,yi) - Height.Get(xi-1,yi))/(2.0*BoxLength/refinementGridHeight);
            d_y = (Height.Get(xi,yi+1) - Height.Get(xi,yi-1))/(2.0*BoxLength/refinementGridHeight);
            d_z = -1.0;
			l_n = sqrt(d_x*d_x + d_y*d_y + d_z*d_z);

			// surface normal is direction of surface tension force
			Normal.At(xi,yi).x = d_x/l_n;
			Normal.At(xi,yi).y = d_y/l_n;
			Normal.At(xi,yi).z = d_z/l_n;
		}
	}
}

// does averaging and smoothing
void BoxAverage(DoubleArray3D& RoughDensity, DoubleArray3D& Density, DoubleArray2D& RoughWallDensity, DoubleArray2D& WallDensity, DoubleArray3D& Filter, int FilterDim, int minx, int maxx, int miny, int maxy, int maxz, int BoxX, int BoxY, int BoxZ)
{
	// need smaller basis for edges
	double avg = 0;
	DoubleArray3D InArray(RoughDensity.Size().x, RoughDensity.Size().y, RoughDensity.Size().z+1);
	DoubleArray3D OutArray(RoughDensity.Size().x, RoughDensity.Size().y, RoughDensity.Size().z+1);
	for (int iset=0;iset<BoxX;iset++)
	{
		for (int jset=0; jset<BoxY; jset++)
		{
			InArray.Set(iset,jset,0,RoughWallDensity.Get(iset,jset));
			for (int kset=1; kset<BoxZ+1; kset++)
			{
				InArray.Set(iset,jset,kset,RoughDensity.Get(iset,jset,kset-1));
			}
		}
	}

	int i,j,k;

	for (int kArray = 0; kArray<maxz+1; kArray++)
	{
		for (int jArray = miny; jArray<maxy; jArray++)
		{
			for (int iArray = minx; iArray<maxx; iArray++)
			{
				avg = 0.0;
				for (int iFilter = -FilterDim; iFilter<=FilterDim; iFilter++)
				{
					for (int jFilter = -FilterDim; jFilter<=FilterDim; jFilter++)
					{
						for (int kFilter = -FilterDim; kFilter<=FilterDim; kFilter++)
						{

							i = iArray+iFilter;
							j = jArray+jFilter;
							k = kArray+kFilter;

							if (k<0)
							{
								if (InArray.Get(i,j,0)>0.00001)
								avg += Filter.Get(iFilter+FilterDim,jFilter+FilterDim,kFilter+FilterDim)*1.0;//InArray.Get(i,j,0);   //*density_threshold; 应当修改的地方
								else
									avg += Filter.Get(iFilter+FilterDim,jFilter+FilterDim,kFilter+FilterDim)*InArray.Get(i,j,0);
							}
							else
								avg += Filter.Get(iFilter+FilterDim,jFilter+FilterDim,kFilter+FilterDim)*InArray.Get(i,j,k);

						}
					}
				}

				//OutArray.Set(iArray,jArray,kArray,avg);
				OutArray.Set(iArray,jArray,kArray,min(1.0,avg/density_threshold));////修改！！！
			}
		}
	}
    for (int iset=0;iset<BoxX;iset++)
	{
		for (int jset=0; jset<BoxY; jset++)
		{
			WallDensity.Set(iset,jset,OutArray.Get(iset,jset,0));
			for (int kset=1; kset<BoxZ+1; kset++)
			{
				Density.Set(iset,jset,kset-1,OutArray.Get(iset,jset,kset));
			}
		}
	}
}


// averages height
void Height_Average(DoubleArray2D& InHeight, DoubleArray2D& OutHeight, int minx, int maxx, int miny, int maxy)
{
	// need smaller basis for edges
	double avg;
	double H0, H, Ht;
	int count;

    for (int jArray = miny*refinementGridHeight; jArray<(maxy+1)*refinementGridHeight; jArray++)
    {
        for (int iArray = minx*refinementGridHeight; iArray<(maxx+1)*refinementGridHeight; iArray++)
        {
			H0 = InHeight.Get(iArray,jArray);
			Ht = H0;
			count = 1;
			for (int iNeighbours = -1; iNeighbours<2; iNeighbours++)
			{
				for (int jNeighbours = -1; jNeighbours<2; jNeighbours++)
				{
					if (!((iNeighbours==0)&&(jNeighbours==0)))
					{
						H = InHeight.Get(iArray+iNeighbours,jArray+jNeighbours);
						if (H<H0)
						{
							Ht += H;
							count++;
						}

					}
				}
			}

			avg = Ht/count;
			OutHeight.Set(iArray,jArray,avg);
		}
	}

}

// smooths any array
void Smooth(DoubleArray2D& InArray, DoubleArray2D& OutArray, DoubleArray2D& Filter, int FilterDim, int minx, int maxx, int miny, int maxy)
{
	// need smaller basis for edges
	double avg, h;

	int i,j;

	for (int jArray = miny*refinementGridHeight; jArray<(maxy+1)*refinementGridHeight; jArray++)
	{
		for (int iArray = minx*refinementGridHeight; iArray<(maxx+1)*refinementGridHeight; iArray++)
		{
			avg = 0.0;

			for (int iFilter = -FilterDim; iFilter<=FilterDim; iFilter++)
			{
				for (int jFilter = -FilterDim; jFilter<=FilterDim; jFilter++)
				{

						i = iArray+iFilter;
						j = jArray+jFilter;

						h = InArray.Get(i,j);
						avg += Filter.Get(iFilter+FilterDim,jFilter+FilterDim)*h;

				}
			}
			OutArray.Set(iArray,jArray,avg);
		}
	}

}
