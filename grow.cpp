#include <math.h>
#include <stdlib.h>
#include <float.h>


#include "Array.h"
#include "Cell.h"
#include "Constants.h"
#include "Forces.h"
#include "grow.h"
#include "tools.h"
#include "UniformGrid.h"
#include "Nutrients.h"


// grow the cell
void grow(double dt, Cell& cell, EnvArray3D& Env, AgaArray2D** Wal, UniformGrid& Grid)
{
	DoubleCoord cm = average(cell.Position);	// center of mass
	double Growth_rate;

	// Get position in uniform grid to access correct index for height
	IntCoord XYAddress = Grid.GetXY(Grid.GetAddress(cm));
    
	// Look up the growth rate in the environment array
	Growth_rate = Env.Get(XYAddress).GrowthRate;
	
    if (XYAddress.z==0)
    {
        double lambda=cm.z/BoxLength-floor(cm.z/BoxLength);
        double y0, y1;
        double Cgr;
        y0=(Wal[0]->Get(XYAddress.x,XYAddress.y).CarbonAgar+Wal[0]->Get(XYAddress.x-1,XYAddress.y).CarbonAgar+Wal[0]->Get(XYAddress.x,XYAddress.y-1).CarbonAgar+Wal[0]->Get(XYAddress.x-1,XYAddress.y-1).CarbonAgar)/4;
        y1=(Env.Get(XYAddress.x,XYAddress.y,XYAddress.z).Carbon+Env.Get(XYAddress.x-1,XYAddress.y,XYAddress.z).Carbon+Env.Get(XYAddress.x,XYAddress.y-1,XYAddress.z).Carbon+Env.Get(XYAddress.x-1,XYAddress.y-1,XYAddress.z).Carbon)/4;
        Cgr = y1*lambda+y0*(1-lambda);
        Cgr = Cgr/(Cgr+KC);
        Growth_rate = maxGrowthRate*max(0.0,Cgr-Maintenance_rate/C_rate);
    }

    if (XYAddress.z>0)
    {
        double lambda=cm.z/BoxLength-floor(cm.z/BoxLength);
        double y0, y1;
        double Cgr;
        y0=(Env.Get(XYAddress.x,XYAddress.y,XYAddress.z-1).Carbon+Env.Get(XYAddress.x-1,XYAddress.y,XYAddress.z-1).Carbon+Env.Get(XYAddress.x,XYAddress.y-1,XYAddress.z-1).Carbon+Env.Get(XYAddress.x-1,XYAddress.y-1,XYAddress.z-1).Carbon)/4;
        y1=(Env.Get(XYAddress.x,XYAddress.y,XYAddress.z).Carbon+Env.Get(XYAddress.x-1,XYAddress.y,XYAddress.z).Carbon+Env.Get(XYAddress.x,XYAddress.y-1,XYAddress.z).Carbon+Env.Get(XYAddress.x-1,XYAddress.y-1,XYAddress.z).Carbon)/4;
        Cgr = y1*lambda+y0*(1-lambda);
        Cgr = Cgr/(Cgr+KC);
        Growth_rate = maxGrowthRate*max(0.0,Cgr-Maintenance_rate/C_rate);
    }

	Growth_rate = Growth_rate*1.5850; // log(3)/log(2)

	DoubleCoord v = diff(cell.Position.q, cell.Position.p);				// vector along segment

	// growth length
	double dL = cell.Length*dt*Growth_rate;

	// growth direction
	DoubleCoord dv = scale(v,dL*0.5/cell.Length);

	cell.Length += dL; // increase the length of the cell

	cell.Position.p = diff(cell.Position.p,dv);
	cell.Position.q = sum(cell.Position.q,dv);
	cell.GrowthRate = Growth_rate;
	if (cell.Ancestor==0)
			printf("wrong growth! Ancestor is zero! \n");

}



// takes in a mother cell and returns mother and daughter after division
void divide(Cell& mother, Cell& daughter, double t){

	double dl = ((float)rand()/RAND_MAX-0.5)*varL;	// random part of the length after division
    double dangle = ((float)rand()/RAND_MAX-0.5)*varAngle;	// random part of the orientation angle after division
    DoubleCoord pq = diff(mother.Position.p,mother.Position.q);
    // random orientation for cell1
    double ddirx1 = (float)rand()/RAND_MAX-0.5;	// random part of the orientation direction after division
    double ddiry1 = (float)rand()/RAND_MAX-0.5;	// random part of the orientation direction after division
    double ddirz1 = (float)rand()/RAND_MAX-0.5;	// random part of the orientation direction after division
    //double ddirz1=0.0;
    DoubleCoord ddir1 = DoubleCoord(ddirx1,ddiry1,ddirz1);
    ddir1 = diff(ddir1,scale(pq,dot(ddir1,pq)/dot(pq,pq)));
    ddir1 = scale(ddir1,1/sqrt(dot(ddir1,ddir1)));
    DoubleCoord dOrien1=scale(ddir1,dangle);
    // random orientation for cell2
    double ddirx2 = (float)rand()/RAND_MAX-0.5;	// random part of the orientation direction after division
    double ddiry2 = (float)rand()/RAND_MAX-0.5;	// random part of the orientation direction after division
    double ddirz2 = (float)rand()/RAND_MAX-0.5;	// random part of the orientation direction after division
    //double ddirz2=0.0;
    DoubleCoord ddir2 = DoubleCoord(ddirx2,ddiry2,ddirz2);
    ddir2 = diff(ddir2,scale(pq,dot(ddir2,pq)/dot(pq,pq)));
    ddir2 = scale(ddir2,1/sqrt(dot(ddir2,ddir2)));
    DoubleCoord dOrien2=scale(ddir2,dangle);

	DoubleCoord uv = diff(mother.Position.p,mother.Position.q); // vector in between two coordinates of the cell
	double L = mother.Length;	// length of original cell
	uv = scale(uv,1/L);	// scale unit vector by length to get unit vector

	DoubleCoord midpoint;
	midpoint = average(mother.Position);

	// calculate Segments for daughter cells
	// p in first daughter cell is the same as mother cell
	// q in second daughter cell is the same as in the mother cell
	daughter.Position.q = mother.Position.q;

	// set new lengths
	daughter.Radius = cellRadius;

	mother.Length = L/2-mother.Radius+dl;

	daughter.Length = L/2+mother.Radius-2*daughter.Radius-dl;
	daughter.GrowthRate = mother.GrowthRate;
	daughter.Velocity = mother.Velocity;

	daughter.AngularVelocity = DoubleCoord(0,0,((float)rand()/RAND_MAX-0.5)*0.001);
	mother.AngularVelocity = DoubleCoord(0,0,0);

	// q in first daughter cell
    mother.Position.q = diff(mother.Position.p,scale(uv,mother.Length));
    daughter.Position.time_q = mother.Position.time_q;
    mother.Position.time_q = t;
    daughter.Position.time_p = t;
    mother.Position.age_p++;
    daughter.Position.age_q = mother.Position.age_q+1;
    //mother.Position.time_q = 0;
    //daughter.Position.time_p = 0;

	// p in second daughter cell
	daughter.Position.p = diff(mother.Position.q, scale(uv,daughter.Radius+mother.Radius));
	daughter.Position.q = diff(daughter.Position.p, scale(uv,daughter.Length));
    
    // add randomness in rotation of new p and q
    mother.Position.p=sum(mother.Position.p,scale(dOrien1,mother.Length));
    mother.Position.q=diff(mother.Position.q,scale(dOrien1,mother.Length));
    daughter.Position.p=sum(daughter.Position.p,scale(dOrien2,daughter.Length));
    daughter.Position.q=diff(daughter.Position.q,scale(dOrien2,daughter.Length));

	// daughter inherits its Type from its mother
	daughter.Type = mother.Type;
	if (mother.Ancestor==0)
		printf("wrong mother! Ancestor is zero! \n");
	daughter.Ancestor = mother.Ancestor;


}

