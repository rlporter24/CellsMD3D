#include <math.h>
#include <float.h>

#include "Array.h"
#include "Cell.h"
#include "Constants.h"
#include "Forces.h"
#include "tools.h"
#include "Neighbours.h"

// Cell-cell forces
// returns force F, position of the contact force r, and distance between cells d
void F_cc(const Cell& cell1, const Cell& cell2, DoubleCoord& F, DoubleCoord& r, double& d)
{
	// the points indicating the segment of closest approach of the two cells
	DoubleCoord c1, c2;
	double gamma_n = 100.0;

	// touching distance
	double sigma = cell1.Radius+cell2.Radius;

	// inertia
	double L1 = cell1.Length + 4.0/3*cell1.Radius;	// Length
	double M1 = L1*cell1.Radius*cell1.Radius*PI;	// Mass (approximate)

	double L2 = cell2.Length + 4.0/3*cell2.Radius;	// Length
	double M2 = L2*cell2.Radius*cell2.Radius*PI;	// Mass (approximate)

	double Meff = M1*M2/(M1+M2);

	// find overlap of cells
	min_distance(cell1,cell2,d,c1,c2);	// minimum distance between cell segments
	double delta = max(sigma-d,0.0);	// overlap

	// normal and tangent unit vectors
	DoubleCoord vn = scale(diff(c2,c1),1.0/d);	// normal vector

	// elastic force, location of overlap
//	double FE_mag = k_cc*pow(delta,3.0/2.0); // Mya's approximation
    double FE_mag = sqrt(cell1.Radius*cell2.Radius/(cell1.Radius+cell2.Radius))*k_cc*pow(delta,3.0/2.0); // Yue's version
	DoubleCoord FE = scale(vn,-FE_mag);	// elastic force

//	r = scale(sum(c1,c2),cell1.Radius/sigma);	// position of the force (if F!=0, this should be on the cell boundary, or close to it), Mya's version
    r = scale(sum(scale(c1,cell2.Radius),scale(c2,cell1.Radius)),1/sigma); // changed from last line because cell1 and cell2 have different radii, Yue's version.

	if (FE_mag>DBL_EPSILON)	// if cells are touching, compute friction
	{

		// angular velocity
		DoubleCoord rcm1 = diff(r,average(cell1.Position));
		DoubleCoord rcm2 = diff(r,average(cell2.Position));

		DoubleCoord va1 = cross(cell1.AngularVelocity,rcm1);
		DoubleCoord va2 = cross(cell2.AngularVelocity,rcm2);

		// total velocity
		DoubleCoord v1 = sum(cell1.Velocity,va1);
		DoubleCoord v2 = sum(cell2.Velocity,va2);

		// difference in velocity
		DoubleCoord dv = diff(v2,v1);

		double gamma = gamma_t;
		double mu = cell_mu;

		// tangential and normal parts of the dissipative force
		DoubleCoord dvn = scale(vn,dot(dv,vn));
		DoubleCoord FDn = scale(dvn,gamma_n*Meff*delta);		// normal dissipation

		DoubleCoord dvt = diff(dv,dvn);					// tangential dissipation
		double v_mag = sqrt(dot(dvt,dvt));

		DoubleCoord FDt(0.0,0.0,0.0);

		if (v_mag>3*DBL_EPSILON)
			FDt = scale(dvt,min(gamma*Meff*sqrt(delta),mu*FE_mag/v_mag));


		// sum total force, elastic + dissipative
		F = sum(sum(FE,FDt),FDn);
	}
	else
		F = FE;
}

// Force between cell and agar (wall) at y=0
// returns forces on each vertex, F1 and F2, positions of the two forces, r1 and r2
void F_cw(const Cell& cell, double Wall_z, DoubleCoord& F1, DoubleCoord& F2, DoubleCoord& r1, DoubleCoord& r2, DoubleCoord& cwStaFric, DoubleCoord& cwDynFric)
{
	double L = cell.Length + 4.0/3*cell.Radius;	// Length
	double M = L*cell.Radius*cell.Radius*PI;	// Mass (approximate)
	double gamma_n = 100.0;

	// **** First Locus ******

	// Elastic

	double sigma = cell.Radius;

	// Location of wall at first locus including effect of roughness
	DoubleCoord wall;
	wall.x = cell.Position.p.x+((float)rand()/RAND_MAX-0.5)*var_pos;
	wall.y = cell.Position.p.y+((float)rand()/RAND_MAX-0.5)*var_pos;
	wall.z = Wall_z;
	r1 = wall; // location of wall force

	// Calc elastic force from overlap
	double d = cell.Position.p.z-wall.z;
	double delta = max(sigma-d,0.0);	// overlap
//	double FE_mag = k_wc*pow(delta,3.0/2.0);	// elastic force, Mya's version
    double FE_mag = sqrt(cell.Radius)*k_wc*pow(delta,3.0/2.0);   // Yue's version

	DoubleCoord FE(0.0,0.0,0.0);
	FE.z = FE_mag;


	// Dissipative
	DoubleCoord rcm, v, va, dv, dvn, dvt, FDt, FDn;
	double v_mag;

	if (FE_mag>DBL_EPSILON)
	{
		rcm = diff(r1,average(cell.Position));
		va = cross(cell.AngularVelocity,rcm);

		// total velocity
		v = sum(cell.Velocity,va);

		// tangential and normal parts of the dissipative force
		dv = scale(v,-1.0);
		dvn = DoubleCoord(0.0,0.0,dv.z);

		dvt = diff(dv,dvn);
		v_mag = sqrt(dot(dvt,dvt));

		if (v_mag>DBL_EPSILON)
			FDt = scale(dvt,min(gamma_t*M*sqrt(delta),wall_mu*FE_mag/v_mag));	// tangential dissipation
		else FDt = DoubleCoord(0,0,0);

		FDn = scale(dvn,gamma_n*M*delta);		// normal dissipation

		// sum total force, elastic + dissipative
		F1 = sum(FE,sum(FDt,FDn));

		cwStaFric = sum(cwStaFric,scale(dvt,gamma_t*M*sqrt(delta))); 
		cwDynFric = sum(cwDynFric,scale(dvt,wall_mu*FE_mag/v_mag));

	}
	else
	{
		F1 = DoubleCoord(0.0,0.0,0.0);
	}

	// **** Second locus ****


	wall.x = cell.Position.q.x;
	wall.y = cell.Position.q.y;
	wall.z = Wall_z;
	r2 = wall; // location of wall force

	d = cell.Position.q.z-wall.z;
	delta = max(sigma-d,0.0);	// overlap
//	FE_mag = k_wc*pow(delta,3.0/2.0);	// elastic force, Mya's version
    FE_mag = sqrt(cell.Radius)*k_wc*pow(delta,3.0/2.0); // Yue's version

	FE.z = FE_mag;


	if (FE_mag>DBL_EPSILON)
	{
		rcm = diff(r2,average(cell.Position));
		va = cross(cell.AngularVelocity,rcm);

		v = sum(cell.Velocity,va);

		// tangential and normal parts of the dissipative force
		dv = scale(v,-1.0);
		dvn = DoubleCoord(0.0,0.0,dv.z);

		dvt = diff(dv,dvn);
		v_mag = sqrt(dot(dvt,dvt));

		if (v_mag>DBL_EPSILON)
			FDt = scale(dvt,min(gamma_t*M*sqrt(delta),wall_mu*FE_mag/v_mag));	// tangential dissipation
		else FDt = DoubleCoord(0,0,0);

		FDn = scale(dvn,gamma_n*M*delta);		// normal dissipation

		// sum total force, elastic + dissipative
		F2 = sum(FE,sum(FDt,FDn));

		cwStaFric = sum(cwStaFric,scale(dvt,gamma_t*M*sqrt(delta)));
                cwDynFric = sum(cwDynFric,scale(dvt,wall_mu*FE_mag/v_mag));
		
	}

	else
	{
		F2 = DoubleCoord(0.0,0.0,0.0);
	}


}

// surface tension force
// takes the Height of the water field, the inward Normal vector
// returns force F and torque T
void F_surf_tension(const Cell& cell, UniformGrid& Grid, const IntCoord& XYAddress, DoubleArray2D& Height, CoordArray2D& Normal, DoubleCoord& F, DoubleCoord& T)
{
	DoubleCoord F1(0,0,0), F2(0,0,0);
	double distance;
	double dh, F0;

//	int xa = XYAddress.x;
//	int ya = XYAddress.y;
    int xa = int((average(cell.Position).x+BoxX/2*BoxLength)/(BoxLength/refinementGridHeight));
    int ya = int((average(cell.Position).y+BoxY/2*BoxLength)/(BoxLength/refinementGridHeight));
    
    int i_ht, j_ht;
    
	// find interpolated water position at the x, y coordinate of the p vertex
    DoubleCoord dxp;// = scale(diff(cell.Position.p,Grid.GetCentroid(XYAddress)),1.0/(BoxLength/refinementGridHeight));
    i_ht = int(cell.Position.p.x/(BoxLength/refinementGridHeight)+BoxX/2*refinementGridHeight);
    j_ht = int(cell.Position.p.y/(BoxLength/refinementGridHeight)+BoxY/2*refinementGridHeight);
    dxp.x = cell.Position.p.x/(BoxLength/refinementGridHeight) - (i_ht-BoxX/2*refinementGridHeight+0.5);
    dxp.y = cell.Position.p.y/(BoxLength/refinementGridHeight) - (j_ht-BoxY/2*refinementGridHeight+0.5);
    
    //	int x0 = ceil(dxp.x)+xa-1;
    //	int y0 = ceil(dxp.y)+ya-1;
    i_ht = ceil(dxp.x)+i_ht-1;
    j_ht = ceil(dxp.y)+j_ht-1;
    
    double height_p = Height.linear_interp(i_ht, j_ht, dxp.x, dxp.y);//Height.linear_interp(x0, y0, dxp.x, dxp.y);
    
    // find interpolated water position at the x, y coordinate of the q vertex
    DoubleCoord dxq;// = scale(diff(cell.Position.q,Grid.GetCentroid(XYAddress)),1.0/BoxLength);
    i_ht = int(cell.Position.q.x/(BoxLength/refinementGridHeight)+BoxX/2*refinementGridHeight);
    j_ht = int(cell.Position.q.y/(BoxLength/refinementGridHeight)+BoxY/2*refinementGridHeight);
    dxq.x = cell.Position.q.x/(BoxLength/refinementGridHeight) - (i_ht-BoxX/2*refinementGridHeight+0.5);
    dxq.y = cell.Position.q.y/(BoxLength/refinementGridHeight) - (j_ht-BoxY/2*refinementGridHeight+0.5);
    
    //	x0 = ceil(dxq.x)+xa-1;
    //	y0 = ceil(dxq.y)+ya-1;
    i_ht = ceil(dxq.x)+i_ht-1;
    j_ht = ceil(dxq.y)+j_ht-1;
    
    double height_q = Height.linear_interp(i_ht, j_ht, dxq.x, dxq.y);//Height.linear_interp(x0, y0, dxq.x, dxq.y);
    
    // find interpolated water position at the x, y coordinate of the center of mass
    DoubleCoord cm = average(cell.Position);
    DoubleCoord dx;// = scale(diff(cm,Grid.GetCentroid(XYAddress)),1.0/BoxLength);
    i_ht = int(cm.x/(BoxLength/refinementGridHeight)+BoxX/2*refinementGridHeight);
    j_ht = int(cm.y/(BoxLength/refinementGridHeight)+BoxY/2*refinementGridHeight);
    dx.x = cm.x/(BoxLength/refinementGridHeight) - (i_ht-BoxX/2*refinementGridHeight+0.5);
    dx.y = cm.y/(BoxLength/refinementGridHeight) - (j_ht-BoxY/2*refinementGridHeight+0.5);
    
    i_ht = ceil(dx.x)+i_ht-1;
    j_ht = ceil(dx.y)+j_ht-1;
    
    double height_cm = Height.linear_interp(i_ht, j_ht, dx.x, dx.y);

	// if the height of the water is far above the cell height, there is no force
	if (height_cm>(cm.z+BoxLength))
	{

		F = DoubleCoord(0,0,0);
		T = DoubleCoord(0,0,0);

	}
	else
	{

		// dh is the relative height of the vertex compared to the water level
		dh = cell.Position.p.z - (height_p-DH/Normal.At(xa,ya).z);

		// if the cell is not above the water height (dh<0) then there is no force
		if (dh<0.0)
		{
			F1 = DoubleCoord(0,0,0);
		}
		else
		{

			distance = -dh*Normal.At(xa,ya).z;

			// slightly ad hoc form of the magnitude of the surface tension force
			F0 = 2*PI*tension * min(distance/(cellRadius/5.0),1.0);
			// F1 is a vector in the normal direction
			F1 = scale(Normal.At(xa,ya),F0);
		}

		// do the same for the q vertex
		dh = cell.Position.q.z - (height_q-DH/Normal.At(xa,ya).z);

		if (dh<0.0)
		{
			F2 = DoubleCoord(0,0,0);
		}
		else
		{

			distance = -dh*Normal.At(xa,ya).z;

			F0 = 2*PI*tension * min(distance/(cellRadius/5.0),1.0);
			F2 = scale(Normal.At(xa,ya),F0);
		}

		F = sum(F1,F2);

		// Calculate torque
		DoubleCoord r = diff(cell.Position.p, cm);
		T = cross(r,F1);

		r = diff(cell.Position.q, cm);
		T = sum(T,cross(r,F2));

	}

}

// viscous force between cell and surrounding liquid damps the cell velocity
void F_v(const Cell& cell, DoubleCoord& F, DoubleCoord& T)
{
//	double M = cell.Length*PI*cell.Radius*cell.Radius + 4.0/3.0*PI*cell.Radius*cell.Radius*cell.Radius;	// Mass (approximate), Mya's version

	// find viscous drag with fluid
//	F = scale(cell.Velocity,-viscosity*M);	// easier to move in direction of cell length; Mya's version
    F = scale(cell.Velocity,-viscosity*6.0*PI*cell.Radius);   // Yue's version
//	T = scale(cell.AngularVelocity,-viscosity*M);  // Mya's version
    T = scale(cell.AngularVelocity,-viscosity*6.0*PI*cell.Radius);   // Yue's version
}

// sum all of the forces on the cell
void sum_forces(const Cell& cell, const Cell* cell_array, const int* neighbours, DoubleCoord& Fnet, DoubleCoord& Tnet, DoubleArray2D& Height, CoordArray2D& Normal, UniformGrid& Grid, const IntCoord& XYAddress, DoubleArray2D& Wall, DoubleCoord& cwStaFric, DoubleCoord& cwDynFric)
{
	Fnet = DoubleCoord(0,0,0);
	Tnet = DoubleCoord(0,0,0);
	DoubleCoord F(0,0,0), F2(0,0,0), T(0,0,0), r, r2;
	DoubleCoord cm = average(cell.Position);
	double d, wall_y;

	int ID;
	int numNeighbours = neighbours[0];

	// loop through neighbours and find the forces
	for (int neighbourID = 1; neighbourID<numNeighbours+1; neighbourID++)
	{
		ID = neighbours[neighbourID];	// the ID of the current neighbour
		F_cc(cell, cell_array[ID], F, r, d);	// contact force between cell and neighbours
		Fnet = sum(Fnet, F);	// net force is the sum of all forces
		r = diff(r, cm);		// r is distance from centre of mass to contact force location
		Tnet = sum(Tnet, cross(r, F));	// net torque = sum(rxF)
	}

	cwStaFric = DoubleCoord(0,0,0);
	cwDynFric = DoubleCoord(0,0,0);
    // calculate cell-wall forces if cell is close to wall
	if (min(cell.Position.q.z,cell.Position.p.z)<1.2*cell.Radius)
	{
		wall_y = Wall.Get(XYAddress.x,XYAddress.y);
		F_cw(cell, wall_y, F, F2, r, r2, cwStaFric, cwDynFric);
		Fnet = sum(Fnet, F);	// net force is the sum of all forces
		r = diff(r, cm);		// r is distance from centre of mass to contact force location
		Tnet = sum(Tnet, cross(r, F));	// net torque = sum(rxF)

		Fnet = sum(Fnet, F2);	// net force is the sum of all forces
		r2 = diff(r2, cm);		// r is distance from centre of mass to contact force location
		Tnet = sum(Tnet, cross(r2, F2));	// net torque = sum(rxF)
	}
    
	// viscous force with fluid
	F_v(cell, F, T);
	Fnet = sum(Fnet, F);
	Tnet = sum(Tnet, T);

	F_surf_tension(cell, Grid, XYAddress, Height, Normal, F, T);
	Fnet = sum(Fnet, F);
	Tnet = sum(Tnet, T);
    

}

