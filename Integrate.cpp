#include <math.h>
#include <iostream>

#include "Cell.h"
#include "Forces.h"
#include "Integrate.h"
#include "tools.h"


// integrates to find the new positions based on the forces due to neighbors
void integrate(double dt, int cellID, const Cell* old_cells, Cell* new_cells, const int* neighbours, DoubleArray2D& Height, CoordArray2D& Normal, UniformGrid& Grid, const IntCoord& XYAddress, DoubleArray2D& Wall)
{
  DoubleCoord F, T, cwStaFric, cwDynFric;

	// cell we're computing forces on
	Cell temp_cell;

	// sum forces from original positions, velocities
	sum_forces(old_cells[cellID], old_cells, neighbours, F, T, Height, Normal, Grid, XYAddress, Wall, cwStaFric, cwDynFric);	// Find F, T

	// update positions of temp cell, velocities updated 1/2 step
	UpdatePositions(dt, F, T, old_cells[cellID], temp_cell);

	// Recalculate dissipative forces
	sum_forces(temp_cell, old_cells, neighbours, F, T, Height, Normal, Grid, XYAddress, Wall, cwStaFric, cwDynFric);

	// Update velocities another 1/2 step
	UpdateVelocities(dt, F, T, temp_cell, new_cells[cellID]);
	new_cells[cellID].DynFric = cwDynFric;
	new_cells[cellID].StaFric = cwStaFric;
}

// Integrate positions forward in time by dt, velocities by dt/2
void UpdatePositions(double dt, const DoubleCoord& F, DoubleCoord& T, const Cell& old_cell, Cell& new_cell)
{
	// sets all properties to be the same
	new_cell = old_cell;

	// Inertia
//	double L = old_cell.Length + 2*old_cell.Radius;	// Length, Mya's version
//	double M = old_cell.Length*PI*old_cell.Radius*old_cell.Radius + 4.0/3.0*PI*old_cell.Radius*old_cell.Radius*old_cell.Radius;	// Mass (approximate), Mya's version

//	double MR2 = M*old_cell.Radius*old_cell.Radius;// Mya's version

//	double Ixy = (3*MR2+M*L*L)*0.5;		// Moment of inertia, Mya's version
//	double Iz = MR2;   // Mya's version
    double R2 = old_cell.Radius*old_cell.Radius; // Yue and Hui's version
    double L2 = old_cell.Length*old_cell.Length; // Yue and Hui's version
    double M = PI*R2*(old_cell.Length + 4.0/3.0*old_cell.Radius); // Yue and Hui's version
    double Ixy = PI*R2*(old_cell.Length*L2/12.0+old_cell.Radius*(8.0*R2/15.0+L2/3.0));  // Yue and Hui's version
    double Iz = PI*R2*R2*(old_cell.Length/2.0+8.0*old_cell.Radius/15.0);   // Yue and Hui's version

	// unit vector along length of cell
	DoubleCoord z = diff(old_cell.Position.q, old_cell.Position.p);	// z axis (through length of cell)
	z = scale(z,1.0/old_cell.Length);

	// angular momentum along cell axis (z)
	DoubleCoord vaz = scale(z,dot(z,old_cell.AngularVelocity));
	DoubleCoord vaxy = diff(old_cell.AngularVelocity,vaz);
	DoubleCoord Tz = scale(z,dot(z,T));
	DoubleCoord Txy = diff(T,Tz);

	// update velocities dt/2
	new_cell.Velocity = sum(old_cell.Velocity,scale(F,dt/(2.0*M)));
	new_cell.AngularVelocity = sum(sum(vaz,scale(Tz,dt/(2.0*Iz))), sum(vaxy,scale(Txy,dt/(2.0*Ixy))));

	// change cell coordinates to reduced coordinates
	DoubleCoord cm = average(old_cell.Position);	// centre of mass

	// vector to nodes
	DoubleCoord rcmp = diff(old_cell.Position.p,cm);
	DoubleCoord rcmq = diff(old_cell.Position.q,cm);

	// angular velocity
	DoubleCoord vap = cross(new_cell.AngularVelocity,rcmp);
	DoubleCoord vaq = cross(new_cell.AngularVelocity,rcmq);

	// total velocity
	DoubleCoord vp = sum(new_cell.Velocity,vap);
	DoubleCoord vq = sum(new_cell.Velocity,vaq);

	// Integrate to find new positions
	new_cell.Position.p = sum(old_cell.Position.p, scale(vp,dt));
	new_cell.Position.q = sum(old_cell.Position.q, scale(vq,dt));

}

// Move velocities forward dt/2
void UpdateVelocities(double dt, const DoubleCoord& F, DoubleCoord& T, const Cell& old_cell, Cell& new_cell)
{
	// sets all properties to be the same
	new_cell = old_cell;

	// Inertia
//	double L = old_cell.Length + 2*old_cell.Radius;	// Length, Mya's version
//	double M = L*PI*old_cell.Radius*old_cell.Radius;	// Mass (approximate), Mya's version
//	double MR2 = M*old_cell.Radius*old_cell.Radius; // Mya's version

//	double Ixy = (3*MR2+M*L*L)*0.5;		// Moment of inertia, Mya's version
//	double Iz = MR2;     // Mya's version
    double R2 = old_cell.Radius*old_cell.Radius; // Yue and Hui's version
    double L2 = old_cell.Length*old_cell.Length; // Yue and Hui's version
    double M = PI*R2*(old_cell.Length + 4.0/3.0*old_cell.Radius); // Yue and Hui's version
    double Ixy = PI*R2*(old_cell.Length*L2/12.0+old_cell.Radius*(8.0*R2/15.0+L2/3.0));  // Yue and Hui's version
    double Iz = PI*R2*R2*(old_cell.Length/2.0+8.0*old_cell.Radius/15.0);   // Yue and Hui's version

	// unit vector along length of cell
	DoubleCoord z = diff(old_cell.Position.q, old_cell.Position.p);	// z axis (through length of cell)
	z = scale(z,1.0/old_cell.Length);

	// angular momentum along cell axis (z)
	DoubleCoord vaz = scale(z,dot(z,old_cell.AngularVelocity));
	DoubleCoord vaxy = diff(old_cell.AngularVelocity,vaz);
	DoubleCoord Tz = scale(z,dot(z,T));
	DoubleCoord Txy = diff(T,Tz);

	// update velocities dt/2
	new_cell.Velocity = sum(old_cell.Velocity,scale(F,dt/(2.0*M)));
	new_cell.AngularVelocity = sum(sum(vaz,scale(Tz,dt/(2.0*Iz))), sum(vaxy,scale(Txy,dt/(2.0*Ixy))));

}



