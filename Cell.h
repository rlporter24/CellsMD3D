/*
 * Cell.h
 *
 *  Created on: Dec 15, 2011
 *      Author: mya
 */

#ifndef CELL_H_
#define CELL_H_

#include "tools.h"


// Cell structure contains information specific to each cell
struct Cell
{
	Segment Position; 	// a segment has two coordinates p and q for the two vertices of the spherocylinder
	double Length;
	double Radius;
	double GrowthRate;
	DoubleCoord Velocity;
	DoubleCoord AngularVelocity; // now a vector
	DoubleCoord DynFric;
	DoubleCoord StaFric;
	int Type;
	int Ancestor;
};

#endif /* CELL_H_ */
