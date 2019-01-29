/*
 * grow.h
 *
 *  Created on: Mar 1, 2010
 *      Author: mya
 */

#ifndef GROW_H_
#define GROW_H_

#include "Array.h"
#include "Cell.h"

void grow(double dt, Cell& cell, EnvArray3D& Env, AgaArray2D** Wal, UniformGrid& Grid);
void divide(Cell& mother, Cell& daughter, double t);

#endif /* GROW_H_ */
