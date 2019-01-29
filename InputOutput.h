#ifndef INPUTOUTPUT_H_

#define INPUTOUTPUT_H_



#include <iostream>

#include "tools.h"
#include "UniformGrid.h"

struct OutputFiles
{
	FILE* cells;
	FILE* restart;
	FILE* roughDensity;
	FILE* roughDensity1;
	FILE* roughDensity2;
	FILE* density;
	FILE* density1;
	FILE* density2;
	FILE* walldensity;
	FILE* walldensity1;
    FILE* walldensity2;
    FILE* roughheight;
	FILE* height;
	FILE* normal;
	FILE* env;
	FILE* aga;
    FILE* wal;
    FILE* lineage;
};

struct Inputs
{
	int ColonyNumber;
	double ColonySeparation;
	double ColonyRadius;// added by YueYan
	int ColonySize;// added by YueYan
	double Fraction;
};

struct Cell;

void CreateOutputFileLineage(int OutputID, OutputFiles& Files, bool append);
void CloseOutputFileLineage(OutputFiles& Files);
void CreateOutputFiles(int OutputID, OutputFiles& Files, bool append);
void CloseOutputFiles(OutputFiles& Files);
int AddFirstCells(Cell* cells, double L_divide, double radius, UniformGrid& Grid, Inputs& Ini);
int LoadCells(char* fname, Cell* cells, UniformGrid& Grid, double& t, double& dt);
void SaveCells(FILE* FID, Cell* cells, int N_cells, double t, double dt);
Inputs ReadParameters(char* fname);
int GetFileLen(FILE* myFile);
char* GetNextString(char*& buffer);
void Output(FILE* FID, int ID, double t, const Cell& cell, const Tensor T);
void Output(FILE* FID, int ID, double t, const Cell& cell, const DoubleCoord F);


#endif /* INPUTOUTPUT_H_ */

 
