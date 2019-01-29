#ifndef NUTRIENTS_H_
#define NUTRIENTS_H_

#include "tools.h"
#include "Constants.h"
//#include "Array.h"

template<typename T>
class Array3D;

struct LocalEnv
{
	LocalEnv()
	{
		Carbon = 0.0;
		GrowthRate = 0.0;
	}

	LocalEnv(double _carbon, double _growthrate): Carbon(_carbon), GrowthRate(_growthrate) {}

	double Carbon;
	double GrowthRate;
};

struct LocalAga
{
	LocalAga()
	{
		CarbonAgar = maxCarbon;
	}

	LocalAga(double _carbonAgar, double _growthrateAgar): CarbonAgar(_carbonAgar) {}

	double CarbonAgar;
};

template<typename T>
class Array2D;





int UpdateEnvArray(Array3D<LocalEnv>* CurrentColony, Array3D<LocalEnv>* PreviousColony, Array3D<LocalAga>** CurrentAgar, Array3D<LocalAga>** PreviousAgar, Array2D<LocalAga>** CurrentWall, Array2D<LocalAga>** PreviousWall, Array3D<double>& DensityShiftP, Array3D<double>& Density1ShiftP, Array3D<double>& Density2ShiftP, Array2D<double>& WallDensityShiftP,Array2D<double>& WallDensity1ShiftP,Array2D<double>& WallDensity2ShiftP, int minX, int maxX, int minY, int maxY, int maxH, Array2D<double>& Height, Array3D<double>& insideColonyDen);
void UpdateEnvironment(Array3D<LocalEnv>* Env, Array3D<LocalEnv>* prevEnv, Array2D<LocalAga>* prevWal, Array3D<double>& DensityShiftP, Array3D<double>& Density1ShiftP, Array3D<double>& Density2ShiftP, Array2D<double>& WallDensityShiftP,IntCoord position, Array3D<double>& insideColonyDen);
void UpdateAgar(Array3D<LocalAga>** Aga, Array3D<LocalAga>** prevAga, Array2D<LocalAga>** Wal, IntCoord position, int level);
void UpdateWall(Array3D<LocalEnv>* Env, Array3D<LocalAga>** prevAga, Array2D<LocalAga>** Wal, Array2D<LocalAga>** prevWal, Array2D<double>& WallDensityShiftP, Array2D<double>& WallDensity1ShiftP, Array2D<double>& WallDensity2ShiftP, IntCoord2D position, int level, Array2D<double>& Height, Array3D<double>& insideColonyDen);
#endif
