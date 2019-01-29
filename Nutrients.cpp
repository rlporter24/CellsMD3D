#include <float.h>
#include <math.h>
#include "Nutrients.h"
#include "Array.h"
#include "Constants.h"
#include <omp.h>

/*************************************************************************************
 * This solves for the concentration of nutrients (C = carbon, O = oxygen, W = waste) within the colony.
 * Diffusion in the agar is NOT computed, but the agar-colony interface is treated like a boundary condition (C = 1, W = 0, dO/dz = 0).
 * UpdateEnvArray loops through the entire environment array and updates. UpdateEnvironment does the reaction-diffusion calculation for
 * each individual cell in the environment array.
 */

int levelM1x(int x)
{
    return x*2-BoxX/2;
}

int levelP1x(int x)
{
    return x/2+BoxX/4;
}

int levelM1y(int y)
{
    return y*2-BoxY/2;
}

int levelP1y(int y)
{
    return y/2+BoxY/4;
}

int levelM1z(int z)
{
    return z*2+1;
}

int levelP1z(int z)
{
    return (z-1)/2;
}

int UpdateEnvArray(Array3D<LocalEnv>* CurrentColony, Array3D<LocalEnv>* PreviousColony, Array3D<LocalAga>** CurrentAgar, Array3D<LocalAga>** PreviousAgar, Array2D<LocalAga>** CurrentWall, Array2D<LocalAga>** PreviousWall, Array3D<double>& DensityShiftP, Array3D<double>& Density1ShiftP, Array3D<double>& Density2ShiftP, Array2D<double>& WallDensityShiftP, Array2D<double>& WallDensity1ShiftP, Array2D<double>& WallDensity2ShiftP, int minX, int maxX, int minY, int maxY, int maxH, Array2D<double>& Height, Array3D<double>& insideColonyDen)
{
	// loop through nutrient array until steady state
	int count = 0;
	IntCoord positionColony; // position in the colony
	IntCoord positionAgar; // position in the agar
	IntCoord2D positionWall; // position on the wall
	int conv = 1;
    int nbk = 1;
    int zmin, zmax, dir;
    int zi, yi, xi;
    double previous_carbon;

	while (count<maxIter && nbk)
	{
		conv = 0;
        if (count%2==0) {
            zmin = -BoxZAgar;
            zmax = maxH+1;
            dir = 1;
        } else {
            zmin = -maxH-1;
            zmax = BoxZAgar;
            dir = -1;
        }
        
#pragma omp parallel for default(shared) private(zi,yi,xi,positionColony,positionAgar,positionWall,previous_carbon) reduction(+:conv)
        for (zi=zmin;zi<=zmax;zi++)
        {
            if (zi*dir>0)
            {
                for (yi=minY; yi<=maxY; yi++)
                {
                    for (xi=minX; xi<=maxX; xi++)
                    {
                    
                        positionColony.x = xi;
                        positionColony.y = yi;
                        positionColony.z = zi*dir-1;
                        
                        previous_carbon = PreviousColony->Get(positionColony).Carbon;
                        
                        if (DensityShiftP.Get(positionColony) > DBL_EPSILON)
                        {
                            UpdateEnvironment(CurrentColony, PreviousColony, PreviousWall[0], DensityShiftP, Density1ShiftP, Density2ShiftP,WallDensityShiftP,positionColony,insideColonyDen);
                        }
                    
                    
                        // check convergence criteria for this grid point
                        conv += int( fabs(CurrentColony->Get(positionColony).Carbon - previous_carbon) > ConvCrit*fabs(CurrentColony->Get(positionColony).Carbon) );
                    }
                }
            }
            else if (zi==0)
            {
                for (int level=0; level<maxLevels; level++)
                {
                    for (yi=0; yi<BoxY; yi++)
                    {
                        for (xi=0; xi<BoxX; xi++)
                        {
                            positionWall.x=xi;
                            positionWall.y=yi;
                            
                            previous_carbon = PreviousWall[0]->Get(positionWall).CarbonAgar;
                            
                            UpdateWall(CurrentColony, PreviousAgar, CurrentWall, PreviousWall, WallDensityShiftP, WallDensity1ShiftP, WallDensity2ShiftP, positionWall, level, Height, insideColonyDen);
                            
                            // check convergence criteria for this grid point
                            conv += int( fabs(CurrentWall[0]->Get(positionWall).CarbonAgar - previous_carbon) > ConvCrit*fabs(CurrentWall[0]->Get(positionWall).CarbonAgar) );
                        }
                    }
                }
            }
            else
            {
                for (int level=0; level<maxLevels; level++)
                {
                    for (yi=0; yi<BoxY; yi++)
                    {
                        for (xi=0; xi<BoxX; xi++)
                        {
                            positionAgar.x = xi;
                            positionAgar.y = yi;
                            positionAgar.z = -(zi*dir+1);
                            
                            previous_carbon = PreviousAgar[0]->Get(positionAgar).CarbonAgar;
                            
                            UpdateAgar(CurrentAgar, PreviousAgar, CurrentWall, positionAgar, level);
                        
                            // check convergence criteria for this grid point
                            conv += int( fabs(CurrentAgar[0]->Get(positionAgar).CarbonAgar - previous_carbon) > ConvCrit*fabs(CurrentAgar[0]->Get(positionAgar).CarbonAgar) );
                        }
                    }
                }
            }
        }
        

		if (NutrientGSI==0)
		{
			Array3D<LocalEnv>* swapColony = PreviousColony;
			PreviousColony = CurrentColony;
			CurrentColony = swapColony;

			Array3D<LocalAga>** swapAgar = PreviousAgar;
			PreviousAgar = CurrentAgar;
			CurrentAgar = swapAgar;

			Array2D<LocalAga>** swapWall = PreviousWall;
			PreviousWall = CurrentWall;
			CurrentWall = swapWall;
		}

		if ((count>minIter)&&(conv==0))
			nbk=0;

		count++;
        
        
	}
    
	return count;

}

void UpdateEnvironment(Array3D<LocalEnv>* Env, Array3D<LocalEnv>* prevEnv, Array2D<LocalAga>* prevWal, Array3D<double>& DensityShiftP, Array3D<double>& Density1ShiftP,Array3D<double>& Density2ShiftP,Array2D<double>& WallDensityShiftP,IntCoord position, Array3D<double>& insideColonyDen)
{
	IntCoord2D positionwall;
	positionwall.x=position.x;
	positionwall.y=position.y;
	if (insideColonyDen.Get(position)>0.5)//((position.z==0 && WallDensityShiftP.Get(positionwall) > 0.1) || DensityShiftP.Get(position)>0.05)	// if density is zero, do not update
	{
	
		double C = prevEnv->Get(position).Carbon;

		// calculate growth rate and consumption rates based on current concentrations
		double fC = C/(C+KC);	// carbon transporter rate

		double Cgr=0;; // redistribute C to vertices, in order to compute growth rate
		if (position.x!=0 && position.y!=0)
		{
			if (position.z==0)
			{
				Cgr=(prevWal->Get(position.x-1,position.y-1).CarbonAgar+prevEnv->Get(position.x-1,position.y-1,position.z).Carbon
						+prevWal->Get(position.x-1,position.y).CarbonAgar+prevEnv->Get(position.x-1,position.y,position.z).Carbon
						+prevWal->Get(position.x,position.y-1).CarbonAgar+prevEnv->Get(position.x,position.y-1,position.z).Carbon
						+prevWal->Get(position.x,position.y).CarbonAgar+prevEnv->Get(position.x,position.y,position.z).Carbon)/8.0;
			}
			else
			{
				Cgr=(prevEnv->Get(position.x-1,position.y-1,position.z-1).Carbon+prevEnv->Get(position.x-1,position.y-1,position.z).Carbon
						+prevEnv->Get(position.x-1,position.y,position.z-1).Carbon+prevEnv->Get(position.x-1,position.y,position.z).Carbon
						+prevEnv->Get(position.x,position.y-1,position.z-1).Carbon+prevEnv->Get(position.x,position.y-1,position.z).Carbon
						+prevEnv->Get(position.x,position.y,position.z-1).Carbon+prevEnv->Get(position.x,position.y,position.z).Carbon)/8.0;
            }
		}

		double fCgr = Cgr/(Cgr+KC);

		double qC, GR;

		qC = C_rate*fC;	// carbon consumption rate
        GR = maxGrowthRate*max(0.0,fCgr-Maintenance_rate/C_rate);
        
		Env->At(position).GrowthRate = GR;

		// second derivatives of C
		double Cxx, Cyy, Czz;

		MyAssert(position.x>2,"Need more boxes");

		// if we're at the bottom boundary
		if (position.z==0)
		{
			Czz = (prevEnv->Get(position.x,position.y,position.z+1).Carbon + prevWal->Get(position.x,position.y).CarbonAgar - 2.0*prevEnv->Get(position.x,position.y,position.z).Carbon)/(BoxLength*BoxLength);
		}
		else
		{
			Czz = (prevEnv->Get(position.x,position.y,position.z+1).Carbon + prevEnv->Get(position.x,position.y,position.z-1).Carbon - 2.0*prevEnv->Get(position.x,position.y,position.z).Carbon)/(BoxLength*BoxLength);
		}

		Cxx = (prevEnv->Get(position.x+1,position.y,position.z).Carbon + prevEnv->Get(position.x-1,position.y,position.z).Carbon - 2.0*prevEnv->Get(position.x,position.y,position.z).Carbon)/(BoxLength*BoxLength);

		Cyy = (prevEnv->Get(position.x,position.y+1,position.z).Carbon + prevEnv->Get(position.x,position.y-1,position.z).Carbon - 2.0*prevEnv->Get(position.x,position.y,position.z).Carbon)/(BoxLength*BoxLength);

		// Calculate new C concentration
        double Cnew;
        
        Cnew = prevEnv->Get(position).Carbon + (DiffColony*(Cxx + Cyy + Czz ) - qC)*Cdt;//*DensityShiftP.Get(position))*Cdt;

        Env->At(position).Carbon = max(0.0,min(Cnew,maxCarbon));	// for stability
		if (NutrientGSI==1)
		{
			prevEnv->At(position).Carbon = max(0.0,min(Cnew,maxCarbon));	// for stability
		}


    }
    else if (insideColonyDen.Get(position)>0.3 && position.z>0)
    {
        int ix = position.x;
        int iy = position.y;
        int iz = position.z;
        
        double Cnew = (prevEnv->Get(ix+1,iy,iz).Carbon*(insideColonyDen.Get(ix+1,iy,iz)>0.5)+prevEnv->Get(ix-1,iy,iz).Carbon*(insideColonyDen.Get(ix-1,iy,iz)>0.5)+prevEnv->Get(ix,iy+1,iz).Carbon*(insideColonyDen.Get(ix,iy+1,iz)>0.5)+prevEnv->Get(ix,iy-1,iz).Carbon*(insideColonyDen.Get(ix,iy-1,iz)>0.5)+prevEnv->Get(ix,iy,iz+1).Carbon*(insideColonyDen.Get(ix,iy,iz+1)>0.5)+prevEnv->Get(ix,iy,iz-1).Carbon*(insideColonyDen.Get(ix,iy,iz-1)>0.5))/((insideColonyDen.Get(ix+1,iy,iz)>0.5)+(insideColonyDen.Get(ix-1,iy,iz)>0.5)+(insideColonyDen.Get(ix,iy+1,iz)>0.5)+(insideColonyDen.Get(ix,iy-1,iz)>0.5)+(insideColonyDen.Get(ix,iy,iz+1)>0.5)+(insideColonyDen.Get(ix,iy,iz-1)>0.5));
        
        Env->At(position).Carbon = max(0.0,min(Cnew, maxCarbon));
        if (NutrientGSI==1)
        {
            prevEnv->At(position).Carbon = max(0.0,min(Cnew,maxCarbon));	// for stability
        }
    }
    else if (insideColonyDen.Get(position)>0.3 && position.z==0)
    {
        int ix = position.x;
        int iy = position.y;
        int iz = position.z;
        
        double Cnew = (prevEnv->Get(ix+1,iy,iz).Carbon*(insideColonyDen.Get(ix+1,iy,iz)>0.5)+prevEnv->Get(ix-1,iy,iz).Carbon*(insideColonyDen.Get(ix-1,iy,iz)>0.5)+prevEnv->Get(ix,iy+1,iz).Carbon*(insideColonyDen.Get(ix,iy+1,iz)>0.5)+prevEnv->Get(ix,iy-1,iz).Carbon*(insideColonyDen.Get(ix,iy-1,iz)>0.5)+prevEnv->Get(ix,iy,iz+1).Carbon*(insideColonyDen.Get(ix,iy,iz+1)>0.5))/((insideColonyDen.Get(ix+1,iy,iz)>0.5)+(insideColonyDen.Get(ix-1,iy,iz)>0.5)+(insideColonyDen.Get(ix,iy+1,iz)>0.5)+(insideColonyDen.Get(ix,iy-1,iz)>0.5)+(insideColonyDen.Get(ix,iy,iz+1)>0.5));
        
        Env->At(position).Carbon = max(0.0,min(Cnew, maxCarbon));
        if (NutrientGSI==1)
        {
            prevEnv->At(position).Carbon = max(0.0,min(Cnew,maxCarbon));	// for stability
        }
    }
}

void UpdateAgar(Array3D<LocalAga>** Aga, Array3D<LocalAga>** prevAga, Array2D<LocalAga>** Wal, IntCoord position, int level)
{
    
    // second derivatives of C
    double Cxx, Cyy, Czz;
    
    
    // if we're at the upper boundary
    if (position.z==0)
    {
        Czz = (prevAga[level]->Get(position.x,position.y,position.z+1).CarbonAgar + Wal[level]->Get(position.x,position.y).CarbonAgar - 2.0*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
    }
    else if (level>0 && position.y>=BoxY/4 && position.y<3*BoxY/4 && position.x>=BoxX/4 && position.x<3*BoxX/4 && position.z<=BoxZAgar/2)
    {
        if (position.z==BoxZAgar/2)
        {
            Czz = (prevAga[level]->Get(position.x,position.y,position.z+1).CarbonAgar + prevAga[level-1]->Get(levelM1x(position.x),levelM1y(position.y),BoxZAgar-1).CarbonAgar - 2.0*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
    }
    else if (level<maxLevels-1 && position.z==BoxZAgar-1)
    {
        int lx=position.x%2, ly=position.y%2;
        Czz = (prevAga[level]->Get(position.x,position.y,position.z-2).CarbonAgar + 1.0/(lx+1)/(ly+1)*prevAga[level+1]->Get(levelP1x(position.x),levelP1y(position.y),BoxZAgar/2).CarbonAgar + lx*1.0/(lx+1)/(ly+1)*prevAga[level+1]->Get((levelP1x(position.x)+1)%BoxX,levelP1y(position.y),BoxZAgar/2).CarbonAgar + ly*1.0/(lx+1)/(ly+1)*prevAga[level+1]->Get(levelP1x(position.x),(levelP1y(position.y)+1)%BoxY,BoxZAgar/2).CarbonAgar + lx*ly*1.0/(lx+1)/(ly+1)*prevAga[level+1]->Get((levelP1x(position.x)+1)%BoxX,(levelP1y(position.y)+1)%BoxY,BoxZAgar/2).CarbonAgar - 2.0*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/4/(BoxLength*BoxLength*pow(4.0,level));
    }
    else if (level==maxLevels-1 && position.z==BoxZAgar-1)
    {
        Czz = 2.0*(prevAga[level]->Get(position.x,position.y,position.z-1).CarbonAgar-prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
    }
    else
    {
        Czz = (prevAga[level]->Get(position.x,position.y,position.z+1).CarbonAgar + prevAga[level]->Get(position.x,position.y,position.z-1).CarbonAgar - 2.0*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
    }
    
    ////////////////////////////////////////////////
    
    if (level>0 && position.y>=BoxY/4 && position.y<3*BoxY/4 && (position.x>=BoxX/4-1) && position.x<=BoxX*3/4 && position.z<=BoxZAgar/2-1)
    {
        if (position.x==BoxX/4-1)
        {
            Cxx = (prevAga[level]->Get(position.x-1,position.y,position.z).CarbonAgar+prevAga[level-1]->Get(0,levelM1y(position.y),levelM1z(position.z)).CarbonAgar-2*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
        else if (position.x==BoxX*3/4)
        {
            Cxx = (prevAga[level]->Get(position.x+1,position.y,position.z).CarbonAgar+prevAga[level-1]->Get(BoxX-2,levelM1y(position.y),levelM1z(position.z)).CarbonAgar-2*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
    }
    else if (level<maxLevels-1 && (position.x==0 || position.x==BoxX-1))
    {
        int ly=position.y%2, lz=(position.z+1)%2;
        if (position.x==0)
        {
            Cxx = (prevAga[level]->Get(position.x+2,position.y,position.z).CarbonAgar+1.0/(ly+1)/(lz+1)*prevAga[level+1]->Get(BoxX/4-1,levelP1y(position.y),levelP1z(position.z)).CarbonAgar+ly*1.0/(ly+1)/(lz+1)*prevAga[level+1]->Get(BoxX/4-1,(levelP1y(position.y)+1)%BoxY,levelP1z(position.z)).CarbonAgar+lz*1.0/(ly+1)/(lz+1)*prevAga[level+1]->Get(BoxX/4-1,levelP1y(position.y),(levelP1z(position.z)+1)%BoxZAgar).CarbonAgar+ly*lz*1.0/(ly+1)/(lz+1)*prevAga[level+1]->Get(BoxX/4-1,(levelP1y(position.y)+1)%BoxY,(levelP1z(position.z)+1)%BoxZAgar).CarbonAgar-2*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/4/(BoxLength*BoxLength*pow(4.0,level));
        }
        else
        {
            Cxx = (prevAga[level]->Get(position.x-1,position.y,position.z).CarbonAgar+1.0/(ly+1)/(lz+1)*prevAga[level+1]->Get(BoxX*3/4,levelP1y(position.y),levelP1z(position.z)).CarbonAgar+ly*1.0/(ly+1)/(lz+1)*prevAga[level+1]->Get(BoxX*3/4,(levelP1y(position.y)+1)%BoxY,levelP1z(position.z)).CarbonAgar+lz*1.0/(ly+1)/(lz+1)*prevAga[level+1]->Get(BoxX*3/4,levelP1y(position.y),(levelP1z(position.z)+1)%BoxZAgar).CarbonAgar+ly*lz*1.0/(ly+1)/(lz+1)*prevAga[level+1]->Get(BoxX*3/4,(levelP1y(position.y)+1)%BoxY,(levelP1z(position.z)+1)%BoxZAgar).CarbonAgar-2*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
    }
    else if (level==maxLevels-1 && (position.x==0 || position.x==BoxX-1))
    {
        if (position.x==0)
        {
            Cxx = (prevAga[level]->Get(position.x+1,position.y,position.z).CarbonAgar+maxCarbon-2.0*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
        else if (position.x==(BoxX-1))
        {
            Cxx = (maxCarbon+prevAga[level]->Get(position.x-1,position.y,position.z).CarbonAgar-2.0*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
    }
    else
    {
        Cxx = (prevAga[level]->Get(position.x-1,position.y,position.z).CarbonAgar+prevAga[level]->Get(position.x+1,position.y,position.z).CarbonAgar-2*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
    }
    
    /////////////////////////////////////
    
    if (level>0 && position.x>=BoxX/4 && position.x<3*BoxX/4 && (position.y>=BoxY/4-1) && position.y<=BoxY*3/4 && position.z<=BoxZAgar/2-1)
    {
        if (position.y==BoxY/4-1)
        {
            Cyy = (prevAga[level]->Get(position.x,position.y-1,position.z).CarbonAgar+prevAga[level-1]->Get(levelM1x(position.x),0,levelM1z(position.z)).CarbonAgar-2*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
        else if (position.y==BoxY*3/4)
        {
            Cyy = (prevAga[level]->Get(position.x,position.y+1,position.z).CarbonAgar+prevAga[level-1]->Get(levelM1x(position.x),BoxY-2,levelM1z(position.z)).CarbonAgar-2*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
    }
    else if (level<maxLevels-1 && (position.y==0 || position.y==BoxY-1))
    {
        int lx=position.x%2, lz=(position.z+1)%2;
        if (position.y==0)
        {
            Cyy = (prevAga[level]->Get(position.x,position.y+2,position.z).CarbonAgar+1.0/(lx+1)/(lz+1)*prevAga[level+1]->Get(levelP1x(position.x),BoxY/4-1,levelP1z(position.z)).CarbonAgar+lx*1.0/(lx+1)/(lz+1)*prevAga[level+1]->Get((levelP1x(position.x)+1)%BoxX,BoxY/4-1,levelP1z(position.z)).CarbonAgar+lz*1.0/(lx+1)/(lz+1)*prevAga[level+1]->Get(levelP1x(position.x),BoxY/4-1,(levelP1z(position.z)+1)%BoxZAgar).CarbonAgar+lx*lz*1.0/(lx+1)/(lz+1)*prevAga[level+1]->Get((levelP1x(position.x)+1)%BoxX,BoxY/4-1,(levelP1z(position.z)+1)%BoxZAgar).CarbonAgar-2*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/4/(BoxLength*BoxLength*pow(4.0,level));
        }
        else
        {
            Cyy = (prevAga[level]->Get(position.x,position.y-1,position.z).CarbonAgar+1.0/(lx+1)/(lz+1)*prevAga[level+1]->Get(levelP1x(position.x),BoxY*3/4,levelP1z(position.z)).CarbonAgar+lx*1.0/(lx+1)/(lz+1)*prevAga[level+1]->Get((levelP1x(position.x)+1)%BoxX,BoxY*3/4,levelP1z(position.z)).CarbonAgar+lz*1.0/(lx+1)/(lz+1)*prevAga[level+1]->Get(levelP1x(position.x),BoxY*3/4,(levelP1z(position.z)+1)%BoxZAgar).CarbonAgar+lx*lz*1.0/(lx+1)/(lz+1)*prevAga[level+1]->Get((levelP1x(position.x)+1)%BoxX,BoxY*3/4,(levelP1z(position.z)+1)%BoxZAgar).CarbonAgar-2*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
    }
    else if (level==maxLevels-1 && (position.y==0 || position.y==BoxY-1))
    {
        if (position.y==0)
        {
            Cyy = (prevAga[level]->Get(position.x,position.y+1,position.z).CarbonAgar+maxCarbon-2.0*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
        else if (position.y==(BoxY-1))
        {
            Cyy = (maxCarbon+prevAga[level]->Get(position.x,position.y-1,position.z).CarbonAgar-2.0*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
    }
    else
    {
        Cyy = (prevAga[level]->Get(position.x,position.y-1,position.z).CarbonAgar+prevAga[level]->Get(position.x,position.y+1,position.z).CarbonAgar-2*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
    }
    
    ///////////////////////////////////////////
    
    
    
    // Calculate new C concentration
    double Cnew = prevAga[level]->Get(position).CarbonAgar + DiffAgar*(Cxx + Cyy + Czz)*Cdt;
    Aga[level]->At(position).CarbonAgar = max(0.0,min(Cnew,maxCarbon));	// for stability
    if (NutrientGSI==1)
    {
        prevAga[level]->At(position).CarbonAgar = max(0.0,min(Cnew,maxCarbon));	// for stability
    }
    
}

void UpdateWall(Array3D<LocalEnv>* Env, Array3D<LocalAga>** prevAga, Array2D<LocalAga>** Wal, Array2D<LocalAga>** prevWal, Array2D<double>& WallDensityShiftP, Array2D<double>& WallDensity1ShiftP, Array2D<double>& WallDensity2ShiftP, IntCoord2D position, int level, Array2D<double>& Height, Array3D<double>& insideColonyDen)
{
    //if ((position.x>=minX) && (position.x<=maxX) && (position.y>=minY) && (position.y<=maxY))
    IntCoord positionWall3D;
    positionWall3D.x = position.x;
    positionWall3D.y = position.y;
    positionWall3D.z = 0;
    if (level==0 && insideColonyDen.Get(positionWall3D)>0.5) //(level==0 && WallDensityShiftP.Get(position) > 0.8 && Height.At(position) > cellRadius+1e-7)
    {
        double C = prevWal[level]->Get(position).CarbonAgar;
        
        // calculate growth rate and consumption rates based on current concentrations
        double fC = C/(C+KC);	// carbon transporter rate
        
        double qC;
        
        qC = C_rate*fC;	// carbon consumption rate
        
        switch (InterfaceCondition)
        {
            case 1:
            {
                double Cnew = (DiffColony*Env->Get(positionWall3D).Carbon + DiffAgar*prevAga[level]->Get(positionWall3D).CarbonAgar)/(DiffColony+DiffAgar);
                Wal[level]->At(position).CarbonAgar = max(0.0,min(Cnew,maxCarbon)); // for stability
                if (NutrientGSI==1)
                {
                    prevWal[level]->At(position).CarbonAgar = max(0.0,min(Cnew,maxCarbon)); // for stability
                }
                break;
            }
            case 2:
            {
                double Cnew = (DiffColony*Env->Get(positionWall3D).Carbon + DiffAgar*prevAga[level]->Get(positionWall3D).CarbonAgar-qC*BoxLength*BoxLength)/(DiffColony+DiffAgar);
                Wal[level]->At(position).CarbonAgar = max(0.0,min(Cnew,maxCarbon)); // for stability
                if (NutrientGSI==1)
                {
                    prevWal[level]->At(position).CarbonAgar = max(0.0,min(Cnew,maxCarbon)); // for stability
                }
                break;
            }
            case 3:
            {
                double Chalf=(prevWal[level]->Get(position).CarbonAgar+Env->Get(position.x,position.y,0).Carbon)/2.0;
                double qChalf=C_rate*Chalf/(Chalf+KC);
                
                double Cnew = (DiffColony*Env->Get(positionWall3D).Carbon + DiffAgar*prevAga[level]->Get(positionWall3D).CarbonAgar-qChalf*BoxLength*BoxLength/2.0)/(DiffColony+DiffAgar);
                Wal[level]->At(position).CarbonAgar = max(0.0,min(Cnew,maxCarbon)); // for stability
                if (NutrientGSI==1)
                {
                    prevWal[level]->At(position).CarbonAgar = max(0.0,min(Cnew,maxCarbon)); // for stability
                }
                break;
            }
            default:
            {
                printf("Wrong switch on interface condition!");
            }
        }
    }
    else
    {
        double Cxx, Cyy, Czz;
        if (level>0 && position.y>=BoxY/4 && position.y<3*BoxY/4 && (position.x>=BoxX/4-1) && position.x<=BoxX*3/4)
        {
            if (position.x==BoxX/4-1)
            {
                Cxx = (prevWal[level]->Get(position.x-1,position.y).CarbonAgar+prevWal[level-1]->Get(0,levelM1y(position.y)).CarbonAgar-2*prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
            }
            else if (position.x==BoxX*3/4)
            {
                Cxx = (prevWal[level]->Get(position.x+1,position.y).CarbonAgar+prevWal[level-1]->Get(BoxX-2,levelM1y(position.y)).CarbonAgar-2*prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
            }
        }
        else if (level<maxLevels-1 && (position.x==0 || position.x==BoxX-1))
        {
            int ly=position.y%2;
            if (position.x==0)
            {
                Cxx = (prevWal[level]->Get(position.x+2,position.y).CarbonAgar+1.0/(ly+1)*prevWal[level+1]->Get(BoxX/4-1,levelP1y(position.y)).CarbonAgar+ly*1.0/(ly+1)*prevWal[level+1]->Get(BoxX/4-1,(levelP1y(position.y)+1)%BoxY).CarbonAgar-2*prevWal[level]->Get(position.x,position.y).CarbonAgar)/4/(BoxLength*BoxLength*pow(4.0,level));
            }
            else
            {
                Cxx = (prevWal[level]->Get(position.x-1,position.y).CarbonAgar+1.0/(ly+1)*prevWal[level+1]->Get(BoxX*3/4,levelP1y(position.y)).CarbonAgar+ly*1.0/(ly+1)*prevWal[level+1]->Get(BoxX*3/4,(levelP1y(position.y)+1)%BoxY).CarbonAgar-2*prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
            }
        }
        else if (level==maxLevels-1 && (position.x==0 || position.x==BoxX-1))
        {
            if (position.x==0)
            {
                Cxx = (prevWal[level]->Get(position.x+1,position.y).CarbonAgar+maxCarbon-2.0*prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
            }
            else if (position.x==(BoxX-1))
            {
                Cxx = (maxCarbon+prevWal[level]->Get(position.x-1,position.y).CarbonAgar-2.0*prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
            }
        }
        else
        {
            Cxx = (prevWal[level]->Get(position.x-1,position.y).CarbonAgar+prevWal[level]->Get(position.x+1,position.y).CarbonAgar-2*prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
        
        /////////////////////////////////////
        
        if (level>0 && position.x>=BoxX/4 && position.x<3*BoxX/4 && (position.y>=BoxY/4-1) && position.y<=BoxY*3/4)
        {
            if (position.y==BoxY/4-1)
            {
                Cyy = (prevWal[level]->Get(position.x,position.y-1).CarbonAgar+prevWal[level-1]->Get(levelM1x(position.x),0).CarbonAgar-2*prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
            }
            else if (position.y==BoxY*3/4)
            {
                Cyy = (prevWal[level]->Get(position.x,position.y+1).CarbonAgar+prevWal[level-1]->Get(levelM1x(position.x),BoxY-2).CarbonAgar-2*prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
            }
        }
        else if (level<maxLevels-1 && (position.y==0 || position.y==BoxX-1))
        {
            int lx=position.x%2;
            if (position.y==0)
            {
                Cyy = (prevWal[level]->Get(position.x,position.y+2).CarbonAgar+1.0/(lx+1)*prevWal[level+1]->Get(levelP1x(position.x),BoxY/4-1).CarbonAgar+lx*1.0/(lx+1)*prevWal[level+1]->Get((levelP1x(position.x)+1)%BoxX,BoxY/4-1).CarbonAgar-2*prevWal[level]->Get(position.x,position.y).CarbonAgar)/4/(BoxLength*BoxLength*pow(4.0,level));
            }
            else
            {
                Cyy = (prevWal[level]->Get(position.x,position.y-1).CarbonAgar+1.0/(lx+1)*prevWal[level+1]->Get(levelP1x(position.x),BoxY*3/4).CarbonAgar+lx*1.0/(lx+1)*prevWal[level+1]->Get((levelP1x(position.x)+1)%BoxX,BoxY*3/4).CarbonAgar-2*prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
            }
        }
        else if (level==maxLevels-1 && (position.y==0 || position.y==BoxY-1))
        {
            if (position.y==0)
            {
                Cyy = (prevWal[level]->Get(position.x,position.y+1).CarbonAgar+maxCarbon-2.0*prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
            }
            else if (position.y==(BoxY-1))
            {
                Cyy = (maxCarbon+prevWal[level]->Get(position.x,position.y-1).CarbonAgar-2.0*prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
            }
        }
        else
        {
            Cyy = (prevWal[level]->Get(position.x,position.y-1).CarbonAgar+prevWal[level]->Get(position.x,position.y+1).CarbonAgar-2*prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
        
        ///////////////////////////////////////////
        
        
        Czz = 2.0*(prevAga[level]->Get(position.x,position.y,0).CarbonAgar-prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
        
        double Cnew = prevWal[level]->Get(position).CarbonAgar + DiffAgar*(Cxx + Cyy+ Czz)*Cdt;
        Wal[level]->At(position).CarbonAgar = max(0.0,min(Cnew,maxCarbon));	// for stability
        if (NutrientGSI==1)
        {
            prevWal[level]->At(position).CarbonAgar = max(0.0,min(Cnew,maxCarbon));	// for stability
        }
        
        
    }
    
}

