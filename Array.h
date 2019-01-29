#ifndef ARRAY_H_
#define ARRAY_H_
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "tools.h"
#include "Nutrients.h"


/*************************************
 * A bunch of functions for array manipulations
 ************************************/
template<typename T>
class Array2D
{
public:
	IntCoord2D m_Size;

	Array2D()
	{
		Data = NULL;
	}

	Array2D(int xSize, int ySize)
			:m_Size(xSize, ySize)
	{
		Data = new T [m_Size.x*m_Size.y];	// 1D array of doubles
	}

	const IntCoord2D& Size() { return m_Size; }

	~Array2D(void)
	{
		if (Data!=NULL)
			delete[] Data;
	}

	T Get(const IntCoord2D& coord)
	{
		return Get(coord.x, coord.y);
	}

	T Get(int xi, int yi)
	{
		return Data[Index(xi, yi)];
	}

	T& At(const IntCoord2D& coord)
	{
		return At(coord.x, coord.y);
	}

	T& At(int xi, int yi)
	{
		return Data[Index(xi, yi)];
	}

	const T& At(const IntCoord2D& coord) const
	{
		return At(coord.x, coord.y);
	}

	const T& At(int xi, int yi) const
	{
		return Data[Index(xi, yi)];
	}

	void Set(int xi, int yi, T value)
	{
		const int index = Index(xi, yi);
		MyAssert(m_Size.x*m_Size.y>index,"Invalid index into array");
		Data[index] = value;
	}

	void SetData(T* DataArray, int Length)
	{
		MyAssert(m_Size.x*m_Size.y==Length,"Incorrect array size");
		memcpy(Data,DataArray,Length*sizeof(T));
	}

	void Initialize(T val)
	{
		for (int xi = 0; xi < m_Size.x; xi++)
		{
			for (int yi = 0; yi < m_Size.y; yi++)
			{
				Data[xi*m_Size.y+yi] = val;
			}
		}
	}

	void Output(char* Mname)//FILE* FID
	{
		FILE* FID = fopen(Mname,"w");
		Output(FID);
		fclose(FID);
	}

	void Output(FILE* FID);
    
    void Append(FILE* FID);

	double linear_interp(const int x0, const int y0, const double dx, const double dy);

private:
	T* Data;

	int Index(const IntCoord2D& coord) { return Index(coord.x, coord.y); }
	int Index(const int x, const int y) { return x*m_Size.y+y; }
	int Index(const int x, const int y) const { return x*m_Size.y+y; }
};

template<typename T>
class Array3D // think about ordering
{
public:
	IntCoord m_Size;

	Array3D()
	{
		Data = NULL;
	}

	Array3D(int xSize, int ySize, int zSize)
			:m_Size(xSize, ySize, zSize)
	{
		Data = new T [m_Size.x*m_Size.y*m_Size.z];
	}

	const IntCoord& Size() { return m_Size; }

	~Array3D(void)
	{
		if (Data!=NULL)
			delete[] Data;
	}

	T Get(const IntCoord& coord)
	{
		return Get(coord.x, coord.y, coord.z);
	}

	T Get(int xi, int yi, int zi)
	{
		return Data[Index(xi, yi, zi)];
	}

	T& At(const IntCoord& coord)
	{
		return At(coord.x, coord.y, coord.z);
	}

	T& At(int xi, int yi, int zi)
	{
		return Data[Index(xi, yi, zi)];
	}

	const T& At(const IntCoord& coord) const
	{
		return At(coord.x, coord.y, coord.z);
	}

	const T& At(int xi, int yi, int zi) const
	{
		return Data[Index(xi, yi, zi)];
	}

	void Set(int xi, int yi, int zi, T value)
	{
		const int index = Index(xi, yi, zi);
		MyAssert(m_Size.x*m_Size.y*m_Size.z>index,"Invalid index into array");
		Data[index] = value;
	}

	void SetData(T* DataArray, int Length)
	{
		MyAssert(m_Size.x * m_Size.y * m_Size.z == Length, "Incorrect array size");
		memcpy( Data, DataArray, Length*sizeof(T) );
	}

	void Initialize(T val)
	{
		for (int ii = 0; ii < m_Size.x*m_Size.y*m_Size.z; ii++)
		{
			Data[ii] = val;
		}
	}

	void Output(char* Mname)//FILE* FID
	{
		FILE* FID = fopen(Mname,"w");
		Output(FID);
		fclose(FID);
	}

	void Output(FILE* FID);
    
    void Output(FILE* FID, int n);
    
    void Append(FILE* FID, int n);


private:
	T* Data;

	int Index(const IntCoord& coord) { return Index(coord.x, coord.y, coord.z); }
	int Index(const int x, const int y, const int z) {	return x*(m_Size.y*m_Size.z) + y*m_Size.z + z; }
	int Index(const int x, const int y, const int z) const {	return x*(m_Size.y*m_Size.z) + y*m_Size.z + z; }
};

template <>
inline void Array2D<double>::Output(FILE* FID)
{
	rewind(FID);
	for (int yi = 0; yi < m_Size.y; yi++)
	{
		for (int xi = 0; xi < m_Size.x; xi++)
		{
			fprintf(FID, "%6E\t", Data[xi*m_Size.y + yi]);
		}
		fprintf(FID, "\n");
	}
	fflush(FID);
};

template <>
inline void Array2D<int>::Output(FILE* FID)
{
	rewind(FID);
	for (int yi = 0; yi < m_Size.y; yi++)
	{
		for (int xi = 0; xi < m_Size.x; xi++)
		{
			fprintf(FID, "%d\t", Data[xi*m_Size.y + yi]);
		}
		fprintf(FID, "\n");
	}
	fflush(FID);
};

template <>
inline double Array2D<double>::linear_interp(const int x0, const int y0, const double dx, const double dy)
{
	float x = dx<0? 1+dx: dx;
	float y = dy<0? 1+dy: dy;

	return Data[Index(x0, y0)]*(1-x)*(1-y) + Data[Index(x0+1, y0)]*x*(1-y) + Data[Index(x0, y0+1)]*(1-x)*y + Data[Index(x0+1, y0+1)]*x*y;
};

template <>
inline void Array2D<DoubleCoord>::Output(FILE* FID)
{
	rewind(FID);
	for (int yi = 0; yi < m_Size.y; yi++)
	{
		for (int xi = 0; xi < m_Size.x; xi++)
		{
			fprintf(FID, "%6E\t", Data[xi*m_Size.y + yi].x);
		}
		fprintf(FID, "\n");
	}
	for (int yi = 0; yi < m_Size.y; yi++)
	{
		for (int xi = 0; xi < m_Size.x; xi++)
		{
			fprintf(FID, "%6E\t", Data[xi*m_Size.y + yi].y);
		}
		fprintf(FID, "\n");
	}
	for (int yi = 0; yi < m_Size.y; yi++)
	{
		for (int xi = 0; xi < m_Size.x; xi++)
		{
			fprintf(FID, "%6E\t", Data[xi*m_Size.y + yi].z);
		}
		fprintf(FID, "\n");
	}
	fflush(FID);
};

template <>
inline void Array2D<LocalEnv>::Output(FILE* FID)
{
	rewind(FID);
	for (int yi = 0; yi < m_Size.y; yi++)
	{
		for (int xi = 0; xi < m_Size.x; xi++)
		{
			fprintf(FID, "%6E\t", Data[xi*m_Size.y + yi].Carbon);
		}
		fprintf(FID, "\n");
	}

	for (int yi = 0; yi < m_Size.y; yi++)
	{
		for (int xi = 0; xi < m_Size.x; xi++)
		{
			fprintf(FID, "%6E\t", Data[xi*m_Size.y + yi].GrowthRate);
		}
		fprintf(FID, "\n");
	}
	fflush(FID);
};

template <>
inline void Array2D<LocalAga>::Output(FILE* FID)
{
	rewind(FID);
	for (int yi = 0; yi < m_Size.y; yi++)
	{
		for (int xi = 0; xi < m_Size.x; xi++)
		{
			fprintf(FID, "%6E\t", Data[xi*m_Size.y + yi].CarbonAgar);
		}
		fprintf(FID, "\n");
	}
	fflush(FID);
};

template <>
inline void Array2D<LocalAga>::Append(FILE* FID)
{
    for (int yi = 0; yi < m_Size.y; yi++)
    {
        for (int xi = 0; xi < m_Size.x; xi++)
        {
            fprintf(FID, "%6E\t", Data[xi*m_Size.y + yi].CarbonAgar);
        }
        fprintf(FID, "\n");
    }
    
    fflush(FID);
};

template <>
inline void Array3D<double>::Output(FILE* FID)
{
	rewind(FID);
	for (int zi = 0; zi < m_Size.z; zi++)
	{
		for (int yi = 0; yi < m_Size.y; yi++)
		{
			for (int xi = 0; xi < m_Size.x; xi++)
			{
				fprintf(FID, "%6E\t", Data[xi*(m_Size.y*m_Size.z) + yi*m_Size.z + zi]);
			}
			fprintf(FID, "\n");
		}
		fprintf(FID, "\n");
	}
	fflush(FID);
};

template <>
inline void Array3D<double>::Output(FILE* FID, int n)
{
    rewind(FID);
    for (int zi = 0; zi < m_Size.z; zi++)
    {
        for (int yi = 0; yi < m_Size.y; yi+=n)
        {
            for (int xi = 0; xi < m_Size.x; xi+=n)
            {
                fprintf(FID, "%6E\t", Data[xi*(m_Size.y*m_Size.z) + yi*m_Size.z + zi]);
            }
            fprintf(FID, "\n");
        }
        fprintf(FID, "\n");
    }
    fflush(FID);
};


template <>
inline void Array3D<int>::Output(FILE* FID)
{
	rewind(FID);
	for (int zi = 0; zi < m_Size.z; zi++)
	{
		for (int yi = 0; yi < m_Size.y; yi++)
		{
			for (int xi = 0; xi < m_Size.x; xi++)
			{
				fprintf(FID, "%d\t", Data[xi*(m_Size.y*m_Size.z) + yi*m_Size.z + zi]);
			}
			fprintf(FID, "\n");
		}
		fprintf(FID, "\n");
	}
	fflush(FID);
};

template <>
inline void Array3D<DoubleCoord>::Output(FILE* FID)
{
	rewind(FID);
	for (int zi = 0; zi < m_Size.z; zi++)
	{
		for (int yi = 0; yi < m_Size.y; yi++)
		{
			for (int xi = 0; xi < m_Size.x; xi++)
			{
				fprintf(FID, "%6E\t", Data[xi*(m_Size.y*m_Size.z) + yi*m_Size.z + zi].x);
			}
			fprintf(FID, "\n");
		}
		fprintf(FID, "\n");
	}
	for (int zi = 0; zi < m_Size.z; zi++)
	{
		for (int yi = 0; yi < m_Size.y; yi++)
		{
			for (int xi = 0; xi < m_Size.x; xi++)
			{
				fprintf(FID, "%6E\t", Data[xi*(m_Size.y*m_Size.z) + yi*m_Size.z + zi].y);
			}
			fprintf(FID, "\n");
		}
		fprintf(FID, "\n");
	}
	for (int zi = 0; zi < m_Size.z; zi++)
	{
		for (int yi = 0; yi < m_Size.y; yi++)
		{
			for (int xi = 0; xi < m_Size.x; xi++)
			{
				fprintf(FID, "%6E\t", Data[xi*(m_Size.y*m_Size.z) + yi*m_Size.z + zi].z);
			}
			fprintf(FID, "\n");
		}
		fprintf(FID, "\n");
	}
	fflush(FID);
};

template <>
inline void Array3D<DoubleCoord>::Output(FILE* FID, int n)
{
    rewind(FID);
    for (int zi = 0; zi < m_Size.z; zi++)
    {
        for (int yi = 0; yi < m_Size.y; yi+=n)
        {
            for (int xi = 0; xi < m_Size.x; xi+=n)
            {
                fprintf(FID, "%6E\t", Data[xi*(m_Size.y*m_Size.z) + yi*m_Size.z + zi].x);
            }
            fprintf(FID, "\n");
        }
        fprintf(FID, "\n");
    }
    for (int zi = 0; zi < m_Size.z; zi++)
    {
        for (int yi = 0; yi < m_Size.y; yi+=n)
        {
            for (int xi = 0; xi < m_Size.x; xi+=n)
            {
                fprintf(FID, "%6E\t", Data[xi*(m_Size.y*m_Size.z) + yi*m_Size.z + zi].y);
            }
            fprintf(FID, "\n");
        }
        fprintf(FID, "\n");
    }
    for (int zi = 0; zi < m_Size.z; zi++)
    {
        for (int yi = 0; yi < m_Size.y; yi+=n)
        {
            for (int xi = 0; xi < m_Size.x; xi+=n)
            {
                fprintf(FID, "%6E\t", Data[xi*(m_Size.y*m_Size.z) + yi*m_Size.z + zi].z);
            }
            fprintf(FID, "\n");
        }
        fprintf(FID, "\n");
    }
    fflush(FID);
};

template <>
inline void Array3D<IntCoord>::Output(FILE* FID)
{
	rewind(FID);
	for (int zi = 0; zi < m_Size.z; zi++)
	{
		for (int yi = 0; yi < m_Size.y; yi++)
		{
			for (int xi = 0; xi < m_Size.x; xi++)
			{
				fprintf(FID, "%d\t", Data[xi*(m_Size.y*m_Size.z) + yi*m_Size.z + zi].x);
			}
			fprintf(FID, "\n");
		}
		fprintf(FID, "\n");
	}
	for (int zi = 0; zi < m_Size.z; zi++)
	{
		for (int yi = 0; yi < m_Size.y; yi++)
		{
			for (int xi = 0; xi < m_Size.x; xi++)
			{
				fprintf(FID, "%d\t", Data[xi*(m_Size.y*m_Size.z) + yi*m_Size.z + zi].y);
			}
			fprintf(FID, "\n");
		}
		fprintf(FID, "\n");
	}
	for (int zi = 0; zi < m_Size.z; zi++)
	{
		for (int yi = 0; yi < m_Size.y; yi++)
		{
			for (int xi = 0; xi < m_Size.x; xi++)
			{
				fprintf(FID, "%d\t", Data[xi*(m_Size.y*m_Size.z) + yi*m_Size.z + zi].z);
			}
			fprintf(FID, "\n");
		}
		fprintf(FID, "\n");
	}
	fflush(FID);
};

template <>
inline void Array3D<LocalEnv>::Output(FILE* FID)
{
	rewind(FID);
	for (int zi = 0; zi < m_Size.z; zi++)
	{
		for (int yi = 0; yi < m_Size.y; yi++)
		{
			for (int xi = 0; xi < m_Size.x; xi++)
			{
				fprintf(FID, "%6E\t", Data[xi*(m_Size.y*m_Size.z) + yi*m_Size.z + zi].Carbon);
			}
			fprintf(FID, "\n");
		}
	}
	for (int zi = 0; zi < m_Size.z; zi++)
	{
		for (int yi = 0; yi < m_Size.y; yi++)
		{
			for (int xi = 0; xi < m_Size.x; xi++)
			{
				fprintf(FID, "%6E\t", Data[xi*(m_Size.y*m_Size.z) + yi*m_Size.z + zi].GrowthRate);
			}
			fprintf(FID, "\n");
		}
	}
	fflush(FID);
};

template <>
inline void Array3D<LocalEnv>::Output(FILE* FID, int n)
{
    rewind(FID);
    for (int zi = 0; zi < m_Size.z; zi++)
    {
        for (int yi = 0; yi < m_Size.y; yi+=n)
        {
            for (int xi = 0; xi < m_Size.x; xi+=n)
            {
                fprintf(FID, "%6E\t", Data[xi*(m_Size.y*m_Size.z) + yi*m_Size.z + zi].Carbon);
            }
            fprintf(FID, "\n");
        }
    }
    for (int zi = 0; zi < m_Size.z; zi++)
    {
        for (int yi = 0; yi < m_Size.y; yi+=n)
        {
            for (int xi = 0; xi < m_Size.x; xi+=n)
            {
                fprintf(FID, "%6E\t", Data[xi*(m_Size.y*m_Size.z) + yi*m_Size.z + zi].GrowthRate);
            }
            fprintf(FID, "\n");
        }
    }
    fflush(FID);
};

template <>
inline void Array3D<LocalAga>::Output(FILE* FID)
{
	rewind(FID);
	for (int zi = 0; zi < m_Size.z; zi++)
	{
		for (int yi = 0; yi < m_Size.y; yi++)
		{
			for (int xi = 0; xi < m_Size.x; xi++)
			{
				fprintf(FID, "%6E\t", Data[xi*(m_Size.y*m_Size.z) + yi*m_Size.z + zi].CarbonAgar);
			}
			fprintf(FID, "\n");
		}
	}
	fflush(FID);
};

template <>
inline void Array3D<LocalAga>::Output(FILE* FID, int n)
{
    rewind(FID);
    for (int zi = 0; zi < m_Size.z; zi++)
    {
        for (int yi = 0; yi < m_Size.y; yi+=n)
        {
            for (int xi = 0; xi < m_Size.x; xi+=n)
            {
                fprintf(FID, "%6E\t", Data[xi*(m_Size.y*m_Size.z) + yi*m_Size.z + zi].CarbonAgar);
            }
            fprintf(FID, "\n");
        }
    }
    fflush(FID);
};

template <>
inline void Array3D<LocalAga>::Append(FILE* FID, int n)
{
    for (int zi = 0; zi < m_Size.z; zi++)
    {
        for (int yi = 0; yi < m_Size.y; yi+=n)
        {
            for (int xi = 0; xi < m_Size.x; xi+=n)
            {
                fprintf(FID, "%6E\t", Data[xi*(m_Size.y*m_Size.z) + yi*m_Size.z + zi].CarbonAgar);
            }
            fprintf(FID, "\n");
        }
    }
    fflush(FID);
};

typedef Array2D<double> DoubleArray2D;
typedef Array2D<DoubleCoord> CoordArray2D;
typedef Array2D<LocalEnv> EnvArray2D;
typedef Array2D<LocalAga> AgaArray2D;

typedef Array3D<double> DoubleArray3D;
typedef Array3D<DoubleCoord> CoordArray3D;
typedef Array3D<LocalEnv> EnvArray3D;
typedef Array3D<LocalAga> AgaArray3D;



#endif /* ARRAY_H_ */
