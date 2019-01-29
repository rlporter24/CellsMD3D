#ifndef TOOLS_H_
#define TOOLS_H_
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

/******************************************************************************
 * Contains helper functions for doing computations on arrays and segments
 ******************************************************************************/

const double PI = 3.141592653589793238462643383279502884197;

#ifdef NDEBUG
inline void MyAssert(bool b, const char* msg){}

#else
inline void MyAssert(bool b, const char* msg)
{
	if (!b)
	{
		printf("%s\n", msg);
		fflush(stdout);
		assert(false);
		exit(-1);
	}
}
#endif

template<typename T>
struct Coord2D
{
	Coord2D(){}
	Coord2D(T _x, T _y): x(_x), y(_y) {}

	T x;
	T y;
};

template<typename T>
struct Coord3D
{
	Coord3D(){}
	Coord3D(T _x, T _y, T _z): x(_x), y(_y), z(_z) {}

	T x;
	T y;
	T z;
};

typedef Coord3D<double> DoubleCoord;
typedef Coord3D<int> IntCoord;
typedef Coord2D<double> DoubleCoord2D;
typedef Coord2D<int> IntCoord2D;


// 3D symmetric tensor for storing the forces on cells
struct Tensor
{
	Tensor(){}
	Tensor(double _xx, double _yy, double _zz, double _rr, double _tt): xx(_xx), yy(_yy), zz(_zz), rr(_rr), tt(_tt) {}

	double xx;
	double yy;
	double zz;
	double rr;
	double tt;
};

// stores two coordinates for the positions of the cell vertices
struct Segment
{
	Segment(){}
	
	Segment(const DoubleCoord& _p, const DoubleCoord& _q)
		: p(_p), q(_q) {}

	Segment(double px, double py, double pz, double qx, double qy, double qz)
		: p(px, py, pz), q(qx, qy, pz) {}

	DoubleCoord p;	// first point
	DoubleCoord q;	// second point
    double time_p;
    double time_q;
    int age_p;
    int age_q;
};

inline double max(double n1, double n2)
{
	return n1>n2? n1: n2;
}

inline int max(int n1, int n2)
{
	return n1>n2? n1: n2;
}

inline double min(double n1, double n2)
{
	return n1<n2? n1: n2;
}

inline double dot(const DoubleCoord& v1, const DoubleCoord& v2)
{
	return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
}

inline DoubleCoord cross(const DoubleCoord& v1, const DoubleCoord& v2)
{
	return DoubleCoord(v1.y*v2.z-v1.z*v2.y, -v1.x*v2.z+v1.z*v2.x, v1.x*v2.y-v1.y*v2.x);
}

inline DoubleCoord diff(const DoubleCoord& a, const DoubleCoord& b)
{
	return DoubleCoord(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline DoubleCoord sum(const DoubleCoord& a, const DoubleCoord& b)
{
	return DoubleCoord(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline DoubleCoord scale(const DoubleCoord& a, double C)
{
	return DoubleCoord(a.x*C, a.y*C, a.z*C);
}

inline Segment sum(const Segment& a, const Segment& b)
{
	return Segment(sum(a.p, b.p), sum(a.q, b.q));
}

inline Segment diff(const Segment& a, const Segment& b)
{
	return Segment(diff(a.p, b.p), diff(a.q, b.q));
}

inline Segment scale(const Segment& a, double C)
{
	return Segment(scale(a.p, C), scale(a.q, C));
}

inline DoubleCoord average(const Segment& a)
{
	return scale(sum(a.p,a.q),0.5);
}

inline double linear_interp(const double y1, const double y2, const double dx)
{
	double slope = (y2-y1);
	if (dx>0)
		return y1+dx*slope;
	else
		return y2+dx*slope;
}

#endif /* TOOLS_H_ */

