#ifndef _UNIFORM_GRID_H_
#define _UNIFORM_GRID_H_

#include <math.h>
#include <iostream>

#include "tools.h"

// Class to store the indices of cells that live within each box of a grid. Important for speeding up calculation of neighbours,
// interfacing between continuum cell positions and discrete fields such as nutrients and densities
class UniformGrid
{
public:
	// address gives the index into the array of boxes
	// a = y*xSize+x, where x and y are the indices in the x and y directions
	struct Address
	{
		Address(){};
		Address(int _a):a(_a){}
		int a;
	};
	
private:
	// A box is one grid location. It can contain maxCells. Actually contains Count pointers to cells that live within it.
	struct Box
	{
		Box(int _maxCells)
		: Count(0), maxCells(_maxCells)	// when it is first created, it has no cells
		{
			Cells = new int[maxCells];
		}
		
		void Add(int cellID)
		{
			MyAssert(Count+1<maxCells,"Tried to add too many cells");
			
			// add cell to first unoccupied position in the array
			Cells[Count] = cellID;
			Count++;
		}
		
		void Remove(int cellID)
		{
			// first find the index of the cell to remove from the array
			int ii = 0;
			while ((Cells[ii]!=cellID)&(ii<Count))
				ii++;
			
			MyAssert(ii!=Count,"Tried to remove a cell that does not exist");
			
			// move the cell from the end of the array into the newly vacant position
			Cells[ii] = Cells[Count-1];
			Count--;
			
		}
		
		~Box() { delete[] Cells; }
		
		int Count;	// the number of cells in the box
		int* Cells;	// an array of cell IDs for cells in the box;
		int maxCells;
	};
	
public:
	// create a grid of boxes of size xSize by ySize, with each box having space for maxCellsPerBox, and
	// having physical dimension of BoxLength
	UniformGrid(int xSize, int ySize, int zSize, int maxCellsPerBox, double BoxLength)
		:m_Size(xSize, ySize, zSize), m_MaxCellsPerBox(maxCellsPerBox), m_BoxLength(BoxLength)
	{
		m_Boxes = new Box*[xSize*ySize*zSize];	// m_Boxes is an array of box pointers
		
		for(int z = 0; z<zSize; z++)
		{
			for(int x=0; x<xSize; x++)
			{
				for(int y=0; y<ySize; y++)
				{
					m_Boxes[x*(ySize*zSize) + y*zSize + z] = new Box(m_MaxCellsPerBox);	// create individual boxes
				}
			}
		}
	}

	int MaxCellsPerBox() { return m_MaxCellsPerBox; }

	// gets the address of the box where a cell should live based on the position of the center of mass
	Address GetAddress(DoubleCoord position)
	{
		// get x and y indexes
		int x_index = int(position.x/m_BoxLength + m_Size.x/2);
		int y_index = int(position.y/m_BoxLength + m_Size.y/2);
		int z_index = int(position.z/m_BoxLength);
		if (!((x_index<m_Size.x)&(y_index<m_Size.y)&(z_index<m_Size.z)&(x_index>=0)&(y_index>=0)&(z_index>=0)))
		  printf("%d %d %d %d %d %d %g %g %g %g\n",x_index, y_index, z_index, m_Size.x, m_Size.y, m_Size.z, position.x, position.y, position.z, m_BoxLength);
		  
		MyAssert((x_index<m_Size.x)&(y_index<m_Size.y)&(z_index<m_Size.z)&(x_index>=0)&(y_index>=0)&(z_index>=0),"Cells have expanded farther than boxes allows");
		
		return Address(x_index*(m_Size.y*m_Size.z) + y_index*m_Size.z + z_index); // check this!!
	}
	
	IntCoord Size()
	{
		return m_Size;
	}

	IntCoord GetXY(Address address)
	{
	    int a = address.a;
	    int x_index = a/(m_Size.y*m_Size.z);
	    a -= x_index*(m_Size.y*m_Size.z);
	    int y_index = a/m_Size.z;
	    int z_index = a%m_Size.z;
		return IntCoord(x_index, y_index, z_index);
	}
	
	DoubleCoord GetCentroid(IntCoord index)
	{

		return DoubleCoord( (float(index.x) - m_Size.x/2.0 + 0.5)*m_BoxLength, (float(index.y) - m_Size.y/2.0 + 0.5)*m_BoxLength, (float(index.z) + 0.5)*m_BoxLength );
	}

	void Add(int cellID, const Address& address)	// add a cell to the correct box
	{
 		m_Boxes[address.a]->Add(cellID);
	}
	
	void Remove(int cellID, const Address& address)	// remove cell from a box
	{
		m_Boxes[address.a]->Remove(cellID);
	}
	
	void Move(int cellID, const Address& oldAddress, const Address& newAddress)
	{
		if (oldAddress.a != newAddress.a)
		{
			m_Boxes[oldAddress.a]->Remove(cellID);
			m_Boxes[newAddress.a]->Add(cellID);
		}
	}
	
	int GetNeighbours(int cellID, const Address& address, int* neighbours, int maxNeighbours)
	{
		int neighboursFound = 1;
		Box* box;
		int neighbourID;
		
		IntCoord a0 = GetXY(address);
		// loop through neighbouring indices to find neighbours
		for (int za = max(a0.z-1,0); za<min(a0.z+2,m_Size.z); za++)
		{
			for (int ya = max(a0.y-1,0); ya<min(a0.y+2,m_Size.y); ya++)
			{
				for (int xa = max(a0.x-1,0); xa<min(a0.x+2,m_Size.x); xa++)
				{
					box = m_Boxes[Index(xa,ya,za)];	// pointer to box as this address
					for (int cella = 0; cella<box->Count; cella++)
					{
						neighbourID = box->Cells[cella];
						if (neighbourID!=cellID)
						{
							neighbours[neighboursFound++] = box->Cells[cella];
							MyAssert(maxNeighbours>neighboursFound,"Too many neighbours!");
						}
					}
				}
			}
		}
		neighbours[0] = neighboursFound-1; // zeroth element is always number of neighbours
		return neighboursFound-1;
	}
	
	int GetCells(int x_index, int y_index, int z_index, int*& cells)
	{
		int a = Index(x_index,y_index,z_index);
		int count = m_Boxes[a]->Count;
		cells = m_Boxes[a]->Cells;
		return count;
	}
	
	int GetNumber(int x_index, int y_index, int z_index)
	{
		int count = m_Boxes[Index(x_index,y_index,z_index)]->Count;
		return count;
	}
	
	~UniformGrid(void)
	{
		for(int z = 0; z<m_Size.z; z++)
		{
			for(int x=0; x<m_Size.x; x++)
			{
				for(int y=0; y<m_Size.y; y++)
				{
					delete m_Boxes[x*(m_Size.y*m_Size.z) + y*m_Size.z + z];	// delete individual boxes
				}
			}
		}

		delete[] m_Boxes;	// delete array of box pointers
	}
	
private:
	int Index(const int x, const int y, const int z) {	return x*(m_Size.y*m_Size.z) + y*m_Size.z + z; }
	IntCoord m_Size;
	int m_MaxCellsPerBox;
	double m_BoxLength;
	Box** m_Boxes;
};


#endif // _UNIFORM_GRID_H_
