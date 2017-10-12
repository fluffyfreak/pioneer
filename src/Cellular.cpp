// Copyright © 2008-2016 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Cellular.h"
#include "Utils.h"
#include "MathUtil.h"
#include <algorithm>



Cellular::Cellular(const int cellSize, const int dimx, const int dimy, const std::vector<vector2d> &s, const std::vector<double> &h) 
	: CELL_DIVISOR(cellSize), INV_CELL_DIVISOR(1.0/double(cellSize))
	, size_x(dimx), size_y(dimy), half_size_x(dimx>>1), half_size_y(dimy>>1)
	, cellsX(dimx / cellSize), cellsY(dimy / cellSize)
	, heights(h) 
{
	PROFILE_SCOPED()
	cells.reset(new Cell[cellsX * cellsY]);
	const size_t count = s.size();
	for( size_t c=0; c<count; c++ )
	{
		// find which Cell to put the values into
		const Uint32 xi = s[c].x * INV_CELL_DIVISOR;
		const Uint32 yi = s[c].y * INV_CELL_DIVISOR;
		const Uint32 index = (yi * cellsX) + xi;
		cells[index].AddEntry(s[c],c);
	}
#ifdef _DEBUG
	Uint32 maxPoint=0;
	Uint32 minPoint=1000;
	Uint32 ptCount[CELL_ENTRIES] = {0,0,0,0};
	for( size_t c=0; c<(cellsX * cellsY); c++ ) {
		maxPoint = std::max(maxPoint, cells[c].count);
		minPoint = std::min(minPoint, cells[c].count);
		ptCount[cells[c].count]++;
	}
	Output("Cellular :: Max Count in Cell %u, and Min %u\n", maxPoint, minPoint);
	for( int i=0; i<CELL_ENTRIES; i++ ) {
		Output("Cellular :: Point count %d, instances %u\n", i, ptCount[i]);
	}
#endif // _DEBUG
	Output("Cellular :: CellsX %u, CellsY %u\n", cellsX, cellsY);

	GenMap();
}

#define CELL_OFFSET 2
size_t Cellular::NearestSite( const int x, const int y ) const
{
	PROFILE_SCOPED()
	size_t site=0xFFFFFFFF;
 	double mindist = DBL_MAX;

	// find this cell
	const int cellX = x * INV_CELL_DIVISOR;
	const int cellY = y * INV_CELL_DIVISOR;
	const int cellXStart = (cellX-CELL_OFFSET); 
	const int cellXEnd = (cellX+CELL_OFFSET);
	const int cellYStart = std::max(cellY-CELL_OFFSET,0); 
	const int cellYEnd = std::min(cellY+CELL_OFFSET,int(cellsY));
	const double xd = x;
	const double yd = y;
	const vector2d xyd(x,y);
	if(cellXStart<0 || cellXEnd>=cellsX) 
	{
		PROFILE_SCOPED_DESC("NearestSiteWrap")
		// search through the NxN grid of cells centred on cellX|cellY;
		for(int yi = cellYStart; yi<cellYEnd; yi++) 
		{
			for(int xi = cellXStart; xi<cellXEnd; xi++) 
			{
				const Uint32 index = (yi * cellsX) + MathUtil::iwrap(xi, cellsX);
				const Cell &cell = cells[index];
				const size_t count = cell.count;
				for (size_t i=0 ; i<count ; ++i)
				{
					const double dist = WrapHorizontalDist(xyd,cell.points[i]);
					if (dist < mindist) 
					{
						mindist = dist;
						site = cell.indices[i];
					}
				}
			}
		}
	}
	else 
	{
		PROFILE_SCOPED_DESC("NearestSiteDist")
		// search through the NxN grid of cells centred on cellX|cellY;
		for(int yi = cellYStart; yi<cellYEnd; yi++) 
		{
			for(int xi = cellXStart; xi<cellXEnd; xi++) 
			{
				const Uint32 index = (yi * cellsX) + MathUtil::iwrap(xi, cellsX);
				const Cell &cell = cells[index];
				const size_t count = cell.count;
				for (size_t i=0 ; i<count ; ++i)
				{
					const double dist = Dist(xyd,cell.points[i]);
					if (dist < mindist) 
					{
						mindist = dist;
						site = cell.indices[i];
					}
				}
			}
		}
	}
 	return site;
}
	
void Cellular::GenMap()
{
	PROFILE_SCOPED()
	buf.reset(new double[size_y * size_x]);
	double *ptr = buf.get();

	for (int i = 0; i < size_y; i++ ) {
		for (int j = 0; j < size_x; j++ ) {
			(*ptr) = heights[NearestSite(j, i)];
			++ptr;
		}
	}
}
