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
		cells[index].points.push_back(s[c]);
		cells[index].indices.push_back(c);
	}
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
	// search through the NxN grid of cells centred on cellX|cellY;
	for(int yi = std::max(cellY-CELL_OFFSET,0); yi<std::min(cellY+CELL_OFFSET,int(cellsY)); yi++) 
	{
		for(int xi = (cellX-CELL_OFFSET); xi<(cellX+CELL_OFFSET); xi++) 
		{
			const Uint32 index = (yi * cellsX) + MathUtil::iwrap(xi, cellsX);
			const std::vector<vector2d> &pts = cells[index].points;
			for (size_t i=0 ; i<pts.size() ; ++i)
			{
				double dist = WrapDist(x,y,pts[i]);
				if (dist < mindist) 
				{
					mindist = dist;
					site = cells[index].indices[i];
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
