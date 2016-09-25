// Copyright © 2008-2016 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#pragma once

#ifndef __CELLULAR_H__
#define __CELLULAR_H__

#include "libs.h"


class Cellular {
public:
	Cellular(const int cellSize, const int dimx, const int dimy, const std::vector<vector2d> &s, const std::vector<double> &h);

	double* CellularMap() const { return buf.get(); }
private:
	const int CELL_DIVISOR;
	const double INV_CELL_DIVISOR;

	const int size_x, size_y, half_size_x, half_size_y;
	const int cellsX, cellsY;
	const std::vector<double> &heights;
	std::unique_ptr<double[]> buf;

	class Cell 
	{
	public:
		// Hardcoded for the size of cell and density of points
		#define CELL_ENTRIES 4
		Cell() : count(0) {}
		__forceinline void AddEntry(const vector2d &point, const size_t index) {
			assert(count<CELL_ENTRIES);
			points[count] = point;
			indices[count] = index;
			++count;
		}
		vector2d	points[CELL_ENTRIES];
		size_t		indices[CELL_ENTRIES];
		Uint32		count;
	};
	std::unique_ptr<Cell[]> cells;
	
	__forceinline double WrapHorizontalDist( const vector2d &xy, const vector2d &p) const
	{
		double dx = fabs(xy.x-p.x);
		const double dy = (xy.y-p.y);
		// only wrap on the horizontal
		dx = (dx > half_size_x) ? (size_x-dx) : dx;
		return (dx*dx) + (dy*dy);
	}
	
	__forceinline double Dist( const vector2d &xy, const vector2d &p) const
	{
		return (xy-p).LengthSqr();
	}

	size_t NearestSite( const int x, const int y ) const;
	
	void GenMap();
};

#endif // __CELLULAR_H__
