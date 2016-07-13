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
	const Uint32 cellsX, cellsY;
	const std::vector<double> &heights;
	std::unique_ptr<double[]> buf;

	class Cell 
	{
	public:
		Cell() {
			PROFILE_SCOPED()
			points.reserve(20);
			indices.reserve(20);
		}
		std::vector<vector2d>	points;
		std::vector<size_t>		indices;
	};
	std::unique_ptr<Cell[]> cells;
	
	inline double WrapDist( int x, int y, const vector2d &p) const
	{
		PROFILE_SCOPED()
		double dx = abs(x-p.x);
		double dy = abs(y-p.y);
		if (dx > half_size_x ) // only wrap on the horizontal
			dx = size_x-dx;
		//if (dy > half_size_y )
		//	dy = size_y-dy;
		// return squared distance
		return dx*dx + dy*dy;
	}

#define CELL_OFFSET 2
	size_t NearestSite( const int x, const int y ) const;
	
	void GenMap();
};

#endif // __CELLULAR_H__
