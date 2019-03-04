// Copyright © 2008-2019 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#pragma once

#ifndef _BILINEAR_INTERPOLATION_H
#define _BILINEAR_INTERPOLATION_H

#include "MathUtil.h"
template< class T >
inline T Bilinear( 
	const float &tx, 
	const float &ty, 
	const T &c00, 
	const T &c10, 
	const T &c01, 
	const T &c11) 
{
	T a = MathUtil::mix(c00, c10, tx);
	T b = MathUtil::mix(c01, c11, tx);
	return MathUtil::mix(a, b, ty);
}

template< class T >
T BilinearSample(const int imgDim, const int imgX, const int imgY, const int gridSize, const T *grid2d)
{
	// convert i,j to grid coordinates
	const float gx = imgX / float(imgDim) * (gridSize - 1); // be careful to interpolate boundaries 
	const float gy = imgY / float(imgDim) * (gridSize - 1); // be careful to interpolate boundaries 
	const int gxi = int(gx); 
	const int gyi = int(gy); 
	const T &c00 = grid2d[gyi * gridSize + gxi]; 
	const T &c10 = grid2d[gyi * gridSize + (gxi + 1)]; 
	const T &c01 = grid2d[(gyi + 1) * gridSize + gxi]; 
	const T &c11 = grid2d[(gyi + 1) * gridSize + (gxi + 1)]; 
	return Bilinear(gx - gxi, gy - gyi, c00, c10, c01, c11); 
}

#endif // _BILINEAR_INTERPOLATION_H
