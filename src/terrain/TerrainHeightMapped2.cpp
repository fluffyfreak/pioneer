// Copyright © 2008-2014 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Terrain.h"
#include "TerrainNoise.h"

using namespace TerrainNoise;
template <>
const char *TerrainHeightFractal<TerrainHeightMapped2>::GetHeightFractalName() const { return "Mapped2"; }

template <>
TerrainHeightFractal<TerrainHeightMapped2>::TerrainHeightFractal(const SystemBody *body) : Terrain(body)
{
}

template <>
double TerrainHeightFractal<TerrainHeightMapped2>::GetHeight(const vector3d &p) const
{

	double latitude = -asin(p.y);
	if (p.y < -1.0) latitude = -0.5*M_PI;
	if (p.y > 1.0) latitude = 0.5*M_PI;

	const double longitude = atan2(p.x, p.z);
	const double px = (((m_heightMapSizeX-1) * (longitude + M_PI)) / (2*M_PI));
	const double py = ((m_heightMapSizeY-1)*(latitude + 0.5*M_PI)) / M_PI;
	const int ix = Clamp(int(floor(px)), 0, m_heightMapSizeX-1);
	const int iy = Clamp(int(floor(py)), 0, m_heightMapSizeY-1);
	const double dx = px-ix;
	const double dy = py-iy;

	// p0,3 p1,3 p2,3 p3,3
	// p0,2 p1,2 p2,2 p3,2
	// p0,1 p1,1 p2,1 p3,1
	// p0,0 p1,0 p2,0 p3,0
	double map[4][4];
	const float *pHMap = m_heightMap.get();
	for (int x=-1; x<3; x++) {
		for (int y=-1; y<3; y++) {
			const int xx = Wrap(ix+x, 0, (m_heightMapSizeX-1));
			const int yy = Wrap(iy+y, 0, (m_heightMapSizeY-1));
			map[x+1][y+1] = pHMap[(yy * m_heightMapSizeX) + xx];
		}
	}

	double v = DoThisThing(dx, dy, map);

	v=v*m_heightScaling+m_minh; // v = v*height scaling+min height
	v*=m_invPlanetRadius;

	v += 0.1;
	double h = 1.5*v*v*v*ridged_octavenoise(16, 4.0*v, 4.0, p);
	h += 30000.0*v*v*v*v*v*v*v*ridged_octavenoise(16, 5.0*v, 20.0*v, p);
	h += v;
	h -= 0.09;

	return (h > 0.0 ? h : 0.0);


}

