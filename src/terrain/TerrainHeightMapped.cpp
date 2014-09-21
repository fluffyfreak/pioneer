// Copyright © 2008-2014 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Terrain.h"
#include "TerrainNoise.h"

using namespace TerrainNoise;

template <>
const char *TerrainHeightFractal<TerrainHeightMapped>::GetHeightFractalName() const { return "Mapped"; }

template <>
TerrainHeightFractal<TerrainHeightMapped>::TerrainHeightFractal(const SystemBody *body) : Terrain(body)
{
	//textures
	if (textures) {
		SetFracDef(0, m_maxHeightInMeters, 10, 10*m_fracmult);
		SetFracDef(1, m_maxHeightInMeters, 25, 10*m_fracmult);
	}
	//small fractal/high detail
	SetFracDef(2-m_fracnum, m_maxHeightInMeters*0.0000005, 50, 20*m_fracmult);//[2]
	//continental/large type fractal
	SetFracDef(3-m_fracnum, m_maxHeightInMeters*0.00005, 1e6, 800*m_fracmult);//[0]
	SetFracDef(4-m_fracnum, m_maxHeightInMeters*0.00005, 1e5, 400*m_fracmult);//[4]
	//medium fractal
	SetFracDef(5-m_fracnum, m_maxHeightInMeters*0.000005, 2e4, 200*m_fracmult);//[5]
	SetFracDef(6-m_fracnum, m_maxHeightInMeters*0.0000005, 5e3, 100*m_fracmult);//[3]
}

template <>
double TerrainHeightFractal<TerrainHeightMapped>::GetHeight(const vector3d &p) const
{
    // This is all used for Earth and Earth alone

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

	v = (v<0 ? 0 : v);
	double h = v;

	//Here's where we add some noise over the heightmap so it doesnt look so boring, we scale by height so values are greater high up
	//large mountainous shapes
	double mountains = h*h*0.001*octavenoise(GetFracDef(3-m_fracnum), 0.5*octavenoise(GetFracDef(5-m_fracnum), 0.45, p),
		p)*ridged_octavenoise(GetFracDef(4-m_fracnum), 0.475*octavenoise(GetFracDef(6-m_fracnum), 0.4, p), p);
	v += mountains;
	//smaller ridged mountains
	if (v < 50.0){
		v += v*v*0.04*ridged_octavenoise(GetFracDef(5-m_fracnum), 0.5, p);
	} else if (v <100.0){
		v += 100.0*ridged_octavenoise(GetFracDef(5-m_fracnum), 0.5, p);
	} else {
		v += (100.0/v)*(100.0/v)*(100.0/v)*(100.0/v)*(100.0/v)*
			100.0*ridged_octavenoise(GetFracDef(5-m_fracnum), 0.5, p);
	}

	//low altitude detail/dunes
	if (v < 10.0){
		v += 2.0*v*dunes_octavenoise(GetFracDef(6-m_fracnum), 0.5, p)
			*octavenoise(GetFracDef(6-m_fracnum), 0.5, p);
	} else if (v <50.0){
		v += 20.0*dunes_octavenoise(GetFracDef(6-m_fracnum), 0.5, p)
			*octavenoise(GetFracDef(6-m_fracnum), 0.5, p);
	} else {
		v += (50.0/v)*(50.0/v)*(50.0/v)*(50.0/v)*(50.0/v)
			*20.0*dunes_octavenoise(GetFracDef(6-m_fracnum), 0.5, p)
			*octavenoise(GetFracDef(6-m_fracnum), 0.5, p);
	}
	if (v <60.0){
		v += (v-40.0)*billow_octavenoise(GetFracDef(5-m_fracnum), 0.5, p);
	} else {
		v += (30.0/v)*(30.0/v)*(30.0/v)*20.0*billow_octavenoise(GetFracDef(5-m_fracnum), 0.5, p);
	}

	//ridges and bumps
	v += h*0.2*voronoiscam_octavenoise(GetFracDef(5-m_fracnum), Clamp(1.0-(h*0.0002), 0.0, 0.6), p)
		* Clamp(1.0-(h*0.0006), 0.0, 1.0);
	//polar ice caps with cracks
	if ((m_icyness*0.5)+(fabs(p.y*p.y*p.y*0.38)) > 0.6) {
		h = Clamp(1.0-(v*10.0), 0.0, 1.0)*voronoiscam_octavenoise(GetFracDef(5-m_fracnum), 0.5, p);
		h *= h*h*2.0;
		h -= 3.0;
		v += h;
	}

	return v<0 ? 0 : (v*m_invPlanetRadius);
}
