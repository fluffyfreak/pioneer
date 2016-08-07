// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Terrain.h"
#include "TerrainNode.h"

#include "SDLWrappers.h"
#include "Easing.h"

static SDLSurfacePtr im;


template <>
const char *TerrainColorFractal<TerrainColourJSON>::GetColorFractalName() const { return "JSON"; }

template <>
TerrainColorFractal<TerrainColourJSON>::TerrainColorFractal(const SystemBody *body, const std::string &JSONColourFile) : Terrain(body)
{
	LoadTerrainJSON(FileSystem::JoinPathBelow("terrain", JSONColourFile), m_terrainSrcs);
	if(!im.Valid()) {
		im = LoadSurfaceFromFile(FileSystem::JoinPathBelow("terrain", "TerrainColor.png"));
	}
	
	m_surfaceEffects |= Terrain::EFFECT_WATER;
}

Color4ub GetPixel(const SDLSurfacePtr& surface, const int x, const int y)
{
    const int bpp = surface->format->BytesPerPixel;
    // Here p is the address to the pixel we want to retrieve
    const Uint8 *p = (Uint8 *)surface->pixels + y * surface->pitch + x * bpp;
    switch(bpp) {
    case 1:
        return Color4ub(*p,*p,*p);
    case 3:
#if (SDL_BYTEORDER == SDL_BIG_ENDIAN)
        return p[0] << 16 | p[1] << 8 | p[2];
#else
        return Color4ub(p[0], p[1], p[2]);
#endif
    case 4:
        return Color4ub(*(Uint32 *)p);
    }
	return Color4ub::BLANK;
}

static double radiusScale = 0.001;
template <>
vector3d TerrainColorFractal<TerrainColourJSON>::GetColor(const vector3d &p, double height, const vector3d &norm) const
{
	if (height <= 0.000001) {
		return vector3d(0.0, 0.0, 1.0);
	}

	// get lat and long
	double latitude = -asin(p.y);
	if (p.y < -1.0) latitude = -0.5*M_PI;
	if (p.y > 1.0) latitude = 0.5*M_PI;
	double longitude = atan2(p.x, p.z);

	// TODO
	// figure out temperature modification from latitude
	// figure out temperature modification from height(altitude)
	// modify temperature +/- based on above
	const double low=0.0, high=1.0;
	const double tempMul = Clamp(Easing::Cubic::EaseOut(1.0-abs(latitude), low, high, high-low), 0.0, 1.0);

	const vector3d posRadius(p*(m_planetRadius * radiusScale));

	double humidity=0.0;
	double temperature=0.0;
	double jitter=0.0;
	for (auto ts : m_terrainSrcs)
	{
		switch ( ts.Type() )
		{
		case TerrainSource::ST_HUMIDITY:
			// initialise the height
			humidity = ts.BaseHeight();

			// for each height source
			for (auto nodex : ts.Nodes())
			{
				// for each node
				humidity += nodex.Call(posRadius,p);
			}
			break;
		case TerrainSource::ST_TEMPERATURE:
			// initialise the height
			temperature = ts.BaseHeight();

			// for each height source
			for (auto nodex : ts.Nodes())
			{
				// for each node
				temperature += nodex.Call(posRadius,p);
			}
			break;
		case TerrainSource::ST_JITTER:
			// initialise the height
			jitter = ts.BaseHeight();

			// for each height source
			for (auto nodex : ts.Nodes())
			{
				// for each node
				jitter += nodex.Call(posRadius,p);
			}
			break;
		default:
			continue;
		}
	}

	const Color4f col = GetPixel(im, ::Clamp(temperature*(tempMul+jitter),0.0,255.0), humidity).ToColor4f();
	return vector3d(col.r, col.g, col.b);
}
