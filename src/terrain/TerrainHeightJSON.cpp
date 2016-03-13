// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Terrain.h"
#include "TerrainNode.h"

static std::vector<TerrainSource> terrainSrcs;

template <>
const char *TerrainHeightFractal<TerrainHeightJSON>::GetHeightFractalName() const { return "JSON"; }

template <>
TerrainHeightFractal<TerrainHeightJSON>::TerrainHeightFractal(const SystemBody *body, const std::string &JSONfilename, const std::string &heightmapFilename) : Terrain(body)
{
	const std::string path("terrain/Terra.json");
	LoadTerrainJSON(path, terrainSrcs);
}

static double heightScale = 1.0 / 15000000.0;
static double radiusScale = 0.001;
template <>
double TerrainHeightFractal<TerrainHeightJSON>::GetHeight(const vector3d &p) const
{
	double n = 0.0;

	//const vector3d posRadius(p * radiusScale);// *(m_planetRadius * radiusScale));
	const vector3d posRadius(p*(m_planetRadius * radiusScale));

	for (auto ts : terrainSrcs)
	{
		if (ts.Type() == TerrainSource::ST_HEIGHT)
		{
			// initialise the height
			n = ts.BaseHeight();

			// for each height source
			for (auto nodex : ts.Nodes())
			{
				// for each node
				n += nodex.Call(posRadius);
			}
			break;
		}
	}

	return n * heightScale;
}
