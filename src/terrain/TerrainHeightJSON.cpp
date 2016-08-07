// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Terrain.h"
#include "TerrainNode.h"

template <>
const char *TerrainHeightFractal<TerrainHeightJSON>::GetHeightFractalName() const { return "JSON"; }

template <>
TerrainHeightFractal<TerrainHeightJSON>::TerrainHeightFractal(const SystemBody *body, const std::string &JSONHeightFile) : Terrain(body)
{
	LoadTerrainJSON(FileSystem::JoinPathBelow("terrain", JSONHeightFile), m_terrainSrcs);
}

static double heightScale = 1.0 / 15000000.0;
static double radiusScale = 0.001;
template <>
double TerrainHeightFractal<TerrainHeightJSON>::GetHeight(const vector3d &p) const
{
	double n = 0.0;

	//const vector3d posRadius(p * radiusScale);// *(m_planetRadius * radiusScale));
	const vector3d posRadius(p*(m_planetRadius * radiusScale));

	for (auto ts : m_terrainSrcs)
	{
		if (ts.Type() == TerrainSource::ST_HEIGHT)
		{
			// initialise the height
			n = ts.BaseHeight();

			// for each height source
			for (auto nodex : ts.Nodes())
			{
				// for each node
				n += nodex.Call(posRadius, p);
			}
			break;
		}
	}

	return n * heightScale;
}
