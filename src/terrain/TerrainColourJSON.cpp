// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Terrain.h"
#include "TerrainNode.h"

static std::vector<TerrainSource> terrainSrcs;

template <>
const char *TerrainColorFractal<TerrainColourJSON>::GetColorFractalName() const { return "JSON"; }

template <>
TerrainColorFractal<TerrainColourJSON>::TerrainColorFractal(const SystemBody *body, const std::string &JSONColourFile) : Terrain(body)
{
	LoadTerrainJSON(FileSystem::JoinPathBelow("terrain", JSONColourFile), m_terrainSrcs);
}

template <>
vector3d TerrainColorFractal<TerrainColourJSON>::GetColor(const vector3d &p, double height, const vector3d &norm) const
{
	double n = 0.0;
	//for (auto ts : terrainSrcs)
	//{
	//	if (ts.Type() == TerrainSource::ST_HEIGHT)
	//	{
	//		// initialise the height
	//		n = ts.BaseHeight();
	//	}
	//}
	if (height > n) {
		return vector3d(1.0, 1.0, 1.0);
	}
	return vector3d(0.0, 0.0, 1.0);

	/*double n = 0.0;

	const vector3d posRadius(p * posScale);

	for (auto ts : terrainSrcs)
	{
	ST_HUMIDITY,
	ST_TEMPERATURE
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

	return n / 3000000.0;//3000000.0;
	*/
}
