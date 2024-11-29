// Copyright © 2008-2024 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Terrain.h"
#include "TerrainNoise.h"

using namespace TerrainNoise;

template <>
const char *TerrainColorFractal<TerrainColorIce>::GetColorFractalName() const { return "Ice"; }

template <>
TerrainColorFractal<TerrainColorIce>::TerrainColorFractal(const SystemBody *body) :
	Terrain(body)
{
}
