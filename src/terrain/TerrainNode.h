// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Terrain.h"
#include "TerrainNoise.h"
#include "TerrainFeature.h"

#include <math.h>
#include "MathUtil.h"

#include "FileSystem.h"
#include "json/json.h"

#include "vector2.h"

class TerrainNodeData
{
public:
	// public data types
	enum EnumTerrainOp {
		TO_ADD = 0,
		TO_SUB,
		TO_MUL,
		TO_DIV
	};

	enum EnumNoiseTypes {
		NT_NOISE = 0,
		NT_NOISE_CELLULAR_SQUARED,
		NT_NOISE_RIDGED,
		NT_NOISE_CUBED
	};

	// public methods
	TerrainNodeData() : 
		m_op(TO_ADD), 
		m_scale(std::make_pair(0, 0)),
		m_octaves(0),
		m_frequency(0.0),
		m_persistence(0.0),
		m_clamp(std::make_pair(DBL_MIN, DBL_MAX)),
		m_noiseType(NT_NOISE)
	{
	}

	void Name(const std::string& str) {
		m_name = str;
	}

	void Frequency(const double freq) {
		m_frequency = freq;
	}

	void Scale(const double lower, const double upper) {
		m_scale.first = lower;
		m_scale.second = upper;
	}

	void Octaves(const int oct) {
		m_octaves = oct;
	}

	void Op(const EnumTerrainOp op) {
		m_op = op;
	}

	void Persistence(const double per) {
		m_persistence = per;
	}

	void NoiseType(const std::string& str) {
		if (str == "noise") {
			m_noiseType = NT_NOISE;
		} else if (str == "noise_cellular_squared") {
			m_noiseType = NT_NOISE_CELLULAR_SQUARED;
		} else if (str == "noise_ridged") {
			m_noiseType = NT_NOISE_RIDGED;
		} else if (str == "noise_cubed") {
			m_noiseType = NT_NOISE_CUBED;
		} 
	}

	void ClampNoise(const double lower, const double upper) {
		m_clamp.first = lower;
		m_clamp.second = upper;
	}

	void AddChild(const TerrainNodeData& child) {
		m_children.push_back(child);
	}

	double Call(const vector3d& p);

private:
	inline double Scale(const double h)
	{
		// convert to 0..1 from -1..1
		return MathUtil::mix(m_scale.first, m_scale.second, ((1.0 + h)*0.5));
	}

	inline double Clamp(const double h)
	{
		return ::Clamp(h, m_clamp.first, m_clamp.second);
	}

	//"name": "Mountains",
	std::string m_name;

	//"op": "mul",
	EnumTerrainOp m_op;

	//"scale": [-13,13],
	std::pair<double, double> m_scale;

	//"octaves": 7,
	int m_octaves;

	//"frequency": 0.002,
	double m_frequency;

	//"persistence": 0.7,
	double m_persistence;

	//"type": "noise_cubed"
	EnumNoiseTypes m_noiseType;

	//"clamp": [ 0, 1 ]
	std::pair<double, double> m_clamp;

	//"children": [ ... ]
	std::vector<TerrainNodeData> m_children;
};

class TerrainSource
{
public:
	enum EnumSourceTypes {
		ST_HEIGHT=0,
		ST_HUMIDITY,
		ST_TEMPERATURE
	};

	void SetType(const EnumSourceTypes type)
	{
		m_type = type;
	}

	void SetBaseHeight(const double base)
	{
		m_baseHeight = base;
	}

	void AddNode(const TerrainNodeData& node)
	{
		m_terrainNodes.push_back(node);
	}

	EnumSourceTypes Type() const { return m_type; }
	double BaseHeight() const { return m_baseHeight; }
	const std::vector<TerrainNodeData>& Nodes() { return m_terrainNodes; }
private:
	EnumSourceTypes m_type;
	double m_baseHeight;
	std::vector<TerrainNodeData> m_terrainNodes;
};

void LoadTerrainJSON(const std::string& path, std::vector<TerrainSource>& sources);
