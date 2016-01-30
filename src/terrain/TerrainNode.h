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
		NT_NOISE_CUBED,
		NT_CONSTANT,
		NT_SQUARED,
		NT_CUBED
	};

	// public methods
	TerrainNodeData() : 
		m_op(TO_ADD), 
		m_scaleHigh(0.0), 
		m_scaleLow(0.0),
		m_bScale(false),
		m_octaves(0),
		m_frequency(0.0),
		m_persistence(0.0),
		m_clamp(std::make_pair(DBL_MIN, DBL_MAX)),
		m_bClamp(false),
		m_noiseType(NT_NOISE)
	{
	}

	void Name(const std::string& str) {
		m_name = str;
	}

	void Frequency(const double freq) {
		m_frequency = freq;
	}

	void ScaleLow(const double low) {
		m_scaleLow = low;
	}

	void ScaleHigh(const double high) {
		m_scaleHigh = high;
		m_bScale = true;
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
		if (str == "noise_cellular_squared") {
			m_noiseType = NT_NOISE_CELLULAR_SQUARED;
		} else if (str == "noise_ridged") {
			m_noiseType = NT_NOISE_RIDGED;
		} else if (str == "noise_cubed") {
			m_noiseType = NT_NOISE_CUBED;
		} else if (str == "constant") {
			m_noiseType = NT_CONSTANT;
		} else if (str == "squared") {
			m_noiseType = NT_SQUARED;
		} else if (str == "cubed") {
			m_noiseType = NT_CUBED;
		}
	}

	void ClampNoise(const double lower, const double upper) {
		m_clamp.first = lower;
		m_clamp.second = upper;
		m_bClamp = true;
	}

	void AddChild(const TerrainNodeData& child) {
		m_children.push_back(child);
	}

	double Call(const vector3d& p);

private:
	inline double Scale(const double h)
	{
		// convert to 0..1 from -1..1
		return (m_bScale) ? MathUtil::mix(m_scaleLow, m_scaleHigh, ((1.0 + h)*0.5)) : h;
	}

	//"name": "Mountains",
	std::string m_name;

	//"op": "mul",
	EnumTerrainOp m_op;

	//"high": 13,
	//"low": -13,
	double m_scaleHigh;
	double m_scaleLow;
	bool m_bScale;

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
	bool m_bClamp;

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
