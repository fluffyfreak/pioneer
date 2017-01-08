// Copyright © 2008-2017 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _TERRAIN_H
#define _TERRAIN_H

#include "libs.h"
#include "galaxy/StarSystem.h"

#include "TerrainNode.h"

#ifdef _MSC_VER
#pragma warning(disable : 4250)			// workaround for MSVC 2008 multiple inheritance bug
#endif

struct fracdef_t {
	fracdef_t() : amplitude(0.0), frequency(0.0), lacunarity(0.0), octaves(0) {}
	double amplitude;
	double frequency;
	double lacunarity;
	int octaves;
};


template <typename,typename> class TerrainGenerator;


class Terrain : public RefCounted {
public:
	// location and intensity of effects are controlled by the colour fractals;
	// it's possible for a Terrain to have a flag set but not actually to exhibit any of that effect
	enum SurfaceEffectFlags {
		EFFECT_LAVA  = 1 << 0,
		EFFECT_WATER = 2
		// can add other effect flags here (e.g., water, snow, ice)
	};

	static Terrain *InstanceTerrain(const SystemBody *body);

	virtual ~Terrain();

	void SetFracDef(const unsigned int index, const double featureHeightMeters, const double featureWidthMeters, const double smallestOctaveMeters = 20.0);
	inline const fracdef_t &GetFracDef(const unsigned int index) const { assert(index>=0 && index<MAX_FRACDEFS); return m_fracdef[index]; }

	virtual double GetHeight(const vector3d &p) const = 0;
	virtual vector3d GetColor(const vector3d &p, double height, const vector3d &norm) const = 0;

	virtual const char *GetHeightFractalName() const = 0;
	virtual const char *GetColorFractalName() const = 0;

	double GetMaxHeight() const { return m_maxHeight; }

	Uint32 GetSurfaceEffects() const { return m_surfaceEffects; }

	void DebugDump() const;

private:
	template <typename HeightFractal, typename ColorFractal>
	static Terrain *InstanceGenerator(const SystemBody *body, const std::string &JSONHeightFile, const std::string &JSONColourFile) { return new TerrainGenerator<HeightFractal,ColorFractal>(body, JSONHeightFile, JSONColourFile); }

	typedef Terrain* (*GeneratorInstancer)(const SystemBody *, const std::string &JSONHeightFile, const std::string &JSONColourFile);


protected:
	Terrain(const SystemBody *body);

	bool textures;
	int m_fracnum;
	double m_fracmult;

	Uint32 m_seed;
	Random m_rand;

	double m_sealevel; // 0 - no water, 1 - 100% coverage
	double m_icyness; // 0 - 1 (0% to 100% cover)
	double m_volcanic;

	Uint32 m_surfaceEffects;

	/** General attributes */
	double m_maxHeight;
	double m_maxHeightInMeters;
	double m_invMaxHeight;
	double m_planetRadius;
	double m_planetEarthRadii;

	double m_entropy[12];

	vector3d m_rockColor[8];
	vector3d m_darkrockColor[8];
	vector3d m_greyrockColor[8];
	vector3d m_plantColor[8];
	vector3d m_darkplantColor[8];
	vector3d m_sandColor[8];
	vector3d m_darksandColor[8];
	vector3d m_dirtColor[8];
	vector3d m_darkdirtColor[8];
	vector3d m_gglightColor[8];
	vector3d m_ggdarkColor[8];

	/* XXX you probably shouldn't increase this. If you are
	   using more than 10 then things will be slow as hell */
	static const Uint32 MAX_FRACDEFS = 10;
	fracdef_t m_fracdef[MAX_FRACDEFS];

	struct MinBodyData {
		MinBodyData(const SystemBody* body) {
			m_radius = body->GetRadius();
			m_aspectRatio = body->GetAspectRatio();
			m_path = body->GetPath();
			m_name = body->GetName();
		}
		double m_radius;
		double m_aspectRatio;
		SystemPath m_path;
		std::string m_name;
	};
	MinBodyData m_minBody;
	
	std::vector<TerrainSource> m_terrainSrcs;
};


template <typename HeightFractal>
class TerrainHeightFractal : virtual public Terrain {
public:
	virtual double GetHeight(const vector3d &p) const;
	virtual const char *GetHeightFractalName() const;
protected:
	TerrainHeightFractal(const SystemBody *body, const std::string &JSONHeightFile);
private:
	TerrainHeightFractal() {}
};

template <typename ColorFractal>
class TerrainColorFractal : virtual public Terrain {
public:
	virtual vector3d GetColor(const vector3d &p, double height, const vector3d &norm) const;
	virtual const char *GetColorFractalName() const;
protected:
	TerrainColorFractal(const SystemBody *body, const std::string &JSONColourFile);
private:
	TerrainColorFractal() {}
};


template <typename HeightFractal, typename ColorFractal>
class TerrainGenerator : public TerrainHeightFractal<HeightFractal>, public TerrainColorFractal<ColorFractal> {
public:
	TerrainGenerator(const SystemBody *body, const std::string &JSONHeightFile, const std::string &JSONColourFile) : Terrain(body), TerrainHeightFractal<HeightFractal>(body, JSONHeightFile), TerrainColorFractal<ColorFractal>(body, JSONColourFile) {}

private:
	TerrainGenerator() {}
};

// CPU side JSON height based generation:
class TerrainHeightJSON;
class TerrainColourJSON;

#ifdef _MSC_VER
#pragma warning(default : 4250)
#endif

#endif /* TERRAIN_H */
