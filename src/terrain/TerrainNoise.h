// Copyright © 2008-2016 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _TERRAINNOISE_H
#define _TERRAINNOISE_H

#include "Terrain.h"
#include "perlin.h"

namespace TerrainNoise {

	// octavenoise functions return range [0,1] if persistence = 0.5
	inline double octavenoise(const fracdef_t &def, const double persistence, const vector3d &p) {
		//assert(persistence <= (1.0 / def.lacunarity));
		double n = 0;
		double amplitude = persistence;
		double maxAmplitude = 0.0;
		double frequency = def.frequency;
		for (int i=0; i<def.octaves; i++) {
			n += amplitude * noise(frequency*p);
			maxAmplitude += amplitude;
			amplitude *= persistence;
			frequency *= def.lacunarity;
		}
		n /= maxAmplitude;
		return (n+1.0)*0.5;
	}

	inline double river_octavenoise(const fracdef_t &def, const double persistence, const vector3d &p) {
		//assert(persistence <= (1.0 / def.lacunarity));
		double n = 0;
		double amplitude = persistence;
		double maxAmplitude = 0.0;
		double frequency = def.frequency;
		for (int i=0; i<def.octaves; i++) {
			n += amplitude * fabs(noise(frequency*p));
			maxAmplitude += amplitude;
			amplitude *= persistence;
			frequency *= def.lacunarity;
		}
		n /= maxAmplitude;
		return fabs(n);
	}

	inline double ridged_octavenoise(const fracdef_t &def, const double persistence, const vector3d &p) {
		//assert(persistence <= (1.0 / def.lacunarity));
		double n = 0;
		double amplitude = persistence;
		double maxAmplitude = 0.0;
		double frequency = def.frequency;
		for (int i=0; i<def.octaves; i++) {
			n += amplitude * noise(frequency*p);
			maxAmplitude += amplitude;
			amplitude *= persistence;
			frequency *= def.lacunarity;
		}
		n /= maxAmplitude;
		n = 1.0 - fabs(n);
		n *= n;
		return n;
	}

	inline double billow_octavenoise(const fracdef_t &def, const double persistence, const vector3d &p) {
		//assert(persistence <= (1.0 / def.lacunarity));
		double n = 0;
		double amplitude = persistence;
		double maxAmplitude = 0.0;
		double frequency = def.frequency;
		for (int i=0; i<def.octaves; i++) {
			n += amplitude * noise(frequency*p);
			maxAmplitude += amplitude;
			amplitude *= persistence;
			frequency *= def.lacunarity;
		}
		n /= maxAmplitude;
		return (2.0 * fabs(n) - 1.0)+1.0;
	}

	inline double voronoiscam_octavenoise(const fracdef_t &def, const double persistence, const vector3d &p) {
		//assert(persistence <= (1.0 / def.lacunarity));
		double n = 0;
		double amplitude = persistence;
		double maxAmplitude = 0.0;
		double frequency = def.frequency;
		for (int i=0; i<def.octaves; i++) {
			n += amplitude * noise(frequency*p);
			maxAmplitude += amplitude;
			amplitude *= persistence;
			frequency *= def.lacunarity;
		}
		n /= maxAmplitude;
		return sqrt(10.0 * fabs(n));
	}

	inline double dunes_octavenoise(const fracdef_t &def, const double persistence, const vector3d &p) {
		//assert(persistence <= (1.0 / def.lacunarity));
		double n = 0;
		double amplitude = persistence;
		double maxAmplitude = 0.0;
		double frequency = def.frequency;
		for (int i=0; i<3; i++) {
			n += amplitude * noise(frequency*p);
			maxAmplitude += amplitude;
			amplitude *= persistence;
			frequency *= def.lacunarity;
		}
		n /= maxAmplitude;
		return 1.0 - fabs(n);
	}

	// XXX merge these with their fracdef versions
	inline double octavenoise(int octaves, const double persistence, const double lacunarity, const vector3d &p) {
		//assert(persistence <= (1.0 / lacunarity));
		double n = 0;
		double amplitude = persistence;
		double maxAmplitude = 0.0;
		double frequency = 1.0;
		while (octaves--) {
			n += amplitude * noise(frequency*p);
			maxAmplitude += amplitude;
			amplitude *= persistence;
			frequency *= lacunarity;
		}
		n /= maxAmplitude;
		return (n+1.0)*0.5;
	}

	inline double river_octavenoise(int octaves, const double persistence, const double lacunarity, const vector3d &p) {
		//assert(persistence <= (1.0 / lacunarity));
		double n = 0;
		double amplitude = persistence;
		double maxAmplitude = 0.0;
		double frequency = 1.0;
		while (octaves--) {
			n += amplitude * fabs(noise(frequency*p));
			maxAmplitude += amplitude;
			amplitude *= persistence;
			frequency *= lacunarity;
		}
		n /= maxAmplitude;
		return n;
	}

	inline double ridged_octavenoise(int octaves, const double persistence, const double lacunarity, const vector3d &p) {
		//assert(persistence <= (1.0 / lacunarity));
		double n = 0;
		double amplitude = persistence;
		double maxAmplitude = 0.0;
		double frequency = 1.0;
		while (octaves--) {
			n += amplitude * noise(frequency*p);
			maxAmplitude += amplitude;
			amplitude *= persistence;
			frequency *= lacunarity;
		}
		n /= maxAmplitude;
		n = 1.0 - fabs(n);
		n *= n;
		return n;
	}

	inline double billow_octavenoise(int octaves, const double persistence, const double lacunarity, const vector3d &p) {
		//assert(persistence <= (1.0 / lacunarity));
		double n = 0;
		double amplitude = persistence;
		double maxAmplitude = 0.0;
		double frequency = 1.0;
		while (octaves--) {
			n += amplitude * noise(frequency*p);
			maxAmplitude += amplitude;
			amplitude *= persistence;
			frequency *= lacunarity;
		}
		n /= maxAmplitude;
		return (2.0 * fabs(n) - 1.0)+1.0;
	}

	inline double voronoiscam_octavenoise(int octaves, const double persistence, const double lacunarity, const vector3d &p) {
		//assert(persistence <= (1.0 / lacunarity));
		double n = 0;
		double amplitude = persistence;
		double maxAmplitude = 0.0;
		double frequency = 1.0;
		while (octaves--) {
			n += amplitude * noise(frequency*p);
			maxAmplitude += amplitude;
			amplitude *= persistence;
			frequency *= lacunarity;
		}
		n /= maxAmplitude;
		return sqrt(10.0 * fabs(n));
	}
}

#endif
