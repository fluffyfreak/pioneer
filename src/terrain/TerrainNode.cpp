// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Terrain.h"
#include "TerrainNode.h"
#include "TerrainNoise.h"

#include <math.h>
#include "MathUtil.h"

#include "Easing.h"

#include "FileSystem.h"
#include "json/json.h"

#include "vector2.h"

using namespace TerrainNoise;

namespace cellywelly
{
	using std::floor;
	using std::sqrt;
	using std::min;
	using std::max;
	// Cellular noise ("Worley noise") in 3D in GLSL.
	// Copyright (c) Stefan Gustavson 2011-04-19. All rights reserved.
	// This code is released under the conditions of the MIT license.
	// See LICENSE file for details.

	vector2d min(const vector2d& a, const vector2d& b) {
		return vector2d(
			min(a.x, b.x),
			min(a.y, b.y));
	}

	vector3d min(const vector3d& a, const vector3d& b) {
		return vector3d(
			min(a.x, b.x),
			min(a.y, b.y),
			min(a.z, b.z));
	}

	vector3d max(const vector3d& a, const vector3d& b) {
		return vector3d(
			max(a.x, b.x),
			max(a.y, b.y),
			max(a.z, b.z));
	}

	vector3d floor(const vector3d& rhs) {
		return vector3d(floor(rhs.x), floor(rhs.y), floor(rhs.z));
	}

	// https://www.opengl.org/sdk/docs/man/html/fract.xhtml
	// fract returns the fractional part of x. This is calculated as x - floor(x).
	vector3d fract(const vector3d& rhs) {
		return vector3d(rhs.x - floor(rhs.x), rhs.y - floor(rhs.y), rhs.z - floor(rhs.z));
	}

	// sqrt returns the square root of x. i.e., the value sqrt(x).
	vector2d sqrt(const vector2d& rhs) {
		return vector2d(sqrt(rhs.x), sqrt(rhs.y));
	}

	// sqrt returns the square root of x. i.e., the value sqrt(x).
	vector3d sqrt(const vector3d& rhs) {
		return vector3d(sqrt(rhs.x), sqrt(rhs.y), sqrt(rhs.z));
	}

	// https://www.opengl.org/sdk/docs/man/html/inversesqrt.xhtml
	// inversesqrt returns the inverse of the square root of x. i.e., the value 1/sqrt(x).
	vector3d inversesqrt(const vector3d& rhs) {
		return vector3d(1.0 / sqrt(rhs.x), 1.0 / sqrt(rhs.y), 1.0 / sqrt(rhs.z));
	}

	// x - y * floor(x / y) for each component in x using the floating point value y.
	vector3d mod(const vector3d& rhs, const double y) {
		return vector3d(
			rhs.x - y * floor(rhs.x / y),
			rhs.y - y * floor(rhs.y / y),
			rhs.z - y * floor(rhs.z / y));
	}

	// Permutation polynomial: (34x^2 + x) mod 289
	vector3d permute(vector3d x) {
		return mod((34.0 * x + 1.0) * x, 289.0);
	}

	// Cellular noise, returning F1 and F2 in a vec2.
	// 3x3x3 search region for good F2 everywhere, but a lot
	// slower than the 2x2x2 version.
	// The code below is a bit scary even to its author,
	// but it has at least half decent performance on a
	// modern GPU. In any case, it beats any software
	// implementation of Worley noise hands down.

	vector2d cellular(const vector3d &P) {
#define K 0.142857142857 // 1/7
#define Ko 0.428571428571 // 1/2-K/2
#define K2 0.020408163265306 // 1/(7*7)
#define Kz 0.166666666667 // 1/6
#define Kzo 0.416666666667 // 1/2-1/6*2
#define jitter 1.0 // smaller jitter gives more regular pattern

		vector3d Pi = mod(floor(P), 289.0);
		vector3d Pf = fract(P) - 0.5;

		vector3d Pfx = Pf.x + vector3d(1.0, 0.0, -1.0);
		vector3d Pfy = Pf.y + vector3d(1.0, 0.0, -1.0);
		vector3d Pfz = Pf.z + vector3d(1.0, 0.0, -1.0);

		vector3d p = permute(Pi.x + vector3d(-1.0, 0.0, 1.0));
		vector3d p1 = permute(p + Pi.y - 1.0);
		vector3d p2 = permute(p + Pi.y);
		vector3d p3 = permute(p + Pi.y + 1.0);

		vector3d p11 = permute(p1 + Pi.z - 1.0);
		vector3d p12 = permute(p1 + Pi.z);
		vector3d p13 = permute(p1 + Pi.z + 1.0);

		vector3d p21 = permute(p2 + Pi.z - 1.0);
		vector3d p22 = permute(p2 + Pi.z);
		vector3d p23 = permute(p2 + Pi.z + 1.0);

		vector3d p31 = permute(p3 + Pi.z - 1.0);
		vector3d p32 = permute(p3 + Pi.z);
		vector3d p33 = permute(p3 + Pi.z + 1.0);

		vector3d ox11 = fract(p11*K) - Ko;
		vector3d oy11 = mod(floor(p11*K), 7.0)*K - Ko;
		vector3d oz11 = floor(p11*K2)*Kz - Kzo; // p11 < 289 guaranteed

		vector3d ox12 = fract(p12*K) - Ko;
		vector3d oy12 = mod(floor(p12*K), 7.0)*K - Ko;
		vector3d oz12 = floor(p12*K2)*Kz - Kzo;

		vector3d ox13 = fract(p13*K) - Ko;
		vector3d oy13 = mod(floor(p13*K), 7.0)*K - Ko;
		vector3d oz13 = floor(p13*K2)*Kz - Kzo;

		vector3d ox21 = fract(p21*K) - Ko;
		vector3d oy21 = mod(floor(p21*K), 7.0)*K - Ko;
		vector3d oz21 = floor(p21*K2)*Kz - Kzo;

		vector3d ox22 = fract(p22*K) - Ko;
		vector3d oy22 = mod(floor(p22*K), 7.0)*K - Ko;
		vector3d oz22 = floor(p22*K2)*Kz - Kzo;

		vector3d ox23 = fract(p23*K) - Ko;
		vector3d oy23 = mod(floor(p23*K), 7.0)*K - Ko;
		vector3d oz23 = floor(p23*K2)*Kz - Kzo;

		vector3d ox31 = fract(p31*K) - Ko;
		vector3d oy31 = mod(floor(p31*K), 7.0)*K - Ko;
		vector3d oz31 = floor(p31*K2)*Kz - Kzo;

		vector3d ox32 = fract(p32*K) - Ko;
		vector3d oy32 = mod(floor(p32*K), 7.0)*K - Ko;
		vector3d oz32 = floor(p32*K2)*Kz - Kzo;

		vector3d ox33 = fract(p33*K) - Ko;
		vector3d oy33 = mod(floor(p33*K), 7.0)*K - Ko;
		vector3d oz33 = floor(p33*K2)*Kz - Kzo;

		vector3d dx11 = Pfx + jitter*ox11;
		vector3d dy11 = Pfy.x + jitter*oy11;
		vector3d dz11 = Pfz.x + jitter*oz11;

		vector3d dx12 = Pfx + jitter*ox12;
		vector3d dy12 = Pfy.x + jitter*oy12;
		vector3d dz12 = Pfz.y + jitter*oz12;

		vector3d dx13 = Pfx + jitter*ox13;
		vector3d dy13 = Pfy.x + jitter*oy13;
		vector3d dz13 = Pfz.z + jitter*oz13;

		vector3d dx21 = Pfx + jitter*ox21;
		vector3d dy21 = Pfy.y + jitter*oy21;
		vector3d dz21 = Pfz.x + jitter*oz21;

		vector3d dx22 = Pfx + jitter*ox22;
		vector3d dy22 = Pfy.y + jitter*oy22;
		vector3d dz22 = Pfz.y + jitter*oz22;

		vector3d dx23 = Pfx + jitter*ox23;
		vector3d dy23 = Pfy.y + jitter*oy23;
		vector3d dz23 = Pfz.z + jitter*oz23;

		vector3d dx31 = Pfx + jitter*ox31;
		vector3d dy31 = Pfy.z + jitter*oy31;
		vector3d dz31 = Pfz.x + jitter*oz31;

		vector3d dx32 = Pfx + jitter*ox32;
		vector3d dy32 = Pfy.z + jitter*oy32;
		vector3d dz32 = Pfz.y + jitter*oz32;

		vector3d dx33 = Pfx + jitter*ox33;
		vector3d dy33 = Pfy.z + jitter*oy33;
		vector3d dz33 = Pfz.z + jitter*oz33;

		vector3d d11 = dx11 * dx11 + dy11 * dy11 + dz11 * dz11;
		vector3d d12 = dx12 * dx12 + dy12 * dy12 + dz12 * dz12;
		vector3d d13 = dx13 * dx13 + dy13 * dy13 + dz13 * dz13;
		vector3d d21 = dx21 * dx21 + dy21 * dy21 + dz21 * dz21;
		vector3d d22 = dx22 * dx22 + dy22 * dy22 + dz22 * dz22;
		vector3d d23 = dx23 * dx23 + dy23 * dy23 + dz23 * dz23;
		vector3d d31 = dx31 * dx31 + dy31 * dy31 + dz31 * dz31;
		vector3d d32 = dx32 * dx32 + dy32 * dy32 + dz32 * dz32;
		vector3d d33 = dx33 * dx33 + dy33 * dy33 + dz33 * dz33;

		// Sort out the two smallest distances (F1, F2)
#if 0
		// Cheat and sort out only F1
		vector3d d1 = min(min(d11, d12), d13);
		vector3d d2 = min(min(d21, d22), d23);
		vector3d d3 = min(min(d31, d32), d33);
		vector3d d = min(min(d1, d2), d3);
		d.x = min(min(d.x, d.y), d.z);
		return sqrt(d.xx); // F1 duplicated, no F2 computed
#else
		// Do it right and sort out both F1 and F2
		vector3d d1a = min(d11, d12);
		d12 = max(d11, d12);
		d11 = min(d1a, d13); // Smallest now not in d12 or d13
		d13 = max(d1a, d13);
		d12 = min(d12, d13); // 2nd smallest now not in d13
		vector3d d2a = min(d21, d22);
		d22 = max(d21, d22);
		d21 = min(d2a, d23); // Smallest now not in d22 or d23
		d23 = max(d2a, d23);
		d22 = min(d22, d23); // 2nd smallest now not in d23
		vector3d d3a = min(d31, d32);
		d32 = max(d31, d32);
		d31 = min(d3a, d33); // Smallest now not in d32 or d33
		d33 = max(d3a, d33);
		d32 = min(d32, d33); // 2nd smallest now not in d33
		vector3d da = min(d11, d21);
		d21 = max(d11, d21);
		d11 = min(da, d31); // Smallest now in d11
		d31 = max(da, d31); // 2nd smallest now not in d31
		d11.xy((d11.x < d11.y) ? d11.xy() : d11.yx());
		d11.xz((d11.x < d11.z) ? d11.xz() : d11.zx()); // d11.x now smallest
		d12 = min(d12, d21); // 2nd smallest now not in d21
		d12 = min(d12, d22); // nor in d22
		d12 = min(d12, d31); // nor in d31
		d12 = min(d12, d32); // nor in d32
		d11.yz(min(d11.yz(), d12.xy())); // nor in d12.yz
		d11.y = min(d11.y, d12.z); // Only two more to go
		d11.y = min(d11.y, d11.z); // Done! (Phew!)
		return sqrt(d11.xy()); // F1, F2
#endif
	}

	//varying vector3d vTexCoord3D;
	//
	//void main(void) {
	//	vec2 F = cellular2x2x2(vTexCoord3D.xyz);
	//	float n = smoothstep(0.4, 0.5, F.x);
	//	gl_FragColor = vector4d(n, n, n, 1.0);
	//}
}

inline double noise(const int octaves, double frequency, const double persistence, const vector3d &position) {
	double total = 0.0;
	double maxAmplitude = 0.0;
	double amplitude = 1.0;
	for (int i = 0; i < octaves; i++) {
		total += noise(position * frequency) * amplitude;
		frequency *= 2.0;
		maxAmplitude += amplitude;
		amplitude *= persistence;
	}
	return total / maxAmplitude;
}

inline double noise_cubed(const int octaves, double frequency, const double persistence, const vector3d &position) {
	double total = 0.0;
	double maxAmplitude = 0.0;
	double amplitude = 1.0;
	for (int i = 0; i < octaves; i++) {
		total += noise(position * frequency) * amplitude;
		frequency *= 2.0;
		maxAmplitude += amplitude;
		amplitude *= persistence;
	}
	return pow(total / maxAmplitude, 3.0);
}

inline double noise_ridged(const int octaves, double frequency, const double persistence, const vector3d &position) {
	double total = 0.0;
	double maxAmplitude = 0.0;
	double amplitude = 1.0;
	for (int i = 0; i < octaves; i++) {
		total += ((1.0 - abs(noise(position * frequency))) * 2.0 - 1.0) * amplitude;
		frequency *= 2.0;
		maxAmplitude += amplitude;
		amplitude *= persistence;
	}
	return total / maxAmplitude;
}

double noise_cellular_squared(const int octaves, double frequency, const double persistence, const vector3d &position) {
	double total = 0.0;
	double maxAmplitude = 0.0;
	double amplitude = 1.0;
	double tmp;
	for (int i = 0; i < octaves; i++) {
		vector2d ff = cellywelly::cellular(position * frequency);
		tmp = ff.y - ff.x;
		total += tmp * tmp * amplitude;
		frequency *= 2.0;
		maxAmplitude += amplitude;
		amplitude *= persistence;
	}
	return total / maxAmplitude;
}


double TerrainNodeData::Call(const vector3d& p, const vector3d& pNorm)
{
	// A good way to fix this is to use another low frequency noise function that is scaled between 0 and 1.
	// If we multiply the output of that noise function with our mountain function, then mountains will exists where it is 1, but not where it is 0.
	// ...
	// Notice that weve added a few things here.
	// We moved the mountain function to be a child of our new function. 
	// "op: mul" will multiply the output of the function with its direct children, and the "clamp" will clamp the output between 0 and 1.

	// process the child nodes as we depend on what values they return
	double childH = 0.0;
	for (auto child = m_children.begin(), end = m_children.end(); child!=end; ++child)
	{
		childH += child->Call(p, pNorm);
	}

	// now our local value
	double localH = 0.0;
	switch (m_nodeType)
	{
		// noise is always scaled and clamped
	case NT_NOISE:					localH = Clamp(Scale(noise(m_octaves, m_frequency, m_persistence, p)));						break;
	case NT_NOISE_CELLULAR_SQUARED:	localH = Clamp(Scale(noise_cellular_squared(m_octaves, m_frequency, m_persistence, p)));	break;
	case NT_NOISE_RIDGED:			localH = Clamp(Scale(noise_ridged(m_octaves, m_frequency, m_persistence, p)));				break;
	case NT_NOISE_CUBED:			localH = Clamp(Scale(noise_cubed(m_octaves, m_frequency, m_persistence, p)));				break;
		// heightmaps are a data source not a noise time
	case NT_HEIGHTMAP:				localH = Clamp(GetHeightMapValue(pNorm));													break;
	}

	// mix the values together
	switch (m_op)
	{
	case TO_ADD: localH += childH; break;
	case TO_SUB: localH -= childH; break;
	case TO_MUL: localH *= childH; break;
	case TO_DIV: localH /= childH; break;
	}

	return localH;
}

static size_t bufread_or_die(void *ptr, size_t size, size_t nmemb, ByteRange &buf)
{
	size_t read_count = buf.read(static_cast<char*>(ptr), size, nmemb);
	if (read_count < nmemb) {
		Output("Error: failed to read file (truncated)\n");
		abort();
	}
	return read_count;
}

bool TerrainNodeData::LoadHeightmap(const std::string &filename)
{
	if (!filename.empty()) 
	{
		RefCountedPtr<FileSystem::FileData> fdata = FileSystem::gameDataFiles.ReadFile(FileSystem::JoinPathBelow("heightmaps", filename));
		if (!fdata) {
			Output("Error: could not open file '%s'\n", filename.c_str());
			abort();
		}

		ByteRange databuf = fdata->AsByteRange();

		Sint16 minHMap = INT16_MAX, maxHMap = INT16_MIN;
		Uint16 minHMapScld = UINT16_MAX, maxHMapScld = 0;

		// XXX unify heightmap types
		Uint16 v;
		bufread_or_die(&v, 2, 1, databuf); m_heightMapSizeX = v;
		bufread_or_die(&v, 2, 1, databuf); m_heightMapSizeY = v;
		const Uint32 heightmapPixelArea = (m_heightMapSizeX * m_heightMapSizeY);

		std::unique_ptr<Sint16[]> heightMap(new Sint16[heightmapPixelArea]);
		bufread_or_die(heightMap.get(), sizeof(Sint16), heightmapPixelArea, databuf);
		assert(m_heightMap==nullptr);
		m_heightMap = new double[heightmapPixelArea];
		double *pHeightMap = m_heightMap;
		for(Uint32 i=0; i<heightmapPixelArea; i++) {
			const Sint16 val = heightMap.get()[i];
			minHMap = std::min(minHMap, val);
			maxHMap = std::max(maxHMap, val);
			// store then increment pointer
			(*pHeightMap) = val;
			++pHeightMap;
		}
		assert(pHeightMap == &m_heightMap[heightmapPixelArea]);
		//Output("minHMap = (%hd), maxHMap = (%hd)\n", minHMap, maxHMap);

		/*case 1: {
			Uint16 v;
			// XXX x and y reversed from above *sigh*
			bufread_or_die(&v, 2, 1, databuf); m_heightMapSizeY = v;
			bufread_or_die(&v, 2, 1, databuf); m_heightMapSizeX = v;
			const Uint32 heightmapPixelArea = (m_heightMapSizeX * m_heightMapSizeY);

			// read height scaling and min height which are doubles
			double te;
			bufread_or_die(&te, 8, 1, databuf);
			m_heightScaling = te;
			bufread_or_die(&te, 8, 1, databuf);
			m_minh = te;

			std::unique_ptr<Uint16[]> heightMapScaled(new Uint16[heightmapPixelArea]);
			bufread_or_die(heightMapScaled.get(), sizeof(Uint16), heightmapPixelArea, databuf);
			m_heightMap.reset(new double[heightmapPixelArea]);
			double *pHeightMap = m_heightMap.get();
			for(Uint32 i=0; i<heightmapPixelArea; i++) {
				const Uint16 val = heightMapScaled[i];
				minHMapScld = std::min(minHMapScld, val);
				maxHMapScld = std::max(maxHMapScld, val);
				// store then increment pointer
				(*pHeightMap) = val;
				++pHeightMap;
			}
			assert(pHeightMap == &m_heightMap[heightmapPixelArea]);
			//Output("minHMapScld = (%hu), maxHMapScld = (%hu)\n", minHMapScld, maxHMapScld);
			break;
		}*/
		return true;
	}
	return false;
}

double TerrainNodeData::GetHeightMapValue(const vector3d& p)
{
	if(!m_heightMap)
		return 0.0;

	double latitude = -asin(p.y);
	if (p.y < -1.0) latitude = -0.5*M_PI;
	if (p.y > 1.0) latitude = 0.5*M_PI;
//	if (!isfinite(latitude)) {
//		// p.y is just n of asin domain [-1,1]
//		latitude = (p.y < 0 ? -0.5*M_PI : M_PI*0.5);
//	}
	double longitude = atan2(p.x, p.z);
	double px = (((m_heightMapSizeX-1) * (longitude + M_PI)) / (2*M_PI));
	double py = ((m_heightMapSizeY-1)*(latitude + 0.5*M_PI)) / M_PI;
	int ix = int(floor(px));
	int iy = int(floor(py));
	ix = ::Clamp(ix, 0, m_heightMapSizeX-1);
	iy = ::Clamp(iy, 0, m_heightMapSizeY-1);
	double dx = px-ix;
	double dy = py-iy;

	// p0,3 p1,3 p2,3 p3,3
	// p0,2 p1,2 p2,2 p3,2
	// p0,1 p1,1 p2,1 p3,1
	// p0,0 p1,0 p2,0 p3,0
	double map[4][4];
	const double *pHMap = m_heightMap;
	for (int x=-1; x<3; x++) {
		for (int y=-1; y<3; y++) {
			map[x+1][y+1] = pHMap[::Clamp(iy+y, 0, m_heightMapSizeY-1)*m_heightMapSizeX + ::Clamp(ix+x, 0, m_heightMapSizeX-1)];
		}
	}

	double c[4];
	for (int j=0; j<4; j++) {
		double d0 = map[0][j] - map[1][j];
		double d2 = map[2][j] - map[1][j];
		double d3 = map[3][j] - map[1][j];
		double a0 = map[1][j];
		double a1 = -(1/3.0)*d0 + d2 - (1/6.0)*d3;
		double a2 = 0.5*d0 + 0.5*d2;
		double a3 = -(1/6.0)*d0 - 0.5*d2 + (1/6.0)*d3;
		c[j] = a0 + a1*dx + a2*dx*dx + a3*dx*dx*dx;
	}

	double h;
	{
		double d0 = c[0] - c[1];
		double d2 = c[2] - c[1];
		double d3 = c[3] - c[1];
		double a0 = c[1];
		double a1 = -(1/3.0)*d0 + d2 - (1/6.0)*d3;
		double a2 = 0.5*d0 + 0.5*d2;
		double a3 = -(1/6.0)*d0 - 0.5*d2 + (1/6.0)*d3;
		double v = a0 + a1*dy + a2*dy*dy + a3*dy*dy*dy;

		//v = (v<0 ? 0 : v); // limit to 0.0???
		h = v;
	}
	return h;
}


void ParseTerrainNode(Json::Value::iterator& j, TerrainNodeData& node)
{
	bool accum = false;
	double tempHigh = 0.0;
	double tempLow = 0.0;
	for (Json::Value::iterator funcTag = (*j).begin(), funcTagEnd = (*j).end(); funcTag != funcTagEnd; ++funcTag)
	{
		const std::string tag = funcTag.key().asString();
		if (tag == "children")
		{
			const bool bIsArray = (*funcTag).isArray();
			for (Json::Value::iterator childIT = (*funcTag).begin(), childEnd = (*funcTag).end(); childIT != childEnd; ++childIT)
			{
				TerrainNodeData child;
				ParseTerrainNode(childIT, child);
				node.AddChild(child);
			}
		}
		else if (tag == "clamp")
		{
			const bool bIsArray = (*funcTag).isArray();
			assert(bIsArray);
			const double low = (*funcTag).get(Json::ArrayIndex(0), 0).asDouble();
			const double high = (*funcTag).get(Json::ArrayIndex(1), 0).asDouble();
			node.ClampNoise(low, high);
		}
		else if (tag == "frequency")
		{
			node.Frequency((*funcTag).asDouble());
		}
		else if (tag == "scale")
		{
			const bool bIsArray = (*funcTag).isArray();
			assert(bIsArray);
			const double low = (*funcTag).get(Json::ArrayIndex(0), 0).asDouble();
			const double high = (*funcTag).get(Json::ArrayIndex(1), 0).asDouble();
			node.Scale(low, high);
		}
		else if (tag == "name")
		{
			node.Name((*funcTag).asString());
		}
		else if (tag == "octaves")
		{
			node.Octaves((*funcTag).asInt());
		}
		else if (tag == "op")
		{
			const std::string opStr = (*funcTag).asString();
			if (opStr == "add") { node.Op(TerrainNodeData::TO_ADD); }
			else if (opStr == "sub") { node.Op(TerrainNodeData::TO_SUB); }
			else if (opStr == "mul") { node.Op(TerrainNodeData::TO_MUL); }
			else if (opStr == "div") { node.Op(TerrainNodeData::TO_DIV); }
		}
		else if (tag == "persistence")
		{
			node.Persistence((*funcTag).asDouble());
		}
		else if (tag == "type")
		{
			node.NoiseType((*funcTag).asString());
		}
		else if (tag == "file")
		{
			node.LoadHeightmap((*funcTag).asString());
		}
		Output("\t\ttag:\"%s\"\n", tag.c_str());
	}
}

void LoadTerrainJSON(const std::string& path, std::vector<TerrainSource>& sources)
{
	Json::Reader reader;
	Json::Value data;

	auto fd = FileSystem::gameDataFiles.ReadFile(path);
	if (!fd) {
		Output("couldn't open json terrain definition '%s'\n", path.c_str());
		return;
	}

	if (!reader.parse(fd->GetData(), fd->GetData() + fd->GetSize(), data)) {
		Output("couldn't read json terrain definition '%s': %s\n", path.c_str(), reader.getFormattedErrorMessages().c_str());
		return;
	}

	Output("\n%s\n", path.c_str());
	for (Json::Value::iterator slot = data.begin(); slot != data.end(); ++slot)
	{
		TerrainSource source;
		const std::string key = slot.key().asString();
		if (key == "baseHeight") {
			source.SetType(TerrainSource::ST_HEIGHT);
		}
		else if (key == "humidity") {
			source.SetType(TerrainSource::ST_HUMIDITY);
		}
		else if (key == "temperature") {
			source.SetType(TerrainSource::ST_TEMPERATURE);
		}
		else if (key == "jitter") {
			source.SetType(TerrainSource::ST_JITTER);
		}

		// get base height
		const double baseHeight = (*slot).get("base", 0).asDouble();
		source.SetBaseHeight(baseHeight);

		Json::Value funcs = (*slot).get("funcs", Json::arrayValue);
		if (!funcs.empty())
		{
			for (Json::Value::iterator j = funcs.begin(), jEnd = funcs.end(); j != jEnd; ++j)
			{
				Json::Value name = (*j).get("name", 0);
				if (!name.empty() && name.isString()) {
					const std::string func(name.asString());
					Output("\tfunc:\"%s\"\n", func.c_str());
				}

				TerrainNodeData tn;
				ParseTerrainNode(j, tn);
				source.AddNode(tn);
			}
		}
		sources.push_back(source);
	}
}