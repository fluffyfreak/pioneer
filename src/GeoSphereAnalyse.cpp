// Copyright © 2008-2016 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "libs.h"
#include "GeoSphereAnalyse.h"
#include "GeoSphere.h"
#include <deque>
#include <algorithm>


#include "FileSystem.h"
#include "PngWriter.h"
#include "MathUtil.h"
#include "perlin.h"
#include "Cellular.h"
#include "Color.h"

double fbm(const vector3d &position, const int octaves, float frequency, const float persistence) 
{
	PROFILE_SCOPED()
	double total = 0.0;
	double maxAmplitude = 0.0;
	double amplitude = 1.0;
	for (int i = 0; i < octaves; i++) 
	{
		total += noise(position * frequency) * amplitude;
		frequency *= 2.0;
		maxAmplitude += amplitude;
		amplitude *= persistence;
	}
	return total / maxAmplitude;
}


static const double NORT = 1.0;
static const double EAST = 1.0;
static const double SOUT = -1.0;
static const double WEST = -1.0;
static const double NONE = 0.0;

static const Uint32 NUM_LAT_WINDS = 7;
static const Uint32 MAX_LAT_WINDS = NUM_LAT_WINDS-1;
static const vector2d latWinds[NUM_LAT_WINDS] = {
	vector2d(NONE,SOUT).Normalized(),
	vector2d(WEST,NONE).Normalized(),
	vector2d(EAST,NORT).Normalized(),
	vector2d(WEST,NONE).Normalized(),
	vector2d(EAST,SOUT).Normalized(),
	vector2d(WEST,NONE).Normalized(),
	vector2d(NONE,NORT).Normalized()
};

void ReadVector2dValues(std::vector<vector2d> &out, FileSystem::FileSource &fs, const std::string &path, const Uint32 scaleX, const Uint32 scaleY)
{
	PROFILE_SCOPED()
	RefCountedPtr<FileSystem::FileData> data = fs.ReadFile(path);
	StringRange buffer = data->AsStringRange();
	buffer = buffer.StripUTF8BOM();

	std::string section_name;

	while (!buffer.Empty()) 
	{
		StringRange line = buffer.ReadLine().StripSpace();

		// if the line is a comment, skip it
		if (line.Empty()) 
			continue;

		const char *kmid = line.FindChar(',');
		// if there's no '=' sign, skip the line
		if (kmid == line.end) {
			Output("WARNING: ignoring invalid line in file:\n   '%.*s'\n", int(line.Size()), line.begin);
			continue;
		}

		StringRange xval(line.begin, kmid);
		StringRange yval(kmid + 1, line.end);

		double x = atof( xval.ToString().c_str() );
		double y = atof( yval.ToString().c_str() );

		out.push_back(vector2d(x*scaleX,y*scaleY));
	}
}

void Analyse(GeoSphere *geo)
{
	PROFILE_SCOPED()

	// resolution of the maps we'll generate
	const Uint32 hmWide = 2048;//1024;//64;
	const Uint32 hmHigh = hmWide>>1;
	const double dWide = hmWide;
	const double dHigh = hmHigh;

	// generate sampling points
	std::vector<double> heights;
	std::vector<vector2d> poisson;
	poisson.reserve(6254);
	ReadVector2dValues(poisson, FileSystem::gameDataFiles, "Poisson.txt", hmWide-1, hmHigh-1);
	heights.resize(	poisson.size() );
	double minH=DBL_MAX;
	double maxH=DBL_MIN;
	for (size_t k = 0; k < poisson.size(); k++ ) 
	{
		PROFILE_SCOPED_DESC("poisson heights")
		const double wx = ((poisson[k].x / dWide) * 2.0) - 1.0;
		const double wy = ((poisson[k].y / dHigh) * 2.0) - 1.0;
		// 2D to polar
		const double lat = asin(wy);
		const double lon = -(wx * M_PI);
		// polar to 3D cartesian (normalised, no radius required/used)
		const double x = cos(lat) * cos(lon);
		const double z = cos(lat) * sin(lon);
		const double y = sin(lat);
		const double height = geo->GetHeight(vector3d(x,y,z).Normalized());
		minH = std::min(minH, height);
		maxH = std::max(maxH, height);
		heights[k] = height;
	}
	const double invMaxH = (is_equal_exact(maxH, 0.0)) ? 1.0 : (1.0 / maxH);
	for (size_t k = 0; k < heights.size(); k++ ) 
	{
		heights[k] *= invMaxH;
	}
	Cellular vor(16, hmWide, hmHigh, poisson, heights);

	// calculate storage for the maps
	const Uint32 bpp = 4; // channels of info... one byte per thing?
	const Uint32 stride = (bpp*hmWide + 3) & ~3;// pad rows to 4 bytes, which is the default row alignment for OpenGL
	std::unique_ptr<Uint8[]> pixels(new Uint8[stride * hmHigh]);
	std::unique_ptr<Uint8[]> pospixels(new Uint8[stride * hmHigh]);
	std::unique_ptr<vector2d[]> winds(new vector2d[hmWide * hmHigh]);
	std::unique_ptr<Uint8[]> windpixels(new Uint8[stride * hmHigh]);
	const double* heightmap = vor.CellularMap();

	// Calculate surface property maps
	Random rng(geo->GetSystemBody()->GetSeed()+4609837U);
	for(Uint32 h = 0; h<hmHigh; h++)
	{
		PROFILE_SCOPED_DESC("Surface properties Y")
		// normalised 2D coordinates (-1..1)
		const double wy = ((double(h) / dHigh) * 2.0) - 1.0;
		// 2D to polar
		const double lat = asin(wy);
		// cache sin & cos latitude values
		const double coslat = cos(lat);
		const double sinlat = sin(lat);
		// normalise from -rad..rad to 0..1 range
		const double normLat = (-lat+1.0) * 0.5;
		// Horizontal / longitude loop
		for(Uint32 w = 0; w<hmWide; w++)
		{
			PROFILE_SCOPED_DESC("Surface properties X")
			// normalised 2D coordinates (-1..1)
			const double wx = ((double(w) / dWide) * 2.0) - 1.0;
			// 2D to polar
			const double lon = -(wx * M_PI);
			// polar to 3D cartesian (normalised, no radius required/used)
			const double x = coslat * cos(lon);
			const double z = coslat * sin(lon);
			const double y = sinlat;
			const vector3d pos(x,y,z);

			// height at point
			const Uint32 hmIndex = (h * hmWide) + (w);
			const double height = heightmap[hmIndex];
			// adjust for water if necessary
			const bool bWater = (height<=0.0);

			// invert latitude since image start is lower-left not top-left
			const double distort = fbm(pos, 5, 2, 0.5) * 0.15;
			const double latWindsIndexD = MAX_LAT_WINDS*Clamp((1.0-normLat)+distort, 0.0, 1.0);
			const Uint32 latWindsIndexLower = Clamp(Uint32(floor(latWindsIndexD)), 0U, MAX_LAT_WINDS);
			const Uint32 latWindsIndexUpper = (latWindsIndexLower<MAX_LAT_WINDS) ? latWindsIndexLower+1 : latWindsIndexLower;
			// Add a little noise/jitter to the wind direction so that it's less uniform
			const vector2d jitter(vector2d((rng.Double()*2.0)-1.0, (rng.Double()*2.0)-1.0) * 0.05);
			const vector2d wind = (MathUtil::mix(latWinds[latWindsIndexLower], latWinds[latWindsIndexUpper], (latWindsIndexD - double(latWindsIndexLower))) + jitter).Normalized();
			winds[hmIndex] = wind;

			// calculate solar heating + height & water contributions
			const double intensity = Clamp((1.0 - abs(sinlat)) + distort, 0.0, 1.0);
			static const double ESun = 1.0;
			const double heatAbsorbtion = (((1.0-(bWater ? 0.5 : height))*intensity*intensity)+(intensity*0.5))*ESun;
			const double evaporation = bWater ? heatAbsorbtion : 0.0;

			// calculate friction
			static const double CLand = 0.9f;
			static const double CWater = 0.2f;
			const double friction = MathUtil::mix(CWater,CLand,height);
			
			const Uint32 index = (h * stride) + (w*bpp);
			{
				pixels[index + 0] = Uint8(Clamp(heatAbsorbtion * 255.0, 0.0, 255.0));
				pixels[index + 1] = Uint8(friction * 255);
				pixels[index + 2] = bWater ? Uint8(Clamp(evaporation * 255.0, 0.0, 255.0)) : 0;
				pixels[index + 3] = 255;

				pospixels[index + 0] = Uint8(((pos.x + 1.0) * 0.5) * 255.0);
				pospixels[index + 1] = Uint8(((pos.y + 1.0) * 0.5) * 255.0);
				pospixels[index + 2] = Uint8(((pos.z + 1.0) * 0.5) * 255.0);
				pospixels[index + 3] = 255;

				windpixels[index + 0] = Uint8(((1.0+wind.x)*0.5) * 255);
				windpixels[index + 1] = Uint8(((1.0+wind.y)*0.5) * 255);
				windpixels[index + 2] = 0;
				windpixels[index + 3] = 255;
			}
		}
	}

	// output images of what we've done
	if(true)
	{
		PROFILE_SCOPED_DESC("Output Images")

		const std::string dir = "planetmaps";
		FileSystem::userFiles.MakeDirectory(dir);
		char res[256];
		sprintf(res, "_%ux%u", hmWide, hmHigh);
		
		const std::string name( geo->GetSystemBody()->GetName() );
		{
			std::string destFile( name + std::string(res) + std::string(".png") );
			const std::string fname = FileSystem::JoinPathBelow(dir, destFile.c_str());
			write_png(FileSystem::userFiles, fname, pixels.get(), hmWide, hmHigh, stride, bpp);
		}
		{
			std::string destFile( name + std::string(res) + std::string("_pos.png") );
			const std::string fname = FileSystem::JoinPathBelow(dir, destFile.c_str());
			write_png(FileSystem::userFiles, fname, pospixels.get(), hmWide, hmHigh, stride, bpp);
		}
		{
			std::string destFile( name + std::string(res) + std::string("_winds.png") );
			const std::string fname = FileSystem::JoinPathBelow(dir, destFile.c_str());
			write_png(FileSystem::userFiles, fname, windpixels.get(), hmWide, hmHigh, stride, bpp);
		}
		std::unique_ptr<Uint8[]> voronoi(new Uint8[stride * hmHigh]);
		const double* buf = vor.CellularMap();
		for(int y=0; y<hmHigh; y++) 
		{
			for(int x=0; x<hmWide; x++) 
			{
				const double h = buf[(y * hmWide + x)];
				const Uint32 index = (y * stride) + (x*bpp);
				voronoi[index + 0] = h * 255;
				voronoi[index + 1] = h * 255;
				voronoi[index + 2] = h * 255;
				voronoi[index + 3] = 255;
			}
		}
		{
			std::string destFile( name + std::string(res) + std::string("_voronoi.png") );
			const std::string fname = FileSystem::JoinPathBelow(dir, destFile.c_str());
			write_png(FileSystem::userFiles, fname, voronoi.get(), hmWide, hmHigh, stride, bpp);
		}
	}
}

