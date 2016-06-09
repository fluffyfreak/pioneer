// Copyright © 2008-2016 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "libs.h"
#include "GeoSphere.h"
#include "GeoPatchContext.h"
#include "GeoPatch.h"
#include "GeoPatchJobs.h"
#include "perlin.h"
#include "Pi.h"
#include "RefCounted.h"
#include "graphics/Material.h"
#include "graphics/Renderer.h"
#include "graphics/Frustum.h"
#include "graphics/Graphics.h"
#include "graphics/Texture.h"
#include "graphics/TextureBuilder.h"
#include "graphics/VertexArray.h"
#include "vcacheopt/vcacheopt.h"
#include <deque>
#include <algorithm>

RefCountedPtr<GeoPatchContext> GeoSphere::s_patchContext;

// must be odd numbers
static const int detail_edgeLen[5] = {
	7, 15, 25, 35, 55
};

static const double gs_targetPatchTriLength(100.0);

#define PRINT_VECTOR(_v) Output("%f,%f,%f\n", (_v).x, (_v).y, (_v).z);

static const int geo_sphere_edge_friends[NUM_PATCHES][4] = {
	{ 3, 4, 1, 2 },
	{ 0, 4, 5, 2 },
	{ 0, 1, 5, 3 },
	{ 0, 2, 5, 4 },
	{ 0, 3, 5, 1 },
	{ 1, 4, 3, 2 }
};

static std::vector<GeoSphere*> s_allGeospheres;

void GeoSphere::Init()
{
	s_patchContext.Reset(new GeoPatchContext(detail_edgeLen[Pi::detail.planets > 4 ? 4 : Pi::detail.planets]));
}

void GeoSphere::Uninit()
{
	assert (s_patchContext.Unique());
	s_patchContext.Reset();
}

static void print_info(const SystemBody *sbody, const Terrain *terrain)
{
	Output(
		"%s:\n"
		"    height fractal: %s\n"
		"    colour fractal: %s\n"
		"    seed: %u\n",
		sbody->GetName().c_str(), terrain->GetHeightFractalName(), terrain->GetColorFractalName(), sbody->GetSeed());
}

// static
void GeoSphere::UpdateAllGeoSpheres()
{
	PROFILE_SCOPED()
	for(std::vector<GeoSphere*>::iterator i = s_allGeospheres.begin(); i != s_allGeospheres.end(); ++i)
	{
		(*i)->Update();
	}
}

// static
void GeoSphere::OnChangeDetailLevel()
{
	s_patchContext.Reset(new GeoPatchContext(detail_edgeLen[Pi::detail.planets > 4 ? 4 : Pi::detail.planets]));

	// reinit the geosphere terrain data
	for(std::vector<GeoSphere*>::iterator i = s_allGeospheres.begin(); i != s_allGeospheres.end(); ++i)
	{
		// clearout anything we don't need
		(*i)->Reset();

		// reinit the terrain with the new settings
		(*i)->m_terrain.Reset(Terrain::InstanceTerrain((*i)->GetSystemBody()));
		print_info((*i)->GetSystemBody(), (*i)->m_terrain.Get());
	}
}

//static
bool GeoSphere::OnAddQuadSplitResult(const SystemPath &path, SQuadSplitResult *res)
{
	// Find the correct GeoSphere via it's system path, and give it the split result
	for(std::vector<GeoSphere*>::iterator i=s_allGeospheres.begin(), iEnd=s_allGeospheres.end(); i!=iEnd; ++i) {
		if( path == (*i)->GetSystemBody()->GetPath() ) {
			(*i)->AddQuadSplitResult(res);
			return true;
		}
	}
	// GeoSphere not found to return the data to, cancel and delete it instead
	if( res ) {
		res->OnCancel();
		delete res;
	}
	return false;
}

//static
bool GeoSphere::OnAddSingleSplitResult(const SystemPath &path, SSingleSplitResult *res)
{
	// Find the correct GeoSphere via it's system path, and give it the split result
	for(std::vector<GeoSphere*>::iterator i=s_allGeospheres.begin(), iEnd=s_allGeospheres.end(); i!=iEnd; ++i) {
		if( path == (*i)->GetSystemBody()->GetPath() ) {
			(*i)->AddSingleSplitResult(res);
			return true;
		}
	}
	// GeoSphere not found to return the data to, cancel and delete it instead
	if( res ) {
		res->OnCancel();
		delete res;
	}
	return false;
}

void GeoSphere::Reset()
{
	{
		std::deque<SSingleSplitResult*>::iterator iter = mSingleSplitResults.begin();
		while(iter!=mSingleSplitResults.end())
		{
			// finally pass SplitResults
			SSingleSplitResult *psr = (*iter);
			assert(psr);

			psr->OnCancel();

			// tidyup
			delete psr;

			// Next!
			++iter;
		}
		mSingleSplitResults.clear();
	}

	{
		std::deque<SQuadSplitResult*>::iterator iter = mQuadSplitResults.begin();
		while(iter!=mQuadSplitResults.end())
		{
			// finally pass SplitResults
			SQuadSplitResult *psr = (*iter);
			assert(psr);

			psr->OnCancel();

			// tidyup
			delete psr;

			// Next!
			++iter;
		}
		mQuadSplitResults.clear();
	}

	for (int p=0; p<NUM_PATCHES; p++) {
		// delete patches
		if (m_patches[p]) {
			m_patches[p].reset();
		}
	}

	CalculateMaxPatchDepth();

	m_initStage = eBuildFirstPatches;
}

#define GEOSPHERE_TYPE	(GetSystemBody()->type)

GeoSphere::GeoSphere(const SystemBody *body) : BaseSphere(body),
	m_hasTempCampos(false), m_tempCampos(0.0), m_tempFrustum(800, 600, 0.5, 1.0, 1000.0),
	m_initStage(eBuildFirstPatches), m_maxDepth(0)
{
	print_info(body, m_terrain.Get());

	s_allGeospheres.push_back(this);

	CalculateMaxPatchDepth();

	//SetUpMaterials is not called until first Render since light count is zero :)
}

GeoSphere::~GeoSphere()
{
	// update thread should not be able to access us now, so we can safely continue to delete
	assert(std::count(s_allGeospheres.begin(), s_allGeospheres.end(), this) == 1);
	s_allGeospheres.erase(std::find(s_allGeospheres.begin(), s_allGeospheres.end(), this));
}

bool GeoSphere::AddQuadSplitResult(SQuadSplitResult *res)
{
	bool result = false;
	assert(res);
	assert(mQuadSplitResults.size()<MAX_SPLIT_OPERATIONS);
	if(mQuadSplitResults.size()<MAX_SPLIT_OPERATIONS) {
		mQuadSplitResults.push_back(res);
		result = true;
	}
	return result;
}

bool GeoSphere::AddSingleSplitResult(SSingleSplitResult *res)
{
	bool result = false;
	assert(res);
	assert(mSingleSplitResults.size()<MAX_SPLIT_OPERATIONS);
	if(mSingleSplitResults.size()<MAX_SPLIT_OPERATIONS) {
		mSingleSplitResults.push_back(res);
		result = true;
	}
	return result;
}

void GeoSphere::ProcessSplitResults()
{
	// now handle the single split results that define the base level of the quad tree
	{
		std::deque<SSingleSplitResult*>::iterator iter = mSingleSplitResults.begin();
		while(iter!=mSingleSplitResults.end())
		{
			// finally pass SplitResults
			SSingleSplitResult *psr = (*iter);
			assert(psr);

			const int32_t faceIdx = psr->face();
			if( m_patches[faceIdx] ) {
				m_patches[faceIdx]->ReceiveHeightmap(psr);
			} else {
				psr->OnCancel();
			}

			// tidyup
			delete psr;

			// Next!
			++iter;
		}
		mSingleSplitResults.clear();
	}

	// now handle the quad split results
	{
		std::deque<SQuadSplitResult*>::iterator iter = mQuadSplitResults.begin();
		while(iter!=mQuadSplitResults.end())
		{
			// finally pass SplitResults
			SQuadSplitResult *psr = (*iter);
			assert(psr);

			const int32_t faceIdx = psr->face();
			if( m_patches[faceIdx] ) {
				m_patches[faceIdx]->ReceiveHeightmaps(psr);
			} else {
				psr->OnCancel();
			}

			// tidyup
			delete psr;

			// Next!
			++iter;
		}
		mQuadSplitResults.clear();
	}
}

#include "FileSystem.h"
#include "PngWriter.h"
#include "MathUtil.h"
#include "perlin.h"
#include "Color.h"

inline int iwrap(int x, int y)
{
    if (x > 0)
        return x % y;
    if (x < 0)
        return (x + 1) % y + y - 1;
    return 0;
}

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

class Cellular {
public:
	Cellular(const int cellSize, const int dimx, const int dimy, const std::vector<vector2d> &s, const std::vector<double> &h) 
		: CELL_DIVISOR(cellSize), INV_CELL_DIVISOR(1.0/double(cellSize))
		, size_x(dimx), size_y(dimy), half_size_x(dimx>>1), half_size_y(dimy>>1)
		, cellsX(dimx / cellSize), cellsY(dimy / cellSize)
		, heights(h) 
	{
		PROFILE_SCOPED()
		cells.reset(new Cell[cellsX * cellsY]);
		const size_t count = s.size();
		for( size_t c=0; c<count; c++ )
		{
			// find which Cell to put the values into
			const Uint32 xi = s[c].x * INV_CELL_DIVISOR;
			const Uint32 yi = s[c].y * INV_CELL_DIVISOR;
			const Uint32 index = (yi * cellsX) + xi;
			cells[index].points.push_back(s[c]);
			cells[index].indices.push_back(c);
		}
		Output("Cellular :: CellsX %u, CellsY %u\n", cellsX, cellsY);

		GenMap();
	}

	double* CellularMap() const { return buf.get(); }
private:
	const int CELL_DIVISOR;
	const double INV_CELL_DIVISOR;

	const int size_x, size_y, half_size_x, half_size_y;
	const Uint32 cellsX, cellsY;
	const std::vector<double> &heights;
	std::unique_ptr<double[]> buf;

	class Cell 
	{
	public:
		Cell() {
			PROFILE_SCOPED()
			points.reserve(20);
			indices.reserve(20);
		}
		std::vector<vector2d>	points;
		std::vector<size_t>		indices;
	};
	std::unique_ptr<Cell[]> cells;
	
	inline double WrapDist( int x, int y, const vector2d &p) const
	{
		PROFILE_SCOPED()
		double dx = abs(x-p.x);
		double dy = abs(y-p.y);
		if (dx > half_size_x ) // only wrap on the horizontal
			dx = size_x-dx;
		//if (dy > half_size_y )
		//	dy = size_y-dy;
		// return squared distance
		return dx*dx + dy*dy;
	}

#define CELL_OFFSET 2
	size_t NearestSite( const int x, const int y ) const
	{
		PROFILE_SCOPED()
		size_t site=0xFFFFFFFF;
 		double mindist = DBL_MAX;

		// find this cell
		const int cellX = x * INV_CELL_DIVISOR;
		const int cellY = y * INV_CELL_DIVISOR;
		// search through the NxN grid of cells centred on cellX|cellY;
		for(int yi = std::max(cellY-CELL_OFFSET,0); yi<std::min(cellY+CELL_OFFSET,int(cellsY)); yi++) 
		{
			for(int xi = (cellX-CELL_OFFSET); xi<(cellX+CELL_OFFSET); xi++) 
			{
				const Uint32 index = (yi * cellsX) + iwrap(xi, cellsX);
				const std::vector<vector2d> &pts = cells[index].points;
				for (size_t i=0 ; i<pts.size() ; ++i)
				{
					double dist = WrapDist(x,y,pts[i]);
					if (dist < mindist) 
					{
						mindist = dist;
						site = cells[index].indices[i];
					}
				}
			}
		}
 		return site;
	}
	
	void GenMap()
	{
		PROFILE_SCOPED()
		buf.reset(new double[size_y * size_x]);
		double *ptr = buf.get();

		for (int i = 0; i < size_y; i++ ) {
			for (int j = 0; j < size_x; j++ ) {
				(*ptr) = heights[NearestSite(j, i)];
				++ptr;
			}
		}
	}
};

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
	const Uint32 hmWide = 1024;//64;
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
		const double lat = -asin(wy);
		const double lon = (wx * M_PI);
		// polar to 3D cartesian (normalised, no radius required/used)
		const double x = cos(lat) * cos(lon);
		const double y = cos(lat) * sin(lon);
		const double z = sin(lat);
		const double height = geo->GetHeight(vector3d(x,y,z));
		minH = std::min(minH, height);
		maxH = std::max(maxH, height);
		heights[k] = height;
	}
	const double invMaxH = (is_equal_exact(maxH, 0.0)) ? 1.0 : (1.0 / maxH);
	for (size_t k = 0; k < heights.size(); k++ ) 
	{
		PROFILE_SCOPED_DESC("poisson heights normalised")
		heights[k] *= invMaxH;
	}
	Cellular vor(8, hmWide, hmHigh, poisson, heights);

	// calculate storage for the maps
	const Uint32 bpp = 4; // channels of info... one byte per thing?
	const Uint32 stride = (bpp*hmWide + 3) & ~3;// pad rows to 4 bytes, which is the default row alignment for OpenGL
	std::unique_ptr<Uint8[]> pixels(new Uint8[stride * hmHigh]);
	std::unique_ptr<vector3d[]> spheremap(new vector3d[hmWide * hmHigh]);
	std::unique_ptr<vector2d[]> winds(new vector2d[hmWide * hmHigh]);
	std::unique_ptr<Uint8[]> windpixels(new Uint8[stride * hmHigh]);
	std::unique_ptr<Uint8[]> accumprecippixels(new Uint8[stride * hmHigh]);
	const double* heightmap = vor.CellularMap();

	// Calculate surface property maps
	Random rng(geo->GetSystemBody()->GetSeed()+4609837U);
	for(Uint32 h = 0; h<hmHigh; h++)
	{
		PROFILE_SCOPED_DESC("Surface properties Y")
		// normalised 2D coordinates (-1..1)
		const double wy = ((double(h) / dHigh) * 2.0) - 1.0;
		// 2D to polar
		const double lat = -asin(wy);
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
			const double lon = (wx * M_PI);
			// polar to 3D cartesian (normalised, no radius required/used)
			const double x = coslat * cos(lon);
			const double y = coslat * sin(lon);
			const double z = sinlat;
			const vector3d pos(x,y,z);

			// height at point
			const Uint32 hmIndex = (h * hmWide) + (w);
			const double height = heightmap[hmIndex] * invMaxH;
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

				windpixels[index + 0] = Uint8(((1.0+wind.x)*0.5) * 255);
				windpixels[index + 1] = Uint8(((1.0+wind.y)*0.5) * 255);
				windpixels[index + 2] = 0;
				windpixels[index + 3] = 255;
			}
		}
	}

	std::unique_ptr<double[]> accumulation(new double[hmWide * hmHigh]);
	std::unique_ptr<double[]> precipitation(new double[hmWide * hmHigh]);
	memset(accumulation.get(), 0, sizeof(double) * hmWide * hmHigh);
	memset(precipitation.get(), 0, sizeof(double) * hmWide * hmHigh);

	// accumulate evaporation
	double minAccum = DBL_MAX;
	double maxAccum = DBL_MIN;
	double minPrecip = DBL_MAX;
	double maxPrecip = DBL_MIN;
	//for(Uint32 simLoop = 0; simLoop<250; simLoop++)
	{
		PROFILE_SCOPED_DESC("simLoop")
		for(Uint32 h = 0; h<hmHigh; h++)
		{
			// Horizontal / longitude loop
			for(Uint32 w = 0; w<hmWide; w++)
			{
				PROFILE_SCOPED_DESC("simLoop inner")
				// indexes
				const Uint32 hmIndex = (h * hmWide) + (w);
				const Uint32 pixIndex = (h * stride) + (w*bpp);

				// values
				const double heat = pixels[pixIndex + 0] / 255.0;
				const double fric = pixels[pixIndex + 1] / 255.0;
				const double evap = pixels[pixIndex + 2] / 255.0;

				// evaporation
				const double accumEvap = accumulation[hmIndex] + evap; // add more water from this source

				// move the water around
				const vector2d &windDir = winds[hmIndex];
				const Sint32 newW = (Sint32(w) + (windDir.x < 0.0 ? -1 : 1));
				const Uint32 xmov = iwrap(newW, (hmWide-1));
				const Uint32 ymov = Clamp(h + (windDir.y < 0.0 ? -1 : 1), 0U, hmHigh-1U);

				const Uint32 xIndex = (h * hmWide) + xmov;
				const Uint32 yIndex = (ymov * hmWide) + w;
				const Uint32 xyIndex = (ymov * hmWide) + xmov;
				const double xPcnt = 1.0 - abs(windDir.x);
				const double yPcnt = 1.0 - abs(windDir.y);
				const double xyPcnt = 1.0 - (xPcnt + yPcnt);
				// actually updating the values
				accumulation[xIndex]	+= accumEvap * xPcnt;
				accumulation[yIndex]	+= accumEvap * yPcnt;
				accumulation[xyIndex]	+= accumEvap * xyPcnt;

				// calculate precipitation for this cell
				const double precip = (fric + (1.0-heat)) * accumEvap;
				// Fade the current amount of precipitation a bit when adding on the new amount
				//precipitation[hmIndex] = Clamp((precipitation[hmIndex]*0.5) + precip, 0.0, 1.0);
				precipitation[hmIndex] = Clamp(precip, 0.0, 1.0);
				minPrecip = std::min(minPrecip, precipitation[hmIndex]);
				maxPrecip = std::max(maxPrecip, precipitation[hmIndex]);
				accumulation[hmIndex] = Clamp(accumulation[hmIndex] - precipitation[hmIndex], 0.0, 1.0);
				minAccum = std::min(minAccum, accumulation[hmIndex]);
				maxAccum = std::max(maxAccum, accumulation[hmIndex]);
			}
		}
	}
	const double invMaxPrecip = 1.0 / maxPrecip;
	const double invMaxAccum = 1.0 / maxAccum;

	for(Uint32 h = 0; h<hmHigh; h++)
	{
		// Horizontal / longitude loop
		for(Uint32 w = 0; w<hmWide; w++)
		{
			// indexes
			const Uint32 hmIndex = (h * hmWide) + (w);
			const Uint32 pixIndex = (h * stride) + (w*bpp);

			// evaporation
			accumprecippixels[pixIndex + 0] = Uint8(Clamp((accumulation[hmIndex]*invMaxAccum) * 255.0, 0.0, 255.0));
			accumprecippixels[pixIndex + 1] = Uint8(Clamp((precipitation[hmIndex]*invMaxPrecip) * 255.0, 0.0, 255.0));
			accumprecippixels[pixIndex + 2] = 0;
			accumprecippixels[pixIndex + 3] = 255;
		}
	}

	// output images of what we've done
	if(false)
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
			std::string destFile( name + std::string(res) + std::string("_winds.png") );
			const std::string fname = FileSystem::JoinPathBelow(dir, destFile.c_str());
			write_png(FileSystem::userFiles, fname, windpixels.get(), hmWide, hmHigh, stride, bpp);
		}
		{
			std::string destFile( name + std::string(res) + std::string("_precip.png") );
			const std::string fname = FileSystem::JoinPathBelow(dir, destFile.c_str());
			write_png(FileSystem::userFiles, fname, accumprecippixels.get(), hmWide, hmHigh, stride, bpp);
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

void GeoSphere::BuildFirstPatches()
{
	assert(!m_patches[0]);
	if(m_patches[0])
		return;

	CalculateMaxPatchDepth();

	// generate root face patches of the cube/sphere
	static const vector3d p1 = (vector3d( 1, 1, 1)).Normalized();
	static const vector3d p2 = (vector3d(-1, 1, 1)).Normalized();
	static const vector3d p3 = (vector3d(-1,-1, 1)).Normalized();
	static const vector3d p4 = (vector3d( 1,-1, 1)).Normalized();
	static const vector3d p5 = (vector3d( 1, 1,-1)).Normalized();
	static const vector3d p6 = (vector3d(-1, 1,-1)).Normalized();
	static const vector3d p7 = (vector3d(-1,-1,-1)).Normalized();
	static const vector3d p8 = (vector3d( 1,-1,-1)).Normalized();

	const uint64_t maxShiftDepth = GeoPatchID::MAX_SHIFT_DEPTH;

	m_patches[0].reset(new GeoPatch(s_patchContext, this, p1, p2, p3, p4, 0, (0ULL << maxShiftDepth)));
	m_patches[1].reset(new GeoPatch(s_patchContext, this, p4, p3, p7, p8, 0, (1ULL << maxShiftDepth)));
	m_patches[2].reset(new GeoPatch(s_patchContext, this, p1, p4, p8, p5, 0, (2ULL << maxShiftDepth)));
	m_patches[3].reset(new GeoPatch(s_patchContext, this, p2, p1, p5, p6, 0, (3ULL << maxShiftDepth)));
	m_patches[4].reset(new GeoPatch(s_patchContext, this, p3, p2, p6, p7, 0, (4ULL << maxShiftDepth)));
	m_patches[5].reset(new GeoPatch(s_patchContext, this, p8, p7, p6, p5, 0, (5ULL << maxShiftDepth)));

	for (int i=0; i<NUM_PATCHES; i++) {
		m_patches[i]->RequestSinglePatch();
	}

	const std::string name( GetSystemBody()->GetName() );
	if(name == "New Hope")
	{
#ifdef PIONEER_PROFILER
		Profiler::reset();
		FileSystem::userFiles.MakeDirectory("GeoSphere");
		const std::string profilerPath = FileSystem::JoinPathBelow(FileSystem::userFiles.GetRoot(), "GeoSphere");
#endif
		{
			PROFILE_SCOPED_DESC("GeoSphere-Analyse")
			Analyse(this);
		}
#ifdef PIONEER_PROFILER
		Profiler::dumphtml(profilerPath.c_str());
#endif
	}

	m_initStage = eRequestedFirstPatches;
}

void GeoSphere::CalculateMaxPatchDepth()
{
	const double circumference = 2.0 * M_PI * m_sbody->GetRadius();
	// calculate length of each edge segment (quad) times 4 due to that being the number around the sphere (1 per side, 4 sides for Root).
	double edgeMetres = circumference / double(s_patchContext->GetEdgeLen() * 8);
	// find out what depth we reach the desired resolution
	while (edgeMetres>gs_targetPatchTriLength && m_maxDepth<GEOPATCH_MAX_DEPTH) {
		edgeMetres *= 0.5;
		++m_maxDepth;
	}
}

void GeoSphere::Update()
{
	switch(m_initStage)
	{
	case eBuildFirstPatches:
		BuildFirstPatches();
		break;
	case eRequestedFirstPatches:
		{
			ProcessSplitResults();
			uint8_t numValidPatches = 0;
			for (int i=0; i<NUM_PATCHES; i++) {
				if(m_patches[i]->HasHeightData()) {
					++numValidPatches;
				}
			}
			m_initStage = (NUM_PATCHES==numValidPatches) ? eReceivedFirstPatches : eRequestedFirstPatches;
		} break;
	case eReceivedFirstPatches:
		{
			for (int i=0; i<NUM_PATCHES; i++) {
				m_patches[i]->NeedToUpdateVBOs();
			}
			m_initStage = eDefaultUpdateState;
		} break;
	case eDefaultUpdateState:
		if(m_hasTempCampos) {
			ProcessSplitResults();
			for (int i=0; i<NUM_PATCHES; i++) {
				m_patches[i]->LODUpdate(m_tempCampos, m_tempFrustum);
			}
			ProcessQuadSplitRequests();
		}
		break;
	}
}

void GeoSphere::AddQuadSplitRequest(double dist, SQuadSplitRequest *pReq, GeoPatch *pPatch)
{
	mQuadSplitRequests.push_back(TDistanceRequest(dist, pReq, pPatch));
}

void GeoSphere::ProcessQuadSplitRequests()
{
	class RequestDistanceSort {
	public:
		bool operator()(const TDistanceRequest &a, const TDistanceRequest &b)
		{
			return a.mDistance < b.mDistance;
		}
	};
	std::sort(mQuadSplitRequests.begin(), mQuadSplitRequests.end(), RequestDistanceSort());

	for(auto iter : mQuadSplitRequests) {
		SQuadSplitRequest *ssrd = iter.mpRequest;
		iter.mpRequester->ReceiveJobHandle(Pi::GetAsyncJobQueue()->Queue(new QuadPatchJob(ssrd)));
	}
	mQuadSplitRequests.clear();
}

void GeoSphere::Render(Graphics::Renderer *renderer, const matrix4x4d &modelView, vector3d campos, const float radius, const std::vector<Camera::Shadow> &shadows)
{
	PROFILE_SCOPED()
	// store this for later usage in the update method.
	m_tempCampos = campos;
	m_hasTempCampos = true;

	if(m_initStage < eDefaultUpdateState)
		return;

	matrix4x4d trans = modelView;
	trans.Translate(-campos.x, -campos.y, -campos.z);
	renderer->SetTransform(trans); //need to set this for the following line to work
	matrix4x4d modv;
	matrix4x4d proj;
	matrix4x4ftod(renderer->GetCurrentModelView(), modv);
	matrix4x4ftod(renderer->GetCurrentProjection(), proj);
	Graphics::Frustum frustum( modv, proj );
	m_tempFrustum = frustum;

	// no frustum test of entire geosphere, since Space::Render does this
	// for each body using its GetBoundingRadius() value

	//First draw - create materials (they do not change afterwards)
	if (!m_surfaceMaterial)
		SetUpMaterials();

	bool bHasAtmosphere = false;
	{
		//Update material parameters
		//XXX no need to calculate AP every frame
		m_materialParameters.atmosphere = GetSystemBody()->CalcAtmosphereParams();
		m_materialParameters.atmosphere.center = trans * vector3d(0.0, 0.0, 0.0);
		m_materialParameters.atmosphere.planetRadius = radius;

		m_materialParameters.shadows = shadows;

		m_materialParameters.maxPatchDepth = GetMaxDepth();

		m_surfaceMaterial->specialParameter0 = &m_materialParameters;
		
		bHasAtmosphere = (m_materialParameters.atmosphere.atmosDensity > 0.0);
		if (bHasAtmosphere) {
			m_atmosphereMaterial->specialParameter0 = &m_materialParameters;

			// make atmosphere sphere slightly bigger than required so
			// that the edges of the pixel shader atmosphere jizz doesn't
			// show ugly polygonal angles
			DrawAtmosphereSurface(renderer, trans, campos,
				m_materialParameters.atmosphere.atmosRadius*1.01,
				m_atmosRenderState, m_atmosphereMaterial );
		}
	}

	// display the terrain height control-mesh
	
	const float rad = m_materialParameters.atmosphere.atmosRadius * 0.99f;
	const matrix4x4d cloudsTrans(trans * matrix4x4d::ScaleMatrix(rad, rad, rad));
	if(bHasAtmosphere)
	{
		// bunny, ball ball!
		if( !m_cloudSphere.get() ) {
			if(!m_cloudMaterial.Valid()) {
				Graphics::MaterialDescriptor matDesc;
				matDesc.effect = Graphics::EFFECT_CLOUD_SPHERE;
				matDesc.textures = 2;
				m_cloudMaterial.Reset(Pi::renderer->CreateMaterial(matDesc));
				m_cloudMaterial->diffuse = Color4f(0.7f, 0.7f, 0.7f, 0.5f);
				m_cloudMaterial->texture0 = Graphics::TextureBuilder::Raw("textures/permTexture.png").GetOrCreateTexture(Pi::renderer, "noise");
				m_cloudMaterial->texture1 = Graphics::TextureBuilder::Raw("textures/gradTexture.png").GetOrCreateTexture(Pi::renderer, "noise");
			}

			//blended
			Graphics::RenderStateDesc rsd;
			rsd.blendMode = Graphics::BLEND_ALPHA;
			rsd.depthWrite = false;
			rsd.cullMode = Graphics::CULL_NONE;
			m_cloudSphere.reset( new Graphics::Drawables::Sphere3D(Pi::renderer, m_cloudMaterial, Pi::renderer->CreateRenderState(rsd), 5, 1.0) );
		}
		m_cloudMaterial->specialParameter0 = &m_materialParameters;
	}

	Color ambient;
	Color &emission = m_surfaceMaterial->emissive;

	// save old global ambient
	const Color oldAmbient = renderer->GetAmbientColor();

	if ((GetSystemBody()->GetSuperType() == SystemBody::SUPERTYPE_STAR) || (GetSystemBody()->GetType() == SystemBody::TYPE_BROWN_DWARF)) {
		// stars should emit light and terrain should be visible from distance
		ambient.r = ambient.g = ambient.b = 51;
		ambient.a = 255;
		emission = StarSystem::starRealColors[GetSystemBody()->GetType()];
		emission.a = 255;
	}

	else {
		// give planet some ambient lighting if the viewer is close to it
		double camdist = campos.Length();
		camdist = 0.1 / (camdist*camdist);
		// why the fuck is this returning 0.1 when we are sat on the planet??
		// JJ: Because campos is relative to a unit-radius planet - 1.0 at the surface
		// XXX oh well, it is the value we want anyway...
		ambient.r = ambient.g = ambient.b = camdist * 255;
		ambient.a = 255;
	}

	renderer->SetAmbientColor(ambient);

	renderer->SetTransform(modelView);

	for (int i=0; i<NUM_PATCHES; i++) {
		m_patches[i]->Render(renderer, campos, modelView, frustum);
	}

	renderer->SetAmbientColor(oldAmbient);
	
	if( bHasAtmosphere && m_cloudSphere )
	{
		// draw it
		renderer->SetTransform(cloudsTrans);
		m_cloudSphere->Draw( Pi::renderer );
	}

	renderer->GetStats().AddToStatCount(Graphics::Stats::STAT_PLANETS, 1);
}

void GeoSphere::SetUpMaterials()
{
	//solid
	Graphics::RenderStateDesc rsd;
	m_surfRenderState = Pi::renderer->CreateRenderState(rsd);

	//blended
	rsd.blendMode = Graphics::BLEND_ALPHA_ONE;
	rsd.cullMode = Graphics::CULL_NONE;
	rsd.depthWrite = false;
	m_atmosRenderState = Pi::renderer->CreateRenderState(rsd);

	// Request material for this star or planet, with or without
	// atmosphere. Separate material for surface and sky.
	Graphics::MaterialDescriptor surfDesc;
	const Uint32 effect_flags = m_terrain->GetSurfaceEffects();
	if (effect_flags & Terrain::EFFECT_LAVA)
		surfDesc.effect = Graphics::EFFECT_GEOSPHERE_TERRAIN_WITH_LAVA;
	else if (effect_flags & Terrain::EFFECT_WATER)
		surfDesc.effect = Graphics::EFFECT_GEOSPHERE_TERRAIN_WITH_WATER;
	else
		surfDesc.effect = Graphics::EFFECT_GEOSPHERE_TERRAIN;

	if ((GetSystemBody()->GetType() == SystemBody::TYPE_BROWN_DWARF) ||
		(GetSystemBody()->GetType() == SystemBody::TYPE_STAR_M)) {
		//dim star (emits and receives light)
		surfDesc.lighting = true;
		surfDesc.quality &= ~Graphics::HAS_ATMOSPHERE;
	}
	else if (GetSystemBody()->GetSuperType() == SystemBody::SUPERTYPE_STAR) {
		//normal star
		surfDesc.lighting = false;
		surfDesc.quality &= ~Graphics::HAS_ATMOSPHERE;
		surfDesc.effect = Graphics::EFFECT_GEOSPHERE_STAR;
	} else {
		//planetoid with or without atmosphere
		const SystemBody::AtmosphereParameters ap(GetSystemBody()->CalcAtmosphereParams());
		surfDesc.lighting = true;
		if(ap.atmosDensity > 0.0) {
			surfDesc.quality |= Graphics::HAS_ATMOSPHERE;
		} else {
			surfDesc.quality &= ~Graphics::HAS_ATMOSPHERE;
		}
	}

	const bool bEnableEclipse = (Pi::config->Int("DisableEclipse") == 0);
	if (bEnableEclipse) {
		surfDesc.quality |= Graphics::HAS_ECLIPSES;
	}
	const bool bEnableDetailMaps = (Pi::config->Int("DisableDetailMaps") == 0);
	if (bEnableDetailMaps) {
		surfDesc.quality |= Graphics::HAS_DETAIL_MAPS;
	}
	m_surfaceMaterial.Reset(Pi::renderer->CreateMaterial(surfDesc));

	m_texHi.Reset( Graphics::TextureBuilder::Model("textures/high.dds").GetOrCreateTexture(Pi::renderer, "model") );
	m_texLo.Reset( Graphics::TextureBuilder::Model("textures/low.dds").GetOrCreateTexture(Pi::renderer, "model") );
	m_surfaceMaterial->texture0 = m_texHi.Get();
	m_surfaceMaterial->texture1 = m_texLo.Get();

	{
		Graphics::MaterialDescriptor skyDesc;
		skyDesc.effect = Graphics::EFFECT_GEOSPHERE_SKY;
		skyDesc.lighting = true;
		if (bEnableEclipse) {
			skyDesc.quality |= Graphics::HAS_ECLIPSES;
		}
		m_atmosphereMaterial.Reset(Pi::renderer->CreateMaterial(skyDesc));
		m_atmosphereMaterial->texture0 = nullptr;
		m_atmosphereMaterial->texture1 = nullptr;
	}
}
