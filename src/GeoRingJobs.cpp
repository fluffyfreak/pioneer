#include "libs.h"
#include "GeoPatchContext.h"
#include "GeoPatchID.h"
#include "GeoRing.h"
#include "GeoRingJobs.h"
#include "perlin.h"
#include "Pi.h"
#include "galaxy/StarSystem.h"
#include "graphics/Frustum.h"
#include "graphics/Renderer.h"
#include "graphics/Material.h"
#include "graphics/Texture.h"
#include "graphics/TextureBuilder.h"

#include "profiler/Profiler.h"

#include "vcacheopt/vcacheopt.h"

#pragma pack(4)
struct VBOVertex
{
	vector3f pos;
	vector3f norm;
	Color4ub col;
	vector2f uv;
};
#pragma pack()

#define VBO_COUNT_HI_EDGE  (3*(edgeLen-1))
#define VBO_COUNT_MID_IDX  (4*3*(edgeLen-3) + 2*(edgeLen-3)*(edgeLen-3)*3)
//                          ^^ serrated teeth bit      ^^^ square inner bit

#define lerp(t, a, b) ( a + t * (b - a) )

// hold the 16 possible terrain edge connections
typedef struct {
	std::vector<unsigned short> v;
} VecShortHolder;

inline void SetColour(Color3ub *r, const vector3d &v) { 
#if 1
	r->r=static_cast<unsigned char>(Clamp(v.x*255.0, 0.0, 255.0)); 
	r->g=static_cast<unsigned char>(Clamp(v.y*255.0, 0.0, 255.0)); 
	r->b=static_cast<unsigned char>(Clamp(v.z*255.0, 0.0, 255.0));
#else
	r->r = 255;
	r->g = 255;
	r->b = 255;
#endif
}

// ********************************************************************************
//
// ********************************************************************************
#define BORDER_SIZE 1

inline vector3d GetSurfacePointCyl(const double x, const double y, const double ang0, const double ang1, const double yoffset, const double halfLength) 
{
	double theta = lerp( x, ang1, ang0 );
		
	const vector3d topEndEdge(sin(theta), yoffset + (halfLength * 0.5), cos(theta));		// vertices at top edge of circle
	const vector3d bottomEndEdge(sin(theta), yoffset - (halfLength * 0.5), cos(theta));	// vertices at bottom edge of circle
		
	const vector3d res = lerp( y, bottomEndEdge, topEndEdge );
	return res;
}

SBaseSplitRequest::SBaseSplitRequest(const double ang0_, const double ang1_, const double yoffset_, const double halfLen_, const vector3d &vb0_, const vector3d &vb1_,
	const uint32_t depth_, const SystemPath &sysPath_, const GeoPlateID &patchID_, const int edgeLen_, const double fracStep_,
	Terrain *pTerrain_)
	: ang0(ang0_), ang1(ang1_), yoffset(yoffset_), halfLen(halfLen_), vbe0(vb0_), vbe1(vb1_), depth(depth_), 
	sysPath(sysPath_), patchID(patchID_), edgeLen(edgeLen_), fracStep(fracStep_), 
	pTerrain(pTerrain_)
{
}


SQuadPlateRequest::SQuadPlateRequest(const double ang0_, const double ang1_, const double yoffset_, const double halfLen_, const vector3d &vb0_, const vector3d &vb1_,
	const uint32_t depth_, const SystemPath &sysPath_, const GeoPlateID &patchID_, const int edgeLen_, const double fracStep_,
	Terrain *pTerrain_)
	: SBaseSplitRequest(ang0_, ang1_, yoffset_, halfLen_, vb0_, vb1_, depth_, sysPath_, patchID_, edgeLen_, fracStep_, pTerrain_)
{
	const int numVerts = NUMVERTICES(edgeLen_);
	for( int i=0 ; i<4 ; ++i )
	{
		heights[i] = new double[numVerts];
		normals[i] = new vector3f[numVerts];
		colors[i] = new Color3ub[numVerts];
	}
	const int numBorderedVerts = NUMVERTICES((edgeLen_*2)+(BORDER_SIZE*2)-1);
	borderHeights.reset(new double[numBorderedVerts]);
	borderVertexs.reset(new vector3d[numBorderedVerts]);
}


SSinglePlateRequest::SSinglePlateRequest(const double ang0_, const double ang1_, const double yoffset_, const double halfLen_, const vector3d &vb0_, const vector3d &vb1_,
	const uint32_t depth_, const SystemPath &sysPath_, const GeoPlateID &patchID_, const int edgeLen_, const double fracStep_,
	Terrain *pTerrain_)
	: SBaseSplitRequest(ang0_, ang1_, yoffset_, halfLen_, vb0_, vb1_, depth_, sysPath_, patchID_, edgeLen_, fracStep_, pTerrain_)
{
	const int numVerts = NUMVERTICES(edgeLen_);
	heights = new double[numVerts];
	normals = new vector3f[numVerts];
	colors = new Color3ub[numVerts];
		
	const int numBorderedVerts = NUMVERTICES(edgeLen_+(BORDER_SIZE*2));
	borderHeights.reset(new double[numBorderedVerts]);
	borderVertexs.reset(new vector3d[numBorderedVerts]);
}

// Generates full-detail vertices, and also non-edge normals and colors 
void SinglePlateJob::GenerateMesh(double *heights, vector3f *normals, Color3ub *colors, 
								double *borderHeights, vector3d *borderVertexs,
								const double ang0,
								const double ang1,
								const double yoffset,
								const double halfLen,
								const int edgeLen,
								const double fracStep,
								const Terrain *pTerrain) const
{
	const int borderedEdgeLen = edgeLen+(BORDER_SIZE*2);
	const int numBorderedVerts = borderedEdgeLen*borderedEdgeLen;

	const vector3d topEnd(0.0, yoffset + halfLen, 0.0);		// vertices at top edge of circle
	const vector3d btmEnd(0.0, yoffset - halfLen, 0.0);		// vertices at bottom edge of circle

	// generate heights plus a 1 unit border
	double *bhts = borderHeights;
	vector3d *vrts = borderVertexs;
	for (int y=-BORDER_SIZE; y<borderedEdgeLen-BORDER_SIZE; y++) 
	{
		const double yfrac = double(y) * fracStep;
		// point along z-axis by zfrac amount
		const vector3d axisPt = lerp( yfrac, btmEnd, topEnd );
		for (int x=-BORDER_SIZE; x<borderedEdgeLen-BORDER_SIZE; x++) 
		{
			const double xfrac = double(x) * fracStep;
			const vector3d p = GetSurfacePointCyl(xfrac, yfrac, ang0, ang1, yoffset, halfLen);
			// vector from axis to point-on-surface
			const vector3d cDir = (p - axisPt).Normalized();
			const double height = pTerrain->GetHeight(p);
			assert(height >= 0.0f && height <= 1.0f);
			*(bhts++) = -height;
			*(vrts++) = p + (cDir * height);//p * (height + 1.0);
		}
	}
	assert(bhts==&borderHeights[numBorderedVerts]);

	// Generate normals & colors for non-edge vertices since they never change
	Color3ub *col = colors;
	vector3f *nrm = normals;
	double *hts = heights;
	vrts = borderVertexs;
	for (int y=BORDER_SIZE; y<borderedEdgeLen-BORDER_SIZE; y++) {
		for (int x=BORDER_SIZE; x<borderedEdgeLen-BORDER_SIZE; x++) {
			// height
			const double height = borderHeights[x + y*borderedEdgeLen];
			assert(hts!=&heights[edgeLen*edgeLen]);
			*(hts++) = height;

			// normal
			const vector3d &x1 = vrts[(x-1) + y*borderedEdgeLen];
			const vector3d &x2 = vrts[(x+1) + y*borderedEdgeLen];
			const vector3d &y1 = vrts[x + (y-1)*borderedEdgeLen];
			const vector3d &y2 = vrts[x + (y+1)*borderedEdgeLen];
			const vector3d n = ((x2-x1).Cross(y2-y1)).Normalized();
			assert(nrm!=&normals[edgeLen*edgeLen]);
			*(nrm++) = vector3f(n);

			// color
			const vector3d p = GetSurfacePointCyl((x-BORDER_SIZE)*fracStep, (y-BORDER_SIZE)*fracStep, ang0, ang1, yoffset, halfLen);
			SetColour(col, pTerrain->GetColor(p, height, n));
			assert(col!=&colors[edgeLen*edgeLen]);
			++col;
		}
	}
	assert(hts==&heights[edgeLen*edgeLen]);
	assert(nrm==&normals[edgeLen*edgeLen]);
	assert(col==&colors[edgeLen*edgeLen]);
}

// ********************************************************************************
// Overloaded PureJob class to handle generating the mesh for each patch
// ********************************************************************************
void SinglePlateJob::OnFinish()  // runs in primary thread of the context
{
	GeoRing::OnAddSinglePlateResult( mData->sysPath, mpResults );
	mpResults = nullptr;
	BasePlateJob::OnFinish();
}

void SinglePlateJob::OnRun()    // RUNS IN ANOTHER THREAD!! MUST BE THREAD SAFE!
{
	BasePlateJob::OnRun();

	const SSinglePlateRequest &srd = *mData;

	// fill out the data
	GenerateMesh(srd.heights, srd.normals, srd.colors, srd.borderHeights.get(), srd.borderVertexs.get(),
		srd.ang0, srd.ang1, srd.yoffset, srd.halfLen, 
		srd.edgeLen, srd.fracStep, srd.pTerrain.Get());
	// add this patches data
	SSinglePlateResult *sr = new SSinglePlateResult(srd.patchID.GetPlateIdx(), srd.depth);
	sr->addResult(srd.heights, srd.normals, srd.colors, 
		srd.ang0, srd.ang1, srd.yoffset, srd.halfLen, srd.vbe0, srd.vbe1,
		srd.patchID.NextPatchID(srd.depth+1, 0));
	// store the result
	mpResults = sr;
}

SinglePlateJob::~SinglePlateJob()
{
	if(mpResults) {
		mpResults->OnCancel();
		delete mpResults;
		mpResults = nullptr;
	}
}

// ********************************************************************************
// Overloaded PureJob class to handle generating the mesh for each patch
// ********************************************************************************
void QuadPlateJob::OnFinish()  // runs in primary thread of the context
{
	GeoRing::OnAddQuadPlateResult( mData->sysPath, mpResults );
	mpResults = nullptr;
	BasePlateJob::OnFinish();
}

void QuadPlateJob::OnRun()    // RUNS IN ANOTHER THREAD!! MUST BE THREAD SAFE!
{
	BasePlateJob::OnRun();

	const SQuadPlateRequest &srd = *mData;

	GenerateBorderedData(srd.borderHeights.get(), srd.borderVertexs.get(),
			srd.ang0, srd.ang1, srd.yoffset, srd.halfLen,
			srd.edgeLen, srd.fracStep, srd.pTerrain.Get());

	
	const vector3d halfVBE = srd.vbe0 + 0.5 * (srd.vbe1 - srd.vbe0);
	const double halfAng = srd.ang0 + 0.5 * (srd.ang1 - srd.ang0);
	const double newLength = (srd.halfLen*0.5);
	const double yoffset0 = srd.yoffset + (newLength*0.5);
	const double yoffset1 = srd.yoffset - (newLength*0.5);
	
	const double vecs[4][4] = {
		{halfAng, srd.ang1, yoffset1, newLength},
		{srd.ang0, halfAng, yoffset1, newLength},
		{srd.ang0, halfAng, yoffset0, newLength},
		{halfAng, srd.ang1, yoffset0, newLength}
	};

	const vector3d vbes[4][2] = {
		{halfVBE, srd.vbe1},
		{srd.vbe0, halfVBE},
		{srd.vbe0, halfVBE},
		{halfVBE, srd.vbe1}
	};

	const int borderedEdgeLen = (srd.edgeLen*2)+(BORDER_SIZE*2)-1;
	const int offxy[4][2] = {
		{0,0},
		{srd.edgeLen-1,0},
		{srd.edgeLen-1,srd.edgeLen-1},
		{0,srd.edgeLen-1}
	};

	SQuadPlateResult *sr = new SQuadPlateResult(srd.patchID.GetPlateIdx(), srd.depth);
	for (int i=0; i<4; i++)
	{
		// fill out the data
		GenerateSubPatchData(srd.heights[i], srd.normals[i], srd.colors[i], srd.borderHeights.get(), srd.borderVertexs.get(),
			vecs[i][0], vecs[i][1], vecs[i][2], vecs[i][3], 
			srd.edgeLen, offxy[i][0], offxy[i][1], 
			borderedEdgeLen, srd.fracStep, srd.pTerrain.Get());

		// add this patches data
		sr->addResult(i, srd.heights[i], srd.normals[i], srd.colors[i], 
			vecs[i][0], vecs[i][1], vecs[i][2], vecs[i][3], vbes[i][0], vbes[i][1],
			srd.patchID.NextPatchID(srd.depth+1, i));
	}
	mpResults = sr;
}

QuadPlateJob::~QuadPlateJob()
{
	if(mpResults) {
		mpResults->OnCancel();
		delete mpResults;
		mpResults = NULL;
	}
}

// Generates full-detail vertices, and also non-edge normals and colors 
void QuadPlateJob::GenerateBorderedData(
	double *borderHeights, vector3d *borderVertexs,
	const double ang0,
	const double ang1,
	const double yoffset,
	const double halfLen,
	const int edgeLen,
	const double fracStep,
	const Terrain *pTerrain) const
{
	const int borderedEdgeLen = (edgeLen * 2) + (BORDER_SIZE * 2) - 1;
	const int numBorderedVerts = borderedEdgeLen*borderedEdgeLen;

	const vector3d topEnd(0.0, yoffset + halfLen, 0.0);		// vertices at top edge of circle
	const vector3d btmEnd(0.0, yoffset - halfLen, 0.0);		// vertices at bottom edge of circle

	// generate heights plus a N=BORDER_SIZE unit border
	double *bhts = borderHeights;
	vector3d *vrts = borderVertexs;
	for ( int y = -BORDER_SIZE; y < (borderedEdgeLen - BORDER_SIZE); y++ ) 
	{
		const double yfrac = double(y) * (fracStep*0.5);
		// point along z-axis by zfrac amount
		const vector3d axisPt = lerp( yfrac, btmEnd, topEnd );
		for ( int x = -BORDER_SIZE; x < (borderedEdgeLen - BORDER_SIZE); x++ ) 
		{
			const double xfrac = double(x) * fracStep;
			const vector3d p = GetSurfacePointCyl(xfrac, yfrac, ang0, ang1, yoffset, halfLen);
			// vector from axis to point-on-surface
			const vector3d cDir = (p - axisPt).Normalized();
			const double height = pTerrain->GetHeight(p);
			assert(height >= 0.0f && height <= 1.0f);
			*(bhts++) = -height;
			*(vrts++) = p + (cDir * height);
		}
	}
	assert(bhts == &borderHeights[numBorderedVerts]);
}

void QuadPlateJob::GenerateSubPatchData(
	double *heights, vector3f *normals, Color3ub *colors, 
	double *borderHeights, vector3d *borderVertexs,
	const double ang0,
	const double ang1,
	const double yoffset,
	const double halfLen,
	const int edgeLen,
	const int xoff, 
	const int yoff,
	const int borderedEdgeLen,
	const double fracStep,
	const Terrain *pTerrain) const
{
	// Generate normals & colors for vertices
	vector3d *vrts = borderVertexs;
	Color3ub *col = colors;
	vector3f *nrm = normals;
	double *hts = heights;

	const vector3d topEnd(0.0, yoffset + halfLen, 0.0);		// vertices at top edge of circle
	const vector3d btmEnd(0.0, yoffset - halfLen, 0.0);		// vertices at bottom edge of circle

	// step over the small square
	for ( int y = 0; y < edgeLen; y++ ) {
		const int by = (y + BORDER_SIZE) + yoff;
		for ( int x = 0; x < edgeLen; x++ ) {
			const int bx = (x + BORDER_SIZE) + xoff;

			// height
			const double height = borderHeights[bx + (by * borderedEdgeLen)];
			assert(hts != &heights[edgeLen * edgeLen]);
			*(hts++) = height;

			// normal
			const vector3d &x1 = vrts[(bx - 1) + (by * borderedEdgeLen)];
			const vector3d &x2 = vrts[(bx + 1) + (by * borderedEdgeLen)];
			const vector3d &y1 = vrts[bx + ((by - 1) * borderedEdgeLen)];
			const vector3d &y2 = vrts[bx + ((by + 1) * borderedEdgeLen)];
			const vector3d n = ((x2 - x1).Cross(y2 - y1)).Normalized();
			assert(nrm != &normals[edgeLen * edgeLen]);
			*(nrm++) = vector3f(n);

			// color
			const vector3d p = GetSurfacePointCyl(x * fracStep, y * fracStep, ang0, ang1, yoffset, halfLen);
			SetColour(col, pTerrain->GetColor(p, -height, -n));
			assert(col != &colors[edgeLen * edgeLen]);
			++col;
		}
	}
	assert(hts == &heights[edgeLen*edgeLen]);
	assert(nrm == &normals[edgeLen*edgeLen]);
	assert(col == &colors[edgeLen*edgeLen]);
}
