// Copyright © 2008-2016 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "libs.h"
#include "GeoPatchJobs.h"
#include "GeoSphere.h"
#include "GeoPatch.h"
#include "perlin.h"
#include "Pi.h"
#include "RefCounted.h"

inline void setColour(Color3ub &r, const vector3d &v) { 
	r.r=static_cast<unsigned char>(Clamp(v.x*255.0, 0.0, 255.0)); 
	r.g=static_cast<unsigned char>(Clamp(v.y*255.0, 0.0, 255.0)); 
	r.b=static_cast<unsigned char>(Clamp(v.z*255.0, 0.0, 255.0));
}

#define NUM_SAMPLES 9 //17 //33
static vector2d spiral[NUM_SAMPLES];
void GenerateSpiral() {
	spiral[0] = vector2d(0.0);
	for ( int i = 1; i < NUM_SAMPLES; i++ ) 
	{
		double s = double(i) / (double(NUM_SAMPLES) - 1.0);
		double a = sqrt(s*512.0);
		double b = sqrt(s);
		spiral[i].x = sin(a)*b;
		spiral[i].y = cos(a)*b;
	}
}
	

// ********************************************************************************
// Overloaded PureJob class to handle generating the mesh for each patch
// ********************************************************************************

#define USE_OCCLUSION
double BasePatchJob::Occlusion(const vector3d &normal, const int x, const int y, const double fracStep, 
	const vector3d &v0, const vector3d &v1, const vector3d &v2, const vector3d &v3,
	const Terrain *pTerrain) const
{
#ifdef USE_OCCLUSION
	const vector3d pos = GetSpherePoint(v0, v1, v2, v3, x*fracStep, y*fracStep);
	double occlusion = 0.0;
	double raylength = 2.0;

	for(int i=1; i<NUM_SAMPLES; i++)
	{
		double srx = spiral[i].x * raylength;
		double sry = spiral[i].y * raylength;
		vector3d sample_pos = GetSpherePoint(v0, v1, v2, v3, (x+srx)*fracStep, (y+sry)*fracStep);
		double h = pTerrain->GetHeight(sample_pos);
		sample_pos *= h;
		vector3d sample_dir = (sample_pos - pos).Normalized();
		double lambert = Clamp(normal.Dot(sample_dir), 0.0, 1.0);
		double dist_factor = 0.23/sqrt((sample_pos - pos).Length());
		occlusion += dist_factor*lambert;
	}
	//return Clamp((occlusion/(double(NUM_SAMPLES) - 1.0)), 0.0, 1.0);
	return 1.0 - (occlusion/(double(NUM_SAMPLES) - 1.0));
#else
	return 1.0;
#endif
}

// Generates full-detail vertices, and also non-edge normals and colors 
void BasePatchJob::GenerateMesh(double *heights, vector3f *normals, Color3ub *colors, 
								double *borderHeights, vector3d *borderVertexs,
								const vector3d &v0,
								const vector3d &v1,
								const vector3d &v2,
								const vector3d &v3,
								const int edgeLen,
								const double fracStep,
								const Terrain *pTerrain) const
{
	const int borderedEdgeLen = edgeLen+2;
	const int numBorderedVerts = borderedEdgeLen*borderedEdgeLen;

	// generate heights plus a 1 unit border
	double *bhts = borderHeights;
	vector3d *vrts = borderVertexs;
	for (int y=-1; y<borderedEdgeLen-1; y++) {
		const double yfrac = double(y) * fracStep;
		for (int x=-1; x<borderedEdgeLen-1; x++) {
			const double xfrac = double(x) * fracStep;
			const vector3d p = GetSpherePoint(v0, v1, v2, v3, xfrac, yfrac);
			const double height = pTerrain->GetHeight(p);
			assert(height >= 0.0f && height <= 1.0f);
			*(bhts++) = height;
			*(vrts++) = p * (height + 1.0);
		}
	}
	assert(bhts==&borderHeights[numBorderedVerts]);

	// Generate normals & colors for non-edge vertices since they never change
	Color3ub *col = colors;
	vector3f *nrm = normals;
	double *hts = heights;
	vrts = borderVertexs;
	for (int y=1; y<borderedEdgeLen-1; y++) {
		for (int x=1; x<borderedEdgeLen-1; x++) {
			// height
			const double height = borderHeights[x + y*borderedEdgeLen];
			assert(hts!=&heights[edgeLen*edgeLen]);
			*(hts++) = height;

			// normal
			const vector3d &x1 = vrts[x-1 + y*borderedEdgeLen];
			const vector3d &x2 = vrts[x+1 + y*borderedEdgeLen];
			const vector3d &y1 = vrts[x + (y-1)*borderedEdgeLen];
			const vector3d &y2 = vrts[x + (y+1)*borderedEdgeLen];
			const vector3d n = ((x2-x1).Cross(y2-y1)).Normalized();
			assert(nrm!=&normals[edgeLen*edgeLen]);
			*(nrm++) = vector3f(n);

			// color
			const vector3d p = GetSpherePoint(v0, v1, v2, v3, (x-1)*fracStep, (y-1)*fracStep);
			const double occ = Occlusion(n, (x-1), (y-1), fracStep, v0, v1, v2, v3,pTerrain);
			setColour(*col, pTerrain->GetColor(p, height, n) * occ);
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
void SinglePatchJob::OnFinish()  // runs in primary thread of the context
{
	GeoSphere::OnAddSingleSplitResult( mData->sysPath, mpResults );
	mpResults = nullptr;
	BasePatchJob::OnFinish();
}

void SinglePatchJob::OnRun()    // RUNS IN ANOTHER THREAD!! MUST BE THREAD SAFE!
{
	BasePatchJob::OnRun();

	const SSingleSplitRequest &srd = *mData;

	// fill out the data
	GenerateMesh(srd.heights, srd.normals, srd.colors, srd.borderHeights.get(), srd.borderVertexs.get(),
		srd.v0, srd.v1, srd.v2, srd.v3, 
		srd.edgeLen, srd.fracStep, srd.pTerrain.Get());
	// add this patches data
	SSingleSplitResult *sr = new SSingleSplitResult(srd.patchID.GetPatchFaceIdx(), srd.depth);
	sr->addResult(srd.heights, srd.normals, srd.colors, 
		srd.v0, srd.v1, srd.v2, srd.v3, 
		srd.patchID.NextPatchID(srd.depth+1, 0));
	// store the result
	mpResults = sr;
}

SinglePatchJob::~SinglePatchJob()
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
void QuadPatchJob::OnFinish()  // runs in primary thread of the context
{
	GeoSphere::OnAddQuadSplitResult( mData->sysPath, mpResults );
	mpResults = nullptr;
	BasePatchJob::OnFinish();
}

void QuadPatchJob::OnRun()    // RUNS IN ANOTHER THREAD!! MUST BE THREAD SAFE!
{
	BasePatchJob::OnRun();

	const SQuadSplitRequest &srd = *mData;
	const vector3d v01	= (srd.v0+srd.v1).Normalized();
	const vector3d v12	= (srd.v1+srd.v2).Normalized();
	const vector3d v23	= (srd.v2+srd.v3).Normalized();
	const vector3d v30	= (srd.v3+srd.v0).Normalized();
	const vector3d cn	= (srd.centroid).Normalized();

	const vector3d vecs[4][4] = {
		{srd.v0,	v01,		cn,			v30},
		{v01,		srd.v1,		v12,		cn},
		{cn,		v12,		srd.v2,		v23},
		{v30,		cn,			v23,		srd.v3}
	};

	SQuadSplitResult *sr = new SQuadSplitResult(srd.patchID.GetPatchFaceIdx(), srd.depth);
	for (int i=0; i<4; i++)
	{
		// fill out the data
		GenerateMesh(srd.heights[i], srd.normals[i], srd.colors[i], srd.borderHeights[i].get(), srd.borderVertexs[i].get(),
			vecs[i][0], vecs[i][1], vecs[i][2], vecs[i][3], 
			srd.edgeLen, srd.fracStep, srd.pTerrain.Get());
		// add this patches data
		sr->addResult(i, srd.heights[i], srd.normals[i], srd.colors[i], 
			vecs[i][0], vecs[i][1], vecs[i][2], vecs[i][3], 
			srd.patchID.NextPatchID(srd.depth+1, i));
	}
	mpResults = sr;
}

QuadPatchJob::~QuadPatchJob()
{
	if(mpResults) {
		mpResults->OnCancel();
		delete mpResults;
		mpResults = NULL;
	}
}
