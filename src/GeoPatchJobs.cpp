// Copyright © 2008-2017 Pioneer Developers. See AUTHORS.txt for details
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

// in patch surface coords, [0,1]
inline vector3d GetSpherePoint(const vector3d &v0, const vector3d &v1, const vector3d &v2, const vector3d &v3, const double x, const double y) {
	return (v0 + x*(1.0-y)*(v1-v0) + x*y*(v2-v0) + (1.0-x)*y*(v3-v0)).Normalized();
}

// ********************************************************************************
// Overloaded PureJob class to handle generating the mesh for each patch
// ********************************************************************************

// Generates full-detail vertices, and also non-edge normals and colors
void SSingleSplitRequest::GenerateMesh() const
{
	const int borderedEdgeLen = edgeLen+(BORDER_SIZE*2);
	const int numBorderedVerts = borderedEdgeLen*borderedEdgeLen;

	// generate heights plus a 1 unit border
	double *bhts = borderHeights.get();
	vector3d *vrts = borderVertexs.get();
	for (int y=-BORDER_SIZE; y<borderedEdgeLen-BORDER_SIZE; y++) {
		const double yfrac = double(y) * fracStep;
		for (int x=-BORDER_SIZE; x<borderedEdgeLen-BORDER_SIZE; x++) {
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
	vrts = borderVertexs.get();
	for (int y=BORDER_SIZE; y<borderedEdgeLen-BORDER_SIZE; y++) {
		for (int x=BORDER_SIZE; x<borderedEdgeLen-BORDER_SIZE; x++) {
			// height
			const double height = borderHeights[x + y*borderedEdgeLen];
			*(hts++) = height;

			// normal
			const vector3d &x1 = vrts[(x-1) + y*borderedEdgeLen];
			const vector3d &x2 = vrts[(x+1) + y*borderedEdgeLen];
			const vector3d &y1 = vrts[x + (y-1)*borderedEdgeLen];
			const vector3d &y2 = vrts[x + (y+1)*borderedEdgeLen];
			const vector3d n = ((x2-x1).Cross(y2-y1)).Normalized();
			*(nrm++) = vector3f(n);

			// color
			const vector3d p = GetSpherePoint(v0, v1, v2, v3, (x-BORDER_SIZE)*fracStep, (y-BORDER_SIZE)*fracStep);
			setColour(*col, pTerrain->GetColor(p, height, n));
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
	timer.Start();

	BasePatchJob::OnRun();

	const SSingleSplitRequest &srd = *mData;

	// fill out the data
	mData->GenerateMesh();

	// add this patches data
	SSingleSplitResult *sr = new SSingleSplitResult(srd.patchID.GetPatchFaceIdx(), srd.depth);
	sr->addResult(srd.heights, srd.normals, srd.colors,
		srd.v0, srd.v1, srd.v2, srd.v3,
		srd.patchID.NextPatchID(srd.depth+1, 0));
	// store the result
	mpResults = sr;

	timer.Stop();
}


double SinglePatchJob::totalTime = 0.0f;
Uint64 SinglePatchJob::totalSinglePatches = 0;


SinglePatchJob::~SinglePatchJob()
{
	const double timing = timer.millicycles();
	totalTime += timing;
	++totalSinglePatches;
	Output("SinglePatchJob took: %lf millicycles, avg %lf ms\n", timing, totalTime / totalSinglePatches);

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
	timer.Start();

	BasePatchJob::OnRun();

	const SQuadSplitRequest &srd = *mData;

	mData->GenerateBorderedData();

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

	const int borderedEdgeLen = (srd.edgeLen*2)+(BORDER_SIZE*2)-1;
	const int offxy[4][2] = {
		{0,0},
		{srd.edgeLen-1,0},
		{srd.edgeLen-1,srd.edgeLen-1},
		{0,srd.edgeLen-1}
	};

	SQuadSplitResult *sr = new SQuadSplitResult(srd.patchID.GetPatchFaceIdx(), srd.depth);
	for (int i=0; i<4; i++)
	{
		// fill out the data
		mData->GenerateSubPatchData(i,
			vecs[i][0], vecs[i][1], vecs[i][2], vecs[i][3],
			srd.edgeLen, offxy[i][0], offxy[i][1],
			borderedEdgeLen);

		// add this patches data
		sr->addResult(i, srd.heights[i], srd.normals[i], srd.colors[i],
			vecs[i][0], vecs[i][1], vecs[i][2], vecs[i][3],
			srd.patchID.NextPatchID(srd.depth+1, i));
	}
	mpResults = sr;

	timer.Stop();
}


double QuadPatchJob::totalTime = 0.0f;
Uint64 QuadPatchJob::totalQuadPatches = 0;


QuadPatchJob::~QuadPatchJob()
{
	const double timing = timer.millicycles();
	totalTime += timing;
	++totalQuadPatches;
	Output("QuadPatchJob took: %lf millicycles, avg %lf ms\n", timing, totalTime / totalQuadPatches);

	if(mpResults) {
		mpResults->OnCancel();
		delete mpResults;
		mpResults = NULL;
	}
}

// Generates full-detail vertices, and also non-edge normals and colors
void SQuadSplitRequest::GenerateBorderedData() const
{
	const int borderedEdgeLen = (edgeLen * 2) + (BORDER_SIZE * 2) - 1;
	const int numBorderedVerts = borderedEdgeLen*borderedEdgeLen;

	// generate heights plus a N=BORDER_SIZE unit border
	double *bhts = borderHeights.get();
	vector3d *vrts = borderVertexs.get();
	for ( int y = -BORDER_SIZE; y < (borderedEdgeLen - BORDER_SIZE); y++ ) {
		const double yfrac = double(y) * (fracStep*0.5);
		for ( int x = -BORDER_SIZE; x < (borderedEdgeLen - BORDER_SIZE); x++ ) {
			const double xfrac = double(x) * (fracStep*0.5);
			const vector3d p = GetSpherePoint(v0, v1, v2, v3, xfrac, yfrac);
			const double height = pTerrain->GetHeight(p);
			assert(height >= 0.0f && height <= 1.0f);
			*(bhts++) = height;
			*(vrts++) = p * (height + 1.0);
		}
	}
	assert(bhts == &borderHeights[numBorderedVerts]);
}

void SQuadSplitRequest::GenerateSubPatchData(const int quadrantIndex,
	const vector3d &v0,
	const vector3d &v1,
	const vector3d &v2,
	const vector3d &v3,
	const int edgeLen,
	const int xoff,
	const int yoff,
	const int borderedEdgeLen) const
{
	// Generate normals & colors for vertices
	vector3d *vrts = borderVertexs.get();
	Color3ub *col = colors[quadrantIndex];
	vector3f *nrm = normals[quadrantIndex];
	double *hts = heights[quadrantIndex];

	// step over the small inner square
	// .........
	// .#######.
	// .#######.
	// .#######.
	// .#######.
	// .#######.
	// .#######.
	// .#######.
	// .#######.
	// .........
	for ( int y = 0; y < edgeLen; y++ ) {
		const int by = (y + BORDER_SIZE) + yoff;
		for ( int x = 0; x < edgeLen; x++ ) {
			const int bx = (x + BORDER_SIZE) + xoff;

			// height
			const double height = borderHeights[bx + (by * borderedEdgeLen)];
			*(hts++) = height;

			// normal
			const vector3d &x1 = vrts[(bx - 1) + (by * borderedEdgeLen)];
			const vector3d &x2 = vrts[(bx + 1) + (by * borderedEdgeLen)];
			const vector3d &y1 = vrts[bx + ((by - 1) * borderedEdgeLen)];
			const vector3d &y2 = vrts[bx + ((by + 1) * borderedEdgeLen)];
			const vector3d n = ((x2 - x1).Cross(y2 - y1)).Normalized();
			*(nrm++) = vector3f(n);

			// color
			const vector3d p = GetSpherePoint(v0, v1, v2, v3, x * fracStep, y * fracStep);
			setColour(*col, pTerrain->GetColor(p, height, n));
			++col;
		}
	}
	assert(hts == &heights[quadrantIndex][edgeLen*edgeLen]);
	assert(nrm == &normals[quadrantIndex][edgeLen*edgeLen]);
	assert(col == &colors[quadrantIndex][edgeLen*edgeLen]);

	/*const int colLen = edgeLen - 1;
	for ( int y = 0; y < colLen; y++ ) {
		for ( int x = 0; x < colLen; x++ ) {

			// height
			const double height = hts[0]; // sample from hts;

			// normal
			const vector3d n = vector3d(nrm[0]); // sample from nrm

			// color
			const vector3d p = GetSpherePoint(v0, v1, v2, v3, x * fracStep, y * fracStep);
			setColour(*col, pTerrain->GetColor(p, height, n));
			assert(col != &colors[colLen * colLen]);
			++col;
		}
	}*/
}
