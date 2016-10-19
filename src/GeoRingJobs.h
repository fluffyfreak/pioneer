#ifndef _GEORINGJOBS_H
#define _GEORINGJOBS_H

#include "vector3.h"
#include "Camera.h"
#include "terrain/Terrain.h"

extern int GEOPATCH_EDGELEN;
#define ATMOSPHERE_RADIUS 1.015

namespace Graphics { 
	class Renderer; 
	class RenderState;
	class Material;
}
class SystemBody;
class GeoPlate;
class GeoPlateHull;
class GeoPlateWall;
class GeoPatchContext;

class SQuadPlateRequest;
class SQuadPlateResult;
class SSinglePlateResult;



class SBaseSplitRequest {
public:
	SBaseSplitRequest(const double ang0_, const double ang1_, const double yoffset_, const double halfLen_, const vector3d &vb0_, const vector3d &vb1_,
		const uint32_t depth_, const SystemPath &sysPath_, const GeoPlateID &patchID_, const int edgeLen_, const double fracStep_,
		Terrain *pTerrain_);

	inline int NUMVERTICES(const int el) const { return el*el; }

	const double ang0, ang1, yoffset, halfLen;
	const vector3d vbe0, vbe1;
	const uint32_t depth;
	const SystemPath sysPath;
	const GeoPlateID patchID;
	const int edgeLen;
	const double fracStep;
	RefCountedPtr<Terrain> pTerrain;

protected:
	// deliberately prevent copy constructor access
	SBaseSplitRequest(const SBaseSplitRequest &r) = delete;
};

class SQuadPlateRequest : public SBaseSplitRequest {
public:
	SQuadPlateRequest(const double ang0_, const double ang1_, const double yoffset_, const double halfLen_, const vector3d &vb0_, const vector3d &vb1_,
		const uint32_t depth_, const SystemPath &sysPath_, const GeoPlateID &patchID_, const int edgeLen_, const double fracStep_,
		Terrain *pTerrain_);

	// these are created with the request and are given to the resulting patches
	vector3f *normals[4];
	Color3ub *colors[4];
	double *heights[4];

	// these are created with the request but are destroyed when the request is finished
	std::unique_ptr<double[]> borderHeights;
	std::unique_ptr<vector3d[]> borderVertexs;

protected:
	// deliberately prevent copy constructor access
	SQuadPlateRequest(const SQuadPlateRequest &r) = delete;
};

class SSinglePlateRequest : public SBaseSplitRequest {
public:
	SSinglePlateRequest(const double ang0_, const double ang1_, const double yoffset_, const double halfLen_, const vector3d &vb0_, const vector3d &vb1_,
		const uint32_t depth_, const SystemPath &sysPath_, const GeoPlateID &patchID_, const int edgeLen_, const double fracStep_,
		Terrain *pTerrain_);

	// these are created with the request and are given to the resulting patches
	vector3f *normals;
	Color3ub *colors;
	double *heights;

	// these are created with the request but are destroyed when the request is finished
	std::unique_ptr<double> borderHeights;
	std::unique_ptr<vector3d> borderVertexs;

protected:
	// deliberately prevent copy constructor access
	SSinglePlateRequest(const SSinglePlateRequest &r) = delete;
};

class SBaseSplitResult {
public:
	struct SSplitResultData {
		SSplitResultData() : patchID(0) {}
		SSplitResultData(double *heights_, vector3f *n_, Color3ub *c_, const double ang0_, const double ang1_, const double yoffset_, const double halfLen_, const vector3d &vbe0_, const vector3d &vbe1_, const GeoPlateID &patchID_) :
			heights(heights_), normals(n_), colors(c_), ang0(ang0_), ang1(ang1_), yoffset(yoffset_), halfLen(halfLen_), vbe0(vbe0_), vbe1(vbe1_), patchID(patchID_)
		{}
		SSplitResultData(const SSplitResultData &r) : 
			normals(r.normals), colors(r.colors), ang0(r.ang0), ang1(r.ang1), yoffset(r.yoffset), halfLen(r.halfLen), vbe0(r.vbe0), vbe1(r.vbe1), patchID(r.patchID)
		{}

		double *heights;
		vector3f *normals;
		Color3ub *colors;
		double ang0, ang1, yoffset, halfLen;
		vector3d vbe0, vbe1;
		GeoPlateID patchID;
	};

	SBaseSplitResult(const uint64_t plate_, const int32_t depth_) : mPlate(plate_), mDepth(depth_) {}
	virtual ~SBaseSplitResult() {}

	inline uint64_t plate() const { return mPlate; }
	inline int32_t depth() const { return mDepth; }

	virtual void OnCancel() = 0;

protected:
	// deliberately prevent copy constructor access
	SBaseSplitResult(const SBaseSplitResult &r) : mPlate(0), mDepth(0) {}

	const uint64_t mPlate;
	const int32_t mDepth;
};

class SQuadPlateResult : public SBaseSplitResult {
	static const int NUM_RESULT_DATA = 4;
public:
	SQuadPlateResult(const uint64_t plate_, const int32_t depth_) : SBaseSplitResult(plate_, depth_)
	{
	}

	void addResult(const int kidIdx, double *h_, vector3f *n_, Color3ub *c_, const double ang0_, const double ang1_, const double yoffset_, const double halfLen_, const vector3d &vbe0, const vector3d &vbe1, const GeoPlateID &patchID_)
	{
		assert(kidIdx>=0 && kidIdx<NUM_RESULT_DATA);
		mData[kidIdx] = (SSplitResultData(h_, n_, c_, ang0_, ang1_, yoffset_, halfLen_, vbe0, vbe1, patchID_));
	}

	inline const SSplitResultData& data(const int32_t idx) const { return mData[idx]; }

	virtual void OnCancel()
	{
		for( int i=0; i<NUM_RESULT_DATA; ++i ) {
			if( mData[i].heights ) {delete [] mData[i].heights;		mData[i].heights = NULL;}
			if( mData[i].normals ) {delete [] mData[i].normals;		mData[i].normals = NULL;}
			if( mData[i].colors ) {delete [] mData[i].colors;		mData[i].colors = NULL;}
		}
	}

protected:
	// deliberately prevent copy constructor access
	SQuadPlateResult(const SQuadPlateResult &r) : SBaseSplitResult(r) {}

	SSplitResultData mData[NUM_RESULT_DATA];
};

class SSinglePlateResult : public SBaseSplitResult {
public:
	SSinglePlateResult(const uint64_t plate_, const int32_t depth_) : SBaseSplitResult(plate_, depth_)
	{
	}

	void addResult(double *h_, vector3f *n_, Color3ub *c_, const double ang0_, const double ang1_, const double yoffset_, const double halfLen_, const vector3d &vbe0, const vector3d &vbe1, const GeoPlateID &patchID_)
	{
		mData = (SSplitResultData(h_, n_, c_, ang0_, ang1_, yoffset_, halfLen_, vbe0, vbe1, patchID_));
	}

	inline const SSplitResultData& data() const { return mData; }

	virtual void OnCancel()
	{
		{
			if( mData.heights ) {delete [] mData.heights;	mData.heights = NULL;}
			if( mData.normals ) {delete [] mData.normals;	mData.normals = NULL;}
			if( mData.colors ) {delete [] mData.colors;		mData.colors = NULL;}
		}
	}

protected:
	// deliberately prevent copy constructor access
	SSinglePlateResult(const SSinglePlateResult &r) : SBaseSplitResult(r) {}

	SSplitResultData mData;
};

// ********************************************************************************
// Overloaded PureJob class to handle generating the mesh for each patch
// ********************************************************************************
class BasePlateJob : public Job
{
public:
	BasePlateJob() {}
	virtual void OnRun() override {}    // RUNS IN ANOTHER THREAD!! MUST BE THREAD SAFE!
	virtual void OnFinish() override {}
	virtual void OnCancel() override {}
};

// ********************************************************************************
// Overloaded PureJob class to handle generating the mesh for each patch
// ********************************************************************************
class SinglePlateJob : public BasePlateJob
{
public:
	SinglePlateJob(SSinglePlateRequest *data) : mData(data), mpResults(NULL) { /* empty */ }
	virtual ~SinglePlateJob();

	virtual void OnRun() override final;      // RUNS IN ANOTHER THREAD!! MUST BE THREAD SAFE!
	virtual void OnFinish() override final;   // runs in primary thread of the context

private:
	// Generates full-detail vertices, and also non-edge normals and colors 
	void GenerateMesh(double *heights, vector3f *normals, Color3ub *colors, double *borderHeights, vector3d *borderVertexs,
		const double ang0, const double ang1, const double yoffset, const double halfLen,
		const int edgeLen, const double fracStep, const Terrain *pTerrain) const;

	std::unique_ptr<SSinglePlateRequest> mData;
	SSinglePlateResult *mpResults;
};

// ********************************************************************************
// Overloaded PureJob class to handle generating the mesh for each patch
// ********************************************************************************
class QuadPlateJob : public BasePlateJob
{
public:
	QuadPlateJob(SQuadPlateRequest *data) : mData(data), mpResults(NULL) { /* empty */ }
	virtual ~QuadPlateJob();

	virtual void OnRun() override final;      // RUNS IN ANOTHER THREAD!! MUST BE THREAD SAFE!
	virtual void OnFinish() override final;   // runs in primary thread of the context

private:
	// Generates full-detail vertices, and also non-edge normals and colors 
	void GenerateBorderedData(double *borderHeights, vector3d *borderVertexs,
		const double ang0, const double ang1, const double yoffset, const double halfLen,
		const int edgeLen, const double fracStep, const Terrain *pTerrain) const;
	
	void GenerateSubPatchData(double *heights, vector3f *normals, Color3ub *colors, double *borderHeights, vector3d *borderVertexs,
		const double ang0, const double ang1, const double yoffset, const double halfLen,
		const int edgeLen, const int xoff, const int yoff, const int borderedEdgeLen, const double fracStep, const Terrain *pTerrain) const;

	std::unique_ptr<SQuadPlateRequest> mData;
	SQuadPlateResult *mpResults;
};

#endif /* _GEORINGJOBS_H */
