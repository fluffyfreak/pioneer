// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _GEOPATCHJOBSGPU_H
#define _GEOPATCHJOBSGPU_H

#include <SDL_stdinc.h>

#include "vector3.h"
#include "Random.h"
#include "galaxy/StarSystem.h"
#include "terrain/Terrain.h"
#include "GeoPatchID.h"
#include "JobQueue.h"
#include "JobQueueGPU.h"

class GeoSphere;

class SBaseRequestGPU {
public:
	SBaseRequestGPU(const vector3d &v0_, const vector3d &v1_, const vector3d &v2_, const vector3d &v3_, const vector3d &cn,
		const uint32_t depth_, const SystemPath &sysPath_, const GeoPatchID &patchID_, const int edgeLen_, const double fracStep_,
		Terrain *pTerrain_)
		: v0(v0_), v1(v1_), v2(v2_), v3(v3_), centroid(cn), depth(depth_), 
		sysPath(sysPath_), patchID(patchID_), edgeLen(edgeLen_), fracStep(fracStep_), 
		pTerrain(pTerrain_)
	{
	}

	inline int NUMVERTICES() const { return edgeLen*edgeLen; }
	inline int NUMVERTICES(const int el) const { return el*el; }

	const vector3d v0, v1, v2, v3;
	const vector3d centroid;
	const uint32_t depth;
	const SystemPath sysPath;
	const GeoPatchID patchID;
	const int edgeLen;
	const double fracStep;
	RefCountedPtr<Terrain> pTerrain;

protected:
	// deliberately prevent copy constructor access
	SBaseRequestGPU(const SBaseRequestGPU &r) : v0(0.0), v1(0.0), v2(0.0), v3(0.0), centroid(0.0), depth(0), 
		patchID(0), edgeLen(0), fracStep(0.0), pTerrain(NULL) { assert(false); }
};

class SQuadRequestGPU : public SBaseRequestGPU {
public:
	SQuadRequestGPU(const vector3d &v0_, const vector3d &v1_, const vector3d &v2_, const vector3d &v3_, const vector3d &cn,
		const uint32_t depth_, const SystemPath &sysPath_, const GeoPatchID &patchID_, const int edgeLen_, const double fracStep_,
		Terrain *pTerrain_)
		: SBaseRequestGPU(v0_, v1_, v2_, v3_, cn, depth_, sysPath_, patchID_, edgeLen_, fracStep_, pTerrain_)
	{
		const int numVerts = NUMVERTICES(edgeLen_);
		const int numBorderedVerts = NUMVERTICES(edgeLen_+2);
		for( int i=0 ; i<4 ; ++i )
		{
			heights[i] = new double[numVerts];
			normals[i] = new vector3f[numVerts];
			colors[i] = new Color3ub[numVerts];

			borderHeights[i].reset(new double[numBorderedVerts]);
			borderVertexs[i].reset(new vector3d[numBorderedVerts]);
		}
	}

	// these are created with the request and are given to the resulting patches
	vector3f *normals[4];
	Color3ub *colors[4];
	double *heights[4];

	// these are created with the request but are destroyed when the request is finished
	std::unique_ptr<double[]> borderHeights[4];
	std::unique_ptr<vector3d[]> borderVertexs[4];

protected:
	// deliberately prevent copy constructor access
	SQuadRequestGPU(const SQuadRequestGPU &r) : SBaseRequestGPU(r)	{ assert(false); }
};

class SSingleRequestGPU : public SBaseRequestGPU {
public:
	SSingleRequestGPU(const vector3d &v0_, const vector3d &v1_, const vector3d &v2_, const vector3d &v3_, const vector3d &cn,
		const uint32_t depth_, const SystemPath &sysPath_, const GeoPatchID &patchID_, const int edgeLen_, const double fracStep_,
		Terrain *pTerrain_)
		: SBaseRequestGPU(v0_, v1_, v2_, v3_, cn, depth_, sysPath_, patchID_, edgeLen_, fracStep_, pTerrain_)
	{
		const int numVerts = NUMVERTICES(edgeLen_);
		heights = new double[numVerts];
		normals = new vector3f[numVerts];
		colors = new Color3ub[numVerts];
		
		const int numBorderedVerts = NUMVERTICES(edgeLen_+2);
		borderHeights.reset(new double[numBorderedVerts]);
		borderVertexs.reset(new vector3d[numBorderedVerts]);
	}

	// these are created with the request and are given to the resulting patches
	vector3f *normals;
	Color3ub *colors;
	double *heights;

	// these are created with the request but are destroyed when the request is finished
	std::unique_ptr<double> borderHeights;
	std::unique_ptr<vector3d> borderVertexs;

protected:
	// deliberately prevent copy constructor access
	SSingleRequestGPU(const SSingleRequestGPU &r) : SBaseRequestGPU(r)	{ assert(false); }
};

class SBaseResultGPU {
public:
	struct SResultDataGPU {
		SResultDataGPU() : patchID(0) {}
		SResultDataGPU(double *heights_, vector3f *n_, Color3ub *c_, const vector3d &v0_, const vector3d &v1_, const vector3d &v2_, const vector3d &v3_, const GeoPatchID &patchID_) :
			heights(heights_), normals(n_), colors(c_), v0(v0_), v1(v1_), v2(v2_), v3(v3_), patchID(patchID_)
		{}
		SResultDataGPU(const SResultDataGPU &r) : 
			normals(r.normals), colors(r.colors), v0(r.v0), v1(r.v1), v2(r.v2), v3(r.v3), patchID(r.patchID)
		{}

		double *heights;
		vector3f *normals;
		Color3ub *colors;
		vector3d v0, v1, v2, v3;
		GeoPatchID patchID;
	};

	SBaseResultGPU(const int32_t face_, const int32_t depth_) : mFace(face_), mDepth(depth_) {}
	virtual ~SBaseResultGPU() {}

	inline int32_t face() const { return mFace; }
	inline int32_t depth() const { return mDepth; }

	virtual void OnCancel() = 0;

protected:
	// deliberately prevent copy constructor access
	SBaseResultGPU(const SBaseResultGPU &r) : mFace(0), mDepth(0) {}

	const int32_t mFace;
	const int32_t mDepth;
};

class SQuadResultGPU : public SBaseResultGPU {
	static const int NUM_RESULT_DATA = 4;
public:
	SQuadResultGPU(const int32_t face_, const int32_t depth_) : SBaseResultGPU(face_, depth_)
	{
	}

	void addResult(const int kidIdx, double *h_, vector3f *n_, Color3ub *c_, const vector3d &v0_, const vector3d &v1_, const vector3d &v2_, const vector3d &v3_, const GeoPatchID &patchID_)
	{
		assert(kidIdx>=0 && kidIdx<NUM_RESULT_DATA);
		mData[kidIdx] = (SResultDataGPU(h_, n_, c_, v0_, v1_, v2_, v3_, patchID_));
	}

	inline const SResultDataGPU& data(const int32_t idx) const { return mData[idx]; }

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
	SQuadResultGPU(const SQuadResultGPU &r) : SBaseResultGPU(r) {}

	SResultDataGPU mData[NUM_RESULT_DATA];
};

class SSingleResultGPU : public SBaseResultGPU {
public:
	SSingleResultGPU(const int32_t face_, const int32_t depth_) : SBaseResultGPU(face_, depth_)
	{
	}

	void addResult(double *h_, vector3f *n_, Color3ub *c_, const vector3d &v0_, const vector3d &v1_, const vector3d &v2_, const vector3d &v3_, const GeoPatchID &patchID_)
	{
		mData = (SResultDataGPU(h_, n_, c_, v0_, v1_, v2_, v3_, patchID_));
	}

	inline const SResultDataGPU& data() const { return mData; }

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
	SSingleResultGPU(const SSingleResultGPU &r) : SBaseResultGPU(r) {}

	SResultDataGPU mData;
};

class GeoPatch;

// ********************************************************************************
// Overloaded PureJob class to handle generating the mesh for each patch
// ********************************************************************************
class BasePatchJobGPU : public JobGPU
{
public:
	BasePatchJobGPU() {}
	virtual void OnRun() = 0;    // RUNS IN ANOTHER THREAD!! MUST BE THREAD SAFE!
	virtual void OnFinish() = 0;

protected:
	// in patch surface coords, [0,1]
	inline vector3d GetSpherePoint(const vector3d &v0, const vector3d &v1, const vector3d &v2, const vector3d &v3, const double x, const double y) const {
		return (v0 + x*(1.0-y)*(v1-v0) + x*y*(v2-v0) + (1.0-x)*y*(v3-v0)).Normalized();
	}

	// Generates full-detail vertices, and also non-edge normals and colors 
	void GenerateMesh(double *heights, vector3f *normals, Color3ub *colors, double *borderHeights, vector3d *borderVertexs,
		const vector3d &v0, const vector3d &v1, const vector3d &v2, const vector3d &v3,
		const int edgeLen, const double fracStep, const Terrain *pTerrain) const;
};

// ********************************************************************************
// Overloaded PureJob class to handle generating the mesh for each patch
// ********************************************************************************
class SinglePatchJobGPU : public BasePatchJobGPU
{
public:
	SinglePatchJobGPU(SSingleRequestGPU *data) : mData(data), mpResults(NULL) { /* empty */ }
	~SinglePatchJobGPU();

	virtual void OnRun();      // RUNS IN ANOTHER THREAD!! MUST BE THREAD SAFE!
	virtual void OnFinish();   // runs in primary thread of the context

private:
	std::unique_ptr<SSingleRequestGPU> mData;
	SSingleResultGPU *mpResults;
};

// ********************************************************************************
// Overloaded PureJob class to handle generating the mesh for each patch
// ********************************************************************************
class QuadPatchJobGPU : public BasePatchJobGPU
{
public:
	QuadPatchJobGPU(SQuadRequestGPU *data) : mData(data), mpResults(NULL) { /* empty */ }
	~QuadPatchJobGPU();

	virtual void OnRun();      // RUNS IN ANOTHER THREAD!! MUST BE THREAD SAFE!
	virtual void OnFinish();   // runs in primary thread of the context

private:
	std::unique_ptr<SQuadRequestGPU> mData;
	SQuadResultGPU *mpResults;
};



// ********************************************************************************
// a quad with reversed winding
class GenFaceQuad {
public:
	GenFaceQuad(Graphics::Renderer *r, const vector2f &size, Graphics::RenderState *state, const Uint32 GGQuality);
	virtual void Draw(Graphics::Renderer *r);

	void SetMaterial(Graphics::Material *mat) { assert(mat); m_material.reset(mat); }
	Graphics::Material* GetMaterial() const { return m_material.get(); }
private:
	std::unique_ptr<Graphics::Material> m_material;
	std::unique_ptr<Graphics::VertexBuffer> m_vertexBuffer;
	Graphics::RenderState *m_renderState;
};

// ********************************************************************************
class SGPUGenRequest {
public:
	SGPUGenRequest(const vector3d *v_, const SystemPath &sysPath_, const Sint32 face_, const Sint32 uvDIMs_, vector3f frequency_, const float planetRadius_, GenFaceQuad* pQuad_, Graphics::Texture *pTex_);

	inline Sint32 Face() const { return face; }
	inline Sint32 UVDims() const { return uvDIMs; }
	inline Graphics::Texture* Texture() const { return m_texture.Get(); }
	inline GenFaceQuad* Quad() const { return pQuad; }
	inline const SystemPath& SysPath() const { return sysPath; }

	void SetupMaterialParams()
	{
		PROFILE_SCOPED()
		m_specialParams.v = corners;
		m_specialParams.fracStep = 1.0f / float(uvDIMs);
		m_specialParams.planetRadius = planetRadius;
		m_specialParams.time = 0.0f;
		m_specialParams.frequency = frequency;
		pQuad->GetMaterial()->specialParameter0 = &m_specialParams;
	}

#ifdef DUMP_PARAMS
	void DumpParams(const int face)
	{
		//
		Output("--face-%d------------\n", face);
		const vector3d *v = corners;
		for (int i = 0; i < 4; i++)
			Output("v(%.4f,%.4f,%.4f)\n", v[i].x, v[i].y, v[i].z);

		Output("fracStep = %.4f\n", 1.0f / float(uvDIMs));
		Output("time = 0.0f\n");

		for (Uint32 i = 0; i<10; i++) {
			const fracdef_t &fd = pTerrain->GetFracDef(i);
			Output("frequency[%u] = %.4f\n", i, fd.frequency);
		}
		// XXX omg hacking galore
	}
#endif

protected:
	// deliberately prevent copy constructor access
	SGPUGenRequest(const SGPUGenRequest &r);

	inline Sint32 NumTexels() const { return uvDIMs*uvDIMs; }

	// this is created with the request and are given to the resulting patches
	RefCountedPtr<Graphics::Texture> m_texture;

	const vector3d *corners;
	const SystemPath sysPath;
	const Sint32 face;
	const Sint32 uvDIMs;
	const vector3f frequency;
	const float planetRadius;
	GenFaceQuad* pQuad;
	Graphics::GenGasGiantColourMaterialParameters m_specialParams;
};

// ********************************************************************************
class SGPUGenResult {
public:
	struct SGPUGenData {
		SGPUGenData() {}
		SGPUGenData(Graphics::Texture *t_, Sint32 uvDims_) : texture(t_), uvDims(uvDims_) {}
		SGPUGenData(const SGPUGenData &r) : texture(r.texture), uvDims(r.uvDims) {}
		RefCountedPtr<Graphics::Texture> texture;
		Sint32 uvDims;
	};

	SGPUGenResult(const int32_t face_) : mFace(face_) {}

	void addResult(Graphics::Texture *t_, Sint32 uvDims_);

	inline const SGPUGenData& data() const { return mData; }
	inline int32_t face() const { return mFace; }

	void OnCancel();

protected:
	// deliberately prevent copy constructor access
	SGPUGenResult(const SGPUGenResult &r);

	const int32_t mFace;
	SGPUGenData mData;
};

// ********************************************************************************
// Overloaded JobGPU class to handle generating the mesh for each patch
// ********************************************************************************
class SingleGPUGenJob : public JobGPU
{
public:
	SingleGPUGenJob(SGPUGenRequest *data);
	virtual ~SingleGPUGenJob();

	virtual void OnRun();
	virtual void OnFinish();
	virtual void OnCancel() {}

private:
	SingleGPUGenJob() {}
	// deliberately prevent copy constructor access
	SingleGPUGenJob(const SingleGPUGenJob &r);

	std::unique_ptr<SGPUGenRequest> mData;
	SGPUGenResult *mpResults;
};

#endif /* _GEOPATCHJOBSGPU_H */
