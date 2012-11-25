// Copyright © 2008-2012 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef __GEOMESH_H__
#define __GEOMESH_H__

#include "GeoPatchID.h"

#include "vector3.h"
#include "mtrand.h"
#include "galaxy/StarSystem.h"
#include "graphics/Frustum.h"
#include "graphics/Material.h"
#include "terrain/Terrain.h"

#include <deque>

static const uint32_t NUM_SIDES = 6;

namespace Graphics { class Renderer; }

class SystemBody;
// fwd decl'
class GeoPatch;
class GeoPatchContext;

struct SSplitRequestDescription {
	SSplitRequestDescription(const vector3f &v0_,
							const vector3f &v1_,
							const vector3f &v2_,
							const vector3f &v3_,
							const vector3f &cn,
							const uint32_t depth_,
							const GeoPatchID &patchID_)
							: v0(v0_), v1(v1_), v2(v2_), v3(v3_), centroid(cn), depth(depth_), patchID(patchID_)
	{
	}

	const vector3f v0;
	const vector3f v1;
	const vector3f v2;
	const vector3f v3;
	const vector3f centroid;
	const uint32_t depth;
	const GeoPatchID patchID;
};

struct SSplitResult {
	struct SSplitResultData {
		SSplitResultData(const GLuint texID_, const vector3f &v0_, const vector3f &v1_, const vector3f &v2_, const vector3f &v3_, const GeoPatchID &patchID_) :
			texID(texID_), v0(v0_), v1(v1_), v2(v2_), v3(v3_), patchID(patchID_)
		{
		}
		const GLuint texID;
		const vector3f v0;
		const vector3f v1;
		const vector3f v2;
		const vector3f v3;
		const GeoPatchID patchID;
	};

	SSplitResult(const int32_t face_, const uint32_t depth_) : face(face_), depth(depth_)
	{
	}

	void addResult(const GLuint tex, const vector3f &v0_, const vector3f &v1_, const vector3f &v2_, const vector3f &v3_, const GeoPatchID &patchID_)
	{
		data.push_back(SSplitResultData(tex, v0_, v1_, v2_, v3_, patchID_));
		assert(data.size()<=4);
	}

	const int32_t face;
	const uint32_t depth;
	std::deque<SSplitResultData> data;
};

class GeoMesh
{
private:
	void BuildFirstPatches();
	void Render(Graphics::Renderer *renderer, const vector3f& campos, const float radius, const float scale);

	static RefCountedPtr<GeoPatchContext> sPatchContext;
	static const uint32_t NUM_PATCHES = 6;
	GeoPatch*			mGeoPatches[NUM_PATCHES];
	const SystemBody*	mSystemBody;

	/* all variables for GetHeight(), GetColor() */
	const Terrain *mTerrain;

	static const uint32_t MAX_SPLIT_REQUESTS = 128;
	std::deque<SSplitRequestDescription*> mSplitRequestDescriptions;
	std::deque<SSplitResult*> mSplitResult;

	void SetUpMaterials();
	ScopedPtr<Graphics::Material> m_surfaceMaterial;
	ScopedPtr<Graphics::Material> m_atmosphereMaterial;
	//special parameters for shaders
	SystemBody::AtmosphereParameters m_atmosphereParameters;

public:
	GeoMesh(const SystemBody *body);
	~GeoMesh();

	static void Init();
	static void Uninit();
	static void OnChangeDetailLevel();

	void Update(const vector3f &campos);
	void Render(const matrix4x4f &ViewMatrix, const matrix4x4f &ModelMatrix, const matrix4x4f &MVP, Graphics::Renderer *renderer, const vector3f& campos, const float radius, const float scale);

	bool AddSplitRequest(SSplitRequestDescription *desc);
	void ProcessSplitRequests();
	void ProcessSplitResults();
};

#endif //__GEOMESH_H__