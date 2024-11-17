// Copyright © 2008-2024 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _GEOPATCH_H
#define _GEOPATCH_H

#include <SDL_stdinc.h>

#include "Color.h"
#include "GeoPatchID.h"
#include "JobQueue.h"
#include "RefCounted.h"
#include "matrix4x4.h"
#include "vector3.h"
#include "graphics/Frustum.h"
#include <deque>
#include <memory>

//#define DEBUG_BOUNDING_SPHERES

#ifdef DEBUG_BOUNDING_SPHERES
#include "graphics/Drawables.h"
namespace Graphics {
	class RenderState;
}
#endif

namespace Graphics {
	class Renderer;
	class MeshObject;
} // namespace Graphics

class GeoPatchContext;
class GeoSphere;
class BasePatchJob;
class SQuadSplitResult;
class SSingleSplitResult;
struct SSplitResultData;

class GeoPatch {
public:
	GeoPatch(const RefCountedPtr<GeoPatchContext> &_ctx, GeoSphere *gs,
		const vector3d &v0_, const vector3d &v1_, const vector3d &v2_, const vector3d &v3_,
		const int depth, const GeoPatchID &ID_);

	~GeoPatch();

	void UpdateVBOs(Graphics::Renderer *renderer);

	int GetChildIdx(const GeoPatch *child) const
	{
		for (int i = 0; i < NUM_KIDS; i++) {
			if (m_kids[i].get() == child) return i;
		}
		abort();
		return -1;
	}

	void Render(Graphics::Renderer *renderer, const vector3d &campos, const matrix4x4d &modelView, const Graphics::Frustum &frustum);
	void RenderImmediate(Graphics::Renderer *renderer, const vector3d &campos, const matrix4x4d &modelView) const;
	void GatherRenderablePatches(std::vector<GeoPatch *> &visiblePatches, Graphics::Renderer *renderer, const vector3d &campos, const Graphics::Frustum &frustum);

	inline bool canBeMerged() const
	{
		bool merge = true;
		if (m_kids[0]) {
			for (int i = 0; i < NUM_KIDS; i++) {
				merge &= m_kids[i]->canBeMerged();
			}
		}
		merge &= !(HasJobRequest());
		return merge;
	}

	void LODUpdate(const vector3d &campos, const Graphics::Frustum &frustum);

	void RequestSinglePatch();
	void ReceiveHeightmaps(SQuadSplitResult *psr);
	void ReceiveHeightmap(const SSingleSplitResult *psr);
	void ReceiveHeightResult(const SSplitResultData &data);
	void ReceiveJobHandle(Job::Handle job);

	inline bool HasHeightData() const { return (m_patchVBOData != nullptr) && (m_patchVBOData->m_heights.get() != nullptr); }

	inline const vector3d &Centroid() const { return m_clipCentroid; }

	// used by GeoSphere so must be public
	inline void SetNeedToUpdateVBOs()
	{
		m_needUpdateVBOs = HasHeightData();
	}

private:

	inline bool NeedToUpdateVBOs() const
	{
		return m_needUpdateVBOs;
	}

	inline void ClearNeedToUpdateVBOs()
	{
		m_needUpdateVBOs = false;
	}

	inline void SetHasJobRequest()
	{
		m_hasJobRequest = true;
	}

	inline bool HasJobRequest() const
	{
		return m_hasJobRequest;
	}

	inline void ClearHasJobRequest()
	{
		m_hasJobRequest = false;
	}


private:
	static const int NUM_KIDS = 4;

	bool IsOverHorizon(const vector3d &camPos);

	RefCountedPtr<GeoPatchContext> m_ctx;
	struct Corners {
		Corners(const vector3d &v0_, const vector3d &v1_, const vector3d &v2_, const vector3d &v3_) :
			m_v0(v0_),
			m_v1(v1_),
			m_v2(v2_),
			m_v3(v3_)
		{}
		Corners() = delete;
		const vector3d m_v0, m_v1, m_v2, m_v3;
	}; 

	struct PatchVBOData {
		PatchVBOData() = delete;
		PatchVBOData(double* h, vector3f* n, Color3ub* c)
		{
			m_heights.reset(h);
			m_normals.reset(n);
			m_colors.reset(c);
		}
		~PatchVBOData() {
			m_heights.reset();
			m_normals.reset();
			m_colors.reset();
			m_patchMesh.reset();
		}
		std::unique_ptr<double[]> m_heights;
		std::unique_ptr<vector3f[]> m_normals;
		std::unique_ptr<Color3ub[]> m_colors;
		std::unique_ptr<Graphics::MeshObject> m_patchMesh;
	};

	std::unique_ptr<Corners> m_corners;
	std::unique_ptr<PatchVBOData> m_patchVBOData;
	std::unique_ptr<GeoPatch> m_kids[NUM_KIDS];

	vector3d m_clipCentroid;
	GeoSphere *m_geosphere;
	double m_splitLength; // rough length, is how near to the camera the m_clipCentroid should be before it must split
	double m_clipRadius;
	const GeoPatchID m_PatchID;
	Job::Handle m_job;
	Sint32 m_depth;
	uint8_t m_patchUpdateState;
	bool m_needUpdateVBOs;
	bool m_hasJobRequest;
#ifdef DEBUG_BOUNDING_SPHERES
	std::unique_ptr<Graphics::Drawables::Sphere3D> m_boundsphere;
#endif
};

#endif /* _GEOPATCH_H */
