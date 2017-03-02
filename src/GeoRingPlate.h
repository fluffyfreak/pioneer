#ifndef _GEORINGPLATE_H
#define _GEORINGPLATE_H

#include "vector3.h"
#include "Camera.h"
#include "MathUtil.h"
#include "terrain/Terrain.h"

extern int GEOPATCH_EDGELEN;
#define ATMOSPHERE_RADIUS 1.015

namespace Graphics { 
	class Renderer; 
	class RenderState;
	class Material;
}
class SystemBody;
class GeoPatchContext;

class SQuadPlateRequest;
class SQuadPlateResult;
class SSinglePlateResult;

#define NUM_VBE 2

// ********************************************************************************
//
// ********************************************************************************
class GeoPlate
{
public:
	RefCountedPtr<GeoPatchContext> ctx;
	double ang[NUM_VBE];
	double m_halfLen;
	double m_yoffset;		// offset from the centre i.e. newyoffset = (m_yoffset + m_halfLen);
	std::unique_ptr<vector3f[]> normals;
	std::unique_ptr<Color3ub[]> colors;
	GeoRing *geoRing;
	std::unique_ptr<Graphics::VertexBuffer> m_vertexBuffer;
	double m_roughLength;
	vector3d clipCentroid;
	double clipRadius;
	vector3d vbe[NUM_VBE];
	std::unique_ptr<GeoPlate> kids[4];
	int m_depth;
	bool m_needUpdateVBOs;
	std::unique_ptr<double[]> heights;

	const GeoPlateID mPatchID;
	Job::Handle m_job;
	bool mHasJobRequest;

	// params
	// v0, v1 - define points on the line describing the loop of the ring/orbital.
	// depth - 0 is the topmost plate with each depth+1 describing it's depth within the tree.
	GeoPlate(const RefCountedPtr<GeoPatchContext> &ctx_, GeoRing *geoRingPtr, const double halfLength, const vector3d &startVBE, const vector3d &endVBE, 
		const double startAng, const double endAng, const double yoffset, const int depth, const GeoPlateID ID);

	virtual ~GeoPlate() {
	}

	bool HasData() const { return heights.get() != nullptr; }

	inline vector3d GetSurfacePointCyl(const double x, const double y, const double halfLength) const {
		double theta = MathUtil::mix(ang[1], ang[0], x);
		
		const vector3d topEndEdge(sin(theta), m_yoffset + (halfLength * 0.5), cos(theta));		// vertices at top edge of circle
		const vector3d bottomEndEdge(sin(theta), m_yoffset - (halfLength * 0.5), cos(theta));	// vertices at bottom edge of circle
		
		const vector3d res = MathUtil::mix( bottomEndEdge, topEndEdge, y );
		return res;
	}

	void NeedToUpdateVBOs()	{ m_needUpdateVBOs = true; }
	void UpdateVBOs();

	bool canBeMerged() const {
		bool merge = true;
		if (kids[0]) {
			for (int i=0; i<4; i++) {
				merge &= kids[i]->canBeMerged();
			}
		}
		merge &= !(mHasJobRequest);
		return merge;
	}

	void ReceiveHeightmaps(SQuadPlateResult *psr);
	void ReceiveHeightmap(const SSinglePlateResult *psr);

	inline void ReceiveJobHandle(Job::Handle job)
	{
		assert(!m_job.HasJob());
		m_job = static_cast<Job::Handle&&>(job);
	}

	void RequestSinglePatch();
	
	void Render(Graphics::Renderer *renderer, const vector3d &campos, const matrix4x4d &modelView, const Graphics::Frustum &frustum);
	void LODUpdate(vector3d &campos);
};

#endif /* _GEORINGPLATE_H */
