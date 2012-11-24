// Copyright © 2008-2012 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef __GEOPATCH_H__
#define __GEOPATCH_H__

#include "GeoPatchID.h"

// fwd decl's
namespace Graphics {
	class VertexBuffer;
};

class GeoPatch
{
private:
	////////////////////////////////////////////////////////////////
	// private members
	static const uint32_t NUM_EDGES = 4;
	GeoPatch *edgeFriend[NUM_EDGES];

	static const uint32_t NUM_KIDS = NUM_EDGES;
	GeoPatch *kids[NUM_KIDS];

	const GeoPatchContext &mContext;
	GeoMesh *mpGeoSphere;
	const vector3f mV0, mV1, mV2, mV3;	// corner vertices for the patch
	const vector3f mClipCentroid;
	vector3f mCentroid;
	const uint32_t mDepth;
	float mClipRadius;
	float mRoughLength;
	const GeoPatchID mPatchID;

	Graphics::VertexBuffer *mVBO;
	GLuint mHeightmap;
	bool mHasSplitRequest;

public:
	////////////////////////////////////////////////////////////////
	// public members
	GeoPatch *parent;

	////////////////////////////////////////////////////////////////
	// public methods

	// constructor
	GeoPatch(const GeoPatchContext &context_, GeoMesh *pGeoSphere_, 
		const vector3f &v0_, const vector3f &v1_, const vector3f &v2_, const vector3f &v3_, 
		const uint32_t depth_, const GeoPatchID &ID_);

	// destructor
	~GeoPatch();

	// in patch surface coords, [0,1]
	inline vector3f GetSpherePoint(const float x, const float y) const {
		return ((mV0 + x*(1.0f-y)*(mV1-mV0) + x*y*(mV2-mV0) + (1.0f-x)*y*(mV3-mV0))).Normalized();
	}

	// Generates full-detail vertices, and also non-edge normals and colors
	void GenerateMesh();

	void ReceiveHeightmaps(const SSplitResult *psr);
	void ReceiveHeightmapTex(const GLuint tex);

	inline void OnEdgeFriendChanged(const int edge, GeoPatch *e) {
		edgeFriend[edge] = e;
	}

	inline void NotifyEdgeFriendSplit(GeoPatch *e) {
		if (!kids[0]) {
			return;
		}
		const int idx = GetEdgeIdxOf(e);
		const int we_are = e->GetEdgeIdxOf(this);
		// match e's new kids to our own... :/
		kids[idx]->OnEdgeFriendChanged(idx, e->kids[(we_are+1)%4]);
		kids[(idx+1)%4]->OnEdgeFriendChanged(idx, e->kids[we_are]);
	}

	inline void NotifyEdgeFriendDeleted(const GeoPatch *e) {
		const int idx = GetEdgeIdxOf(e);
		assert(idx>=0 && idx<4);
		edgeFriend[idx] = nullptr;
	}

	inline int GetEdgeIdxOf(const GeoPatch *e) const {
		for (int i=0; i<4; i++) {
			if (edgeFriend[i] == e) {
				return i;
			}
		}
		assert(false);
		return -1;
	}

	inline GeoPatch *GetEdgeFriendForKid(const int kid, const int edge) const {
		const GeoPatch *e = edgeFriend[edge];
		if (!e) {
			return nullptr;
		}
		const int we_are = e->GetEdgeIdxOf(this);
		assert(we_are>=0 && we_are<4);
		// neighbour patch has not split yet (is at depth of this patch), so kids of this patch do
		// not have same detail level neighbours yet
		if (edge == kid) {
			return e->kids[(we_are+1)%4];
		} 
		return e->kids[we_are];
	}

	inline GLuint determineIndexbuffer() const {
		return // index buffers are ordered by edge resolution flags
			(edgeFriend[0] ? 1u : 0u) |
			(edgeFriend[1] ? 2u : 0u) |
			(edgeFriend[2] ? 4u : 0u) |
			(edgeFriend[3] ? 8u : 0u);
	}

	void Render(const vector3f &campos, const Graphics::Frustum &frustum);

	void LODUpdate(const vector3f &campos);
};

#endif //__GEOPATCH_H__