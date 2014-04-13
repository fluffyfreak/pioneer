// Copyright © 2008-2014 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "libs.h"
#include "GeoPatchContext.h"
#include "GeoPatch.h"
#include "GeoPatchJobs.h"
#include "GeoSphere.h"
#include "perlin.h"
#include "Pi.h"
#include "RefCounted.h"
#include "graphics/Material.h"
#include "graphics/Renderer.h"
#include "graphics/Frustum.h"
#include "graphics/Graphics.h"
#include "graphics/VertexArray.h"
#include "MathUtil.h"
#include "vcacheopt/vcacheopt.h"
#include <deque>
#include <algorithm>

////////////////////////////////////////////////////////////////////////////////////
// cHaltonSequence2
//

namespace
{
   static const float kOneOverThree = float(1.0 / 3.0);
   static const float kOneOverFive  = float(1.0 / 5.0);
}

/// This calculates the Halton sequence incrementally
/// for a base 2/3 doublet.
class cHaltonSequence2
{
private:
//    cFXVector3 mPoint;  ///< Current sample point.
    float mX;
    float mY;
    
    Uint32 mBase2;
    Uint32 mBase3;
    
public:
    cHaltonSequence2() : 
        mBase2(0),
        mBase3(0),
        mX(0.0f),
        mY(0.0f)
    {}
    
    ///< Advance to next point in the sequence. Returns the index of this point. 
    int inc(vector2f &out)
	{
		/////////////////////////////////////
		// base 2
		Uint32 oldBase2 = mBase2;
		++mBase2;
		Uint32 diff = mBase2 ^ oldBase2;

		// bottom bit always changes, higher bits
		// change less frequently.
		float s = 0.5f;

		// diff will be of the form 0*1+, i.e. one bits up until the last carry.
		// expected iterations = 1 + 0.5 + 0.25 + ... = 2
		do
		{
			if (oldBase2 & 1) {
				mX -= s;
			} else {
				mX += s;
			}
        
			s *= 0.5f;
        
			diff = diff >> 1;
			oldBase2 = oldBase2 >> 1;
		}
		while (diff);

    
		/////////////////////////////////////
		// base 3: use 2 bits for each base 3 digit.
    
		Uint32 mask = 0x3;  // also the max base 3 digit
		Uint32 add  = 0x1;  // amount to add to force carry once digit==3
		s = kOneOverThree;

		++mBase3;

		// expected iterations: 1.5
		while (1)
		{
			if ((mBase3 & mask) == mask)
			{
				mBase3 += add;          // force carry into next 2-bit digit
				mY -= 2 * s;
            
				mask = mask << 2;
				add  = add  << 2;
            
				s *= kOneOverThree;
			}
			else 
			{
				mY += s;     // we know digit n has gone from a to a + 1
				break;
			}
		};

		out.x = mX;
		out.y = mY;

		return mBase2; // return the index of this sequence point
	}
	
	///< Move back to first point in the sequence (i.e. the origin.)
	void reset()
	{
		mBase2 = 0;
		mBase3 = 0;
		mX = 0.0f;
		mY = 0.0f;
	}
};

//for station waypoint interpolation
vector2f lerp(const vector2f& v1, const vector2f& v2, const float t)
{
	return t*v2 + (1.0-t)*v1;
}

// tri edge lengths
static const double GEOPATCH_SUBDIVIDE_AT_CAMDIST = 5.0;
#define GEOPATCH_MAX_DEPTH  15 + (2*Pi::detail.fracmult) //15

GeoPatch::GeoPatch(const RefCountedPtr<GeoPatchContext> &ctx_, GeoSphere *gs,
	const vector3d &v0_, const vector3d &v1_, const vector3d &v2_, const vector3d &v3_,
	const int depth, const GeoPatchID &ID_)
	: ctx(ctx_), v0(v0_), v1(v1_), v2(v2_), v3(v3_),
	heights(nullptr), normals(nullptr), colors(nullptr),
	m_numInstances(0), parent(nullptr), geosphere(gs),
	m_depth(depth), mPatchID(ID_),
	mHasJobRequest(false)
{
	for (int i=0; i<NUM_KIDS; ++i) {
		edgeFriend[i]	= NULL;
	}

	clipCentroid = (v0+v1+v2+v3) * 0.25;
	centroid = clipCentroid.Normalized();
	clipRadius = 0.0;
	clipRadius = std::max(clipRadius, (v0-clipCentroid).Length());
	clipRadius = std::max(clipRadius, (v1-clipCentroid).Length());
	clipRadius = std::max(clipRadius, (v2-clipCentroid).Length());
	clipRadius = std::max(clipRadius, (v3-clipCentroid).Length());
	double distMult;
	if (geosphere->GetSystemBody()->GetType() < SystemBody::TYPE_PLANET_ASTEROID) {
 		distMult = 10.0 / Clamp(depth, 1, 10);
 	} else {
 		distMult = 5.0 / Clamp(depth, 1, 5);
 	}
	m_roughLength = GEOPATCH_SUBDIVIDE_AT_CAMDIST / pow(2.0, depth) * distMult;
	m_needUpdateVBOs = false;
}

GeoPatch::~GeoPatch() {
	mHasJobRequest = false;

	for (int i=0; i<NUM_KIDS; i++) {
		if (edgeFriend[i]) edgeFriend[i]->NotifyEdgeFriendDeleted(this);
	}
	for (int i=0; i<NUM_KIDS; i++) {
		kids[i].reset();
	}
	heights.reset();
	normals.reset();
	colors.reset();
}
#pragma optimize("",off)
void GeoPatch::_UpdateVBOs(Graphics::Renderer *renderer) 
{
	if (m_needUpdateVBOs) {
		assert(renderer);
		m_needUpdateVBOs = false;

		//create buffer and upload data
		Graphics::VertexBufferDesc vbd;
		vbd.attrib[0].semantic = Graphics::ATTRIB_POSITION;
		vbd.attrib[0].format   = Graphics::ATTRIB_FORMAT_FLOAT3;
		vbd.attrib[1].semantic = Graphics::ATTRIB_NORMAL;
		vbd.attrib[1].format   = Graphics::ATTRIB_FORMAT_FLOAT3;
		vbd.attrib[2].semantic = Graphics::ATTRIB_DIFFUSE;
		vbd.attrib[2].format   = Graphics::ATTRIB_FORMAT_UBYTE4;
		vbd.numVertices = ctx->NUMVERTICES();
		vbd.usage = Graphics::BUFFER_USAGE_STATIC;
		m_vertexBuffer.reset(renderer->CreateVertexBuffer(vbd));

		GeoPatchContext::VBOVertex* vtxPtr = m_vertexBuffer->Map<GeoPatchContext::VBOVertex>(Graphics::BUFFER_MAP_WRITE);
		assert(m_vertexBuffer->GetDesc().stride == sizeof(GeoPatchContext::VBOVertex));

		const Sint32 edgeLen = ctx->edgeLen;
		const double frac = ctx->frac;
		const double *pHts = heights.get();
		const vector3f *pNorm = normals.get();
		const Color3ub *pColr = colors.get();
		for (Sint32 y=0; y<edgeLen; y++) {
			for (Sint32 x=0; x<edgeLen; x++) {
				const double height = *pHts;
				const double xFrac = double(x)*frac;
				const double yFrac = double(y)*frac;
				const vector3d p((GetSpherePoint(xFrac, yFrac) * (height + 1.0)) - clipCentroid);
				clipRadius = std::max(clipRadius, p.Length());
				vtxPtr->pos = vector3f(p);
				++pHts;	// next height

				const vector3f norma(pNorm->Normalized());
				vtxPtr->norm = norma;
				++pNorm; // next normal

				vtxPtr->col[0] = pColr->r;
				vtxPtr->col[1] = pColr->g;
				vtxPtr->col[2] = pColr->b;
				vtxPtr->col[3] = 255;
				++pColr; // next colour

				++vtxPtr; // next vertex
			}
		}
		m_vertexBuffer->Unmap();

		if( (GEOPATCH_MAX_DEPTH - m_depth) < 3 )
		{
			// allocate space for instances
			const Uint32 MaxInstances = ((edgeLen-1) * (edgeLen-1)) * 5;
			assert(MaxInstances > 0);
			instances.reset( new vector3f[MaxInstances] );
			m_numInstances = 0;

			// populate
			cHaltonSequence2 seq;
			vector2f outVec;
			pHts = heights.get();
			pNorm = normals.get();
			pColr = colors.get();
			for (int y=0; y<edgeLen-1; y++) {
				for (int x=0; x<edgeLen-1; x++) {
					// for each quad
					const double h0 = pHts[(x   + (y   * edgeLen))];
					const double h1 = pHts[(x+1 + (y   * edgeLen))];
					const double h2 = pHts[(x   + (y+1 * edgeLen))];
					const double h3 = pHts[(x+1 + (y+1 * edgeLen))];

					/*const vector3f& n0 = pNorm[(x   + (y   * edgeLen))];
					const vector3f& n1 = pNorm[(x+1 + (y   * edgeLen))];
					const vector3f& n2 = pNorm[(x   + (y+1 * edgeLen))];
					const vector3f& n3 = pNorm[(x+1 + (y+1 * edgeLen))];

					const Color3ub c0 = pColr[(x   + (y   * edgeLen))];
					const Color3ub c1 = pColr[(x+1 + (y   * edgeLen))];
					const Color3ub c2 = pColr[(x   + (y+1 * edgeLen))];
					const Color3ub c3 = pColr[(x+1 + (y+1 * edgeLen))];*/

					// p0--p1
					// |    |
					// p2--p3
					const vector3d p0 = (GetSpherePoint(x*frac,   y*frac)   * (h0 + 1.0)) - clipCentroid;
					const vector3d p1 = (GetSpherePoint(x+1*frac, y*frac)   * (h1 + 1.0)) - clipCentroid;
					const vector3d p2 = (GetSpherePoint(x*frac,   y+1*frac) * (h2 + 1.0)) - clipCentroid;
					const vector3d p3 = (GetSpherePoint(x+1*frac, y+1*frac) * (h3 + 1.0)) - clipCentroid;

					for( Uint32 iHal=0; iHal<5; ++iHal) {
						seq.inc(outVec);
						// bi-linear interpolation
						// x-axis first
						const vector3d i0 = MathUtil::mix<vector3d, double>(p0, p1, double(outVec.x));
						const vector3d i1 = MathUtil::mix<vector3d, double>(p2, p3, double(outVec.x));
						// y-axis and final position
						const vector3d pos = MathUtil::mix<vector3d, double>(i0, i1, double(outVec.y));
						instances.get()[m_numInstances++] = vector3f(pos);
					}
				}
			}
			assert(MaxInstances == m_numInstances);
		}
	}
}

void GeoPatch::Render(Graphics::Renderer *renderer, const vector3d &campos, const matrix4x4d &modelView, const Graphics::Frustum &frustum) {
	if (kids[0]) {
		for (int i=0; i<NUM_KIDS; i++) kids[i]->Render(renderer, campos, modelView, frustum);
	} else if (heights) {
		_UpdateVBOs(renderer);

		if (!frustum.TestPoint(clipCentroid, clipRadius))
			return;

		Graphics::Material *mat = geosphere->GetSurfaceMaterial();
		Graphics::RenderState *rs = geosphere->GetSurfRenderState();

		const vector3d relpos = clipCentroid - campos;
		renderer->SetTransform(modelView * matrix4x4d::Translation(relpos));

		Pi::statSceneTris += 2*(ctx->edgeLen-1)*(ctx->edgeLen-1);

		renderer->DrawBufferIndexed(m_vertexBuffer.get(), ctx->indices_list[determineIndexbuffer()].Get(), rs, mat);
	}
}

void GeoPatch::LODUpdate(const vector3d &campos) {
	// there should be no LODUpdate'ing when we have active split requests
	if(mHasJobRequest)
		return;

	bool canSplit = true;
	bool canMerge = bool(kids[0]);

	// always split at first level
	if (parent) {
		for (int i=0; i<NUM_EDGES; i++) {
			if (!edgeFriend[i]) {
				canSplit = false;
				break;
			} else if (edgeFriend[i]->m_depth < m_depth) {
				canSplit = false;
				break;
			}
		}
		const float centroidDist = (campos - centroid).Length();
		const bool errorSplit = (centroidDist < m_roughLength);
		if( !(canSplit && (m_depth < GEOPATCH_MAX_DEPTH) && errorSplit) ) {
			canSplit = false;
		}
	}

	if (canSplit) {
		if (!kids[0]) {
            assert(!mHasJobRequest);
            assert(!m_job.HasJob());
			mHasJobRequest = true;

			SQuadSplitRequest *ssrd = new SQuadSplitRequest(v0, v1, v2, v3, centroid.Normalized(), m_depth,
						geosphere->GetSystemBody()->GetPath(), mPatchID, ctx->edgeLen,
						ctx->frac, geosphere->GetTerrain());
			m_job = Pi::Jobs()->Queue(new QuadPatchJob(ssrd));
		} else {
			for (int i=0; i<NUM_KIDS; i++) {
				kids[i]->LODUpdate(campos);
			}
		}
	} else if (canMerge) {
		for (int i=0; i<NUM_KIDS; i++) {
			canMerge &= kids[i]->canBeMerged();
		}
		if( canMerge ) {
			for (int i=0; i<NUM_KIDS; i++) {
				kids[i].reset();
			}
		}
	}
}

void GeoPatch::RequestSinglePatch()
{
	if( !heights ) {
        assert(!mHasJobRequest);
        assert(!m_job.HasJob());
		mHasJobRequest = true;
		SSingleSplitRequest *ssrd = new SSingleSplitRequest(v0, v1, v2, v3, centroid.Normalized(), m_depth,
					geosphere->GetSystemBody()->GetPath(), mPatchID, ctx->edgeLen, ctx->frac, geosphere->GetTerrain());
		m_job = Pi::Jobs()->Queue(new SinglePatchJob(ssrd));
	}
}

void GeoPatch::ReceiveHeightmaps(SQuadSplitResult *psr)
{
	assert(NULL!=psr);
	if (m_depth<psr->depth()) {
		// this should work because each depth should have a common history
		const Uint32 kidIdx = psr->data(0).patchID.GetPatchIdx(m_depth+1);
		if( kids[kidIdx] ) {
			kids[kidIdx]->ReceiveHeightmaps(psr);
		} else {
			psr->OnCancel();
		}
	} else {
		assert(mHasJobRequest);
		const int nD = m_depth+1;
		for (int i=0; i<NUM_KIDS; i++)
		{
			assert(!kids[i]);
			const SQuadSplitResult::SSplitResultData& data = psr->data(i);
			assert(i==data.patchID.GetPatchIdx(nD));
			assert(0==data.patchID.GetPatchIdx(nD+1));
			kids[i].reset(new GeoPatch(ctx, geosphere,
				data.v0, data.v1, data.v2, data.v3,
				nD, data.patchID));
		}

		// hm.. edges. Not right to pass this
		// edgeFriend...
		kids[0]->edgeFriend[0] = GetEdgeFriendForKid(0, 0);
		kids[0]->edgeFriend[1] = kids[1].get();
		kids[0]->edgeFriend[2] = kids[3].get();
		kids[0]->edgeFriend[3] = GetEdgeFriendForKid(0, 3);
		kids[1]->edgeFriend[0] = GetEdgeFriendForKid(1, 0);
		kids[1]->edgeFriend[1] = GetEdgeFriendForKid(1, 1);
		kids[1]->edgeFriend[2] = kids[2].get();
		kids[1]->edgeFriend[3] = kids[0].get();
		kids[2]->edgeFriend[0] = kids[1].get();
		kids[2]->edgeFriend[1] = GetEdgeFriendForKid(2, 1);
		kids[2]->edgeFriend[2] = GetEdgeFriendForKid(2, 2);
		kids[2]->edgeFriend[3] = kids[3].get();
		kids[3]->edgeFriend[0] = kids[0].get();
		kids[3]->edgeFriend[1] = kids[2].get();
		kids[3]->edgeFriend[2] = GetEdgeFriendForKid(3, 2);
		kids[3]->edgeFriend[3] = GetEdgeFriendForKid(3, 3);
		kids[0]->parent = kids[1]->parent = kids[2]->parent = kids[3]->parent = this;

		for (int i=0; i<NUM_KIDS; i++)
		{
			const SQuadSplitResult::SSplitResultData& data = psr->data(i);
			kids[i]->heights.reset(data.heights);
			kids[i]->normals.reset(data.normals);
			kids[i]->colors.reset(data.colors);
		}
		for (int i=0; i<NUM_EDGES; i++) { if(edgeFriend[i]) edgeFriend[i]->NotifyEdgeFriendSplit(this); }
		for (int i=0; i<NUM_KIDS; i++) {
			kids[i]->UpdateVBOs();
		}
		mHasJobRequest = false;
	}
}

void GeoPatch::ReceiveHeightmap(const SSingleSplitResult *psr)
{
	assert(nullptr == parent);
	assert(nullptr != psr);
	assert(mHasJobRequest);
	{
		const SSingleSplitResult::SSplitResultData& data = psr->data();
		heights.reset(data.heights);
		normals.reset(data.normals);
		colors.reset(data.colors);
	}
	mHasJobRequest = false;
}
