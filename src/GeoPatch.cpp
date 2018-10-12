// Copyright © 2008-2018 Pioneer Developers. See AUTHORS.txt for details
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
#include "scenegraph/Model.h"
#include "MathUtil.h"
#include "Sphere.h"
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

GeoPatch::GeoPatch(const RefCountedPtr<GeoPatchContext> &ctx_, GeoSphere *gs,
	const vector3d &v0_, const vector3d &v1_, const vector3d &v2_, const vector3d &v3_,
	const int depth, const GeoPatchID &ID_)
	: ctx(ctx_), v0(v0_), v1(v1_), v2(v2_), v3(v3_),
	heights(nullptr), normals(nullptr), colors(nullptr),
	m_numInstances(0), parent(nullptr), geosphere(gs),
	m_depth(depth), mPatchID(ID_),
	mHasJobRequest(false)
{

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
		kids[i].reset();
	}
	heights.reset();
	normals.reset();
	colors.reset();
}

void GeoPatch::UpdateVBOs(Graphics::Renderer *renderer)
{
	PROFILE_SCOPED()
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
		vbd.attrib[3].semantic = Graphics::ATTRIB_UV0;
		vbd.attrib[3].format   = Graphics::ATTRIB_FORMAT_FLOAT2;
		vbd.numVertices = ctx->NUMVERTICES();
		vbd.usage = Graphics::BUFFER_USAGE_STATIC;
		m_vertexBuffer.reset(renderer->CreateVertexBuffer(vbd));

		GeoPatchContext::VBOVertex* VBOVtxPtr = m_vertexBuffer->Map<GeoPatchContext::VBOVertex>(Graphics::BUFFER_MAP_WRITE);
		assert(m_vertexBuffer->GetDesc().stride == sizeof(GeoPatchContext::VBOVertex));

		const Sint32 edgeLen = ctx->GetEdgeLen();
		const double frac = ctx->GetFrac();
		const double *pHts = heights.get();
		const vector3f *pNorm = normals.get();
		const Color3ub *pColr = colors.get();

		double minh = DBL_MAX;

		// ----------------------------------------------------
		// inner loops
		for (Sint32 y = 1; y<edgeLen-1; y++) {
			for (Sint32 x = 1; x<edgeLen-1; x++) {
				const double height = *pHts;
				minh = std::min(height, minh);
				const double xFrac = double(x - 1) * frac;
				const double yFrac = double(y - 1) * frac;
				const vector3d p((GetSpherePoint(xFrac, yFrac) * (height + 1.0)) - clipCentroid);
				clipRadius = std::max(clipRadius, p.Length());

				GeoPatchContext::VBOVertex* vtxPtr = &VBOVtxPtr[x + (y*edgeLen)];
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

				// uv coords
				vtxPtr->uv.x = 1.0f - xFrac;
				vtxPtr->uv.y = yFrac;

				++vtxPtr; // next vertex
			}
		}
		const double minhScale = (minh + 1.0) * 0.999995;
		// ----------------------------------------------------
		const Sint32 innerLeft = 1;
		const Sint32 innerRight = edgeLen - 2;
		const Sint32 outerLeft = 0;
		const Sint32 outerRight = edgeLen - 1;
		// vertical edges
		// left-edge
		for (Sint32 y = 1; y < edgeLen - 1; y++) {
			const Sint32 x = innerLeft-1;
			const double xFrac = double(x - 1) * frac;
			const double yFrac = double(y - 1) * frac;
			const vector3d p((GetSpherePoint(xFrac, yFrac) * minhScale) - clipCentroid);

			GeoPatchContext::VBOVertex* vtxPtr = &VBOVtxPtr[outerLeft + (y*edgeLen)];
			GeoPatchContext::VBOVertex* vtxInr = &VBOVtxPtr[innerLeft + (y*edgeLen)];
			vtxPtr->pos = vector3f(p);
			vtxPtr->norm = vtxInr->norm;
			vtxPtr->col = vtxInr->col;
			vtxPtr->uv = vtxInr->uv;
		}
		// right-edge
		for (Sint32 y = 1; y < edgeLen - 1; y++) {
			const Sint32 x = innerRight+1;
			const double xFrac = double(x - 1) * frac;
			const double yFrac = double(y - 1) * frac;
			const vector3d p((GetSpherePoint(xFrac, yFrac) * minhScale) - clipCentroid);

			GeoPatchContext::VBOVertex* vtxPtr = &VBOVtxPtr[outerRight + (y*edgeLen)];
			GeoPatchContext::VBOVertex* vtxInr = &VBOVtxPtr[innerRight + (y*edgeLen)];
			vtxPtr->pos = vector3f(p);
			vtxPtr->norm = vtxInr->norm;
			vtxPtr->col = vtxInr->col;
			vtxPtr->uv = vtxInr->uv;
		}
		// ----------------------------------------------------
		const Sint32 innerTop = 1;
		const Sint32 innerBottom = edgeLen - 2;
		const Sint32 outerTop = 0;
		const Sint32 outerBottom = edgeLen - 1;
		// horizontal edges
		// top-edge
		for (Sint32 x = 1; x < edgeLen - 1; x++)
		{
			const Sint32 y = innerTop-1;
			const double xFrac = double(x - 1) * frac;
			const double yFrac = double(y - 1) * frac;
			const vector3d p((GetSpherePoint(xFrac, yFrac) * minhScale) - clipCentroid);

			GeoPatchContext::VBOVertex* vtxPtr = &VBOVtxPtr[x + (outerTop*edgeLen)];
			GeoPatchContext::VBOVertex* vtxInr = &VBOVtxPtr[x + (innerTop*edgeLen)];
			vtxPtr->pos = vector3f(p);
			vtxPtr->norm = vtxInr->norm;
			vtxPtr->col = vtxInr->col;
			vtxPtr->uv = vtxInr->uv;
		}
		// bottom-edge
		for (Sint32 x = 1; x < edgeLen - 1; x++)
		{
			const Sint32 y = innerBottom+1;
			const double xFrac = double(x - 1) * frac;
			const double yFrac = double(y - 1) * frac;
			const vector3d p((GetSpherePoint(xFrac, yFrac) * minhScale) - clipCentroid);

			GeoPatchContext::VBOVertex* vtxPtr = &VBOVtxPtr[x + (outerBottom * edgeLen)];
			GeoPatchContext::VBOVertex* vtxInr = &VBOVtxPtr[x + (innerBottom * edgeLen)];
			vtxPtr->pos = vector3f(p);
			vtxPtr->norm = vtxInr->norm;
			vtxPtr->col = vtxInr->col;
			vtxPtr->uv = vtxInr->uv;
		}
		// ----------------------------------------------------
		// corners
		{
			// top left
			GeoPatchContext::VBOVertex* tarPtr = &VBOVtxPtr[0];
			GeoPatchContext::VBOVertex* srcPtr = &VBOVtxPtr[1];
			(*tarPtr) = (*srcPtr);
		}
		{
			// top right
			GeoPatchContext::VBOVertex* tarPtr = &VBOVtxPtr[(edgeLen - 1)];
			GeoPatchContext::VBOVertex* srcPtr = &VBOVtxPtr[(edgeLen - 2)];
			(*tarPtr) = (*srcPtr);
		}
		{
			// bottom left
			GeoPatchContext::VBOVertex* tarPtr = &VBOVtxPtr[(edgeLen - 1) * edgeLen];
			GeoPatchContext::VBOVertex* srcPtr = &VBOVtxPtr[(edgeLen - 2) * edgeLen];
			(*tarPtr) = (*srcPtr);
		}
		{
			// bottom right
			GeoPatchContext::VBOVertex* tarPtr = &VBOVtxPtr[(edgeLen - 1) + ((edgeLen - 1) * edgeLen)];
			GeoPatchContext::VBOVertex* srcPtr = &VBOVtxPtr[(edgeLen - 1) + ((edgeLen - 2) * edgeLen)];
			(*tarPtr) = (*srcPtr);
		}

		// ----------------------------------------------------
		// end of mapping
		m_vertexBuffer->Unmap();

		// Don't need this anymore so throw it away
		normals.reset();
		colors.reset();

		if( (geosphere->GetMaxDepth() - m_depth) < 12 )
		{
			// allocate space for instances
			instances.reset( new vector3f[(edgeLen-2)*(edgeLen-2)] );
			m_numInstances = 0;

			const double *pHts = heights.get();
			const vector3f *pNorm = normals.get();
			const Color3ub *pColr = colors.get();
			for (Sint32 y=1; y<edgeLen-1; y++) {
				for (Sint32 x=1; x<edgeLen-1; x++) {
					const double height = *pHts;
					const double xFrac = double(x - 1) * frac;
					const double yFrac = double(y - 1) * frac;
					const vector3d p((GetSpherePoint(xFrac, yFrac) * (height + 1.0)) - clipCentroid);
					clipRadius = std::max(clipRadius, p.Length());
					instances.get()[m_numInstances++] = vector3f(p);
					++pHts;	// next height
				}
			}

			// populate
			/*cHaltonSequence2 seq;
			vector2f outVec;
			pHts = heights.get();
			for( Uint32 iHal=0; iHal<5; ++iHal) 
			{
				seq.inc(outVec);
				const Sint32 x = Sint32(outVec.x * (edgeLen-1));
				const Sint32 y = Sint32(outVec.y * (edgeLen-1));

				// for each quad
				const double h0 = pHts[(x   + (y   * edgeLen))];
				const double h1 = pHts[(x+1 + (y   * edgeLen))];
				const double h2 = pHts[(x   + (y+1 * edgeLen))];
				const double h3 = pHts[(x+1 + (y+1 * edgeLen))];

				// p0--p1
				// |    |
				// p2--p3
				const vector3d p0 = (GetSpherePoint(x*frac,   y*frac)   * (h0 + 1.0)) - clipCentroid;
				const vector3d p1 = (GetSpherePoint(x+1*frac, y*frac)   * (h1 + 1.0)) - clipCentroid;
				const vector3d p2 = (GetSpherePoint(x*frac,   y+1*frac) * (h2 + 1.0)) - clipCentroid;
				const vector3d p3 = (GetSpherePoint(x+1*frac, y+1*frac) * (h3 + 1.0)) - clipCentroid;

				// bi-linear interpolation
				// x-axis first
				const vector3d i0 = MathUtil::mix<vector3d, double>(p0, p1, double(outVec.x));
				const vector3d i1 = MathUtil::mix<vector3d, double>(p2, p3, double(outVec.x));
				// y-axis and final position
				const vector3d pos = MathUtil::mix<vector3d, double>(i0, i1, double(outVec.y));
				instances.get()[m_numInstances++] = vector3f(pos);
			}*/
			assert(MaxInstances == m_numInstances);
		}
#ifdef DEBUG_BOUNDING_SPHERES
		RefCountedPtr<Graphics::Material> mat(Pi::renderer->CreateMaterial(Graphics::MaterialDescriptor()));
		m_boundsphere.reset( new Graphics::Drawables::Sphere3D(Pi::renderer, mat, Pi::renderer->CreateRenderState(Graphics::RenderStateDesc()), 0, clipRadius) );
#endif
	}
}

// the default sphere we do the horizon culling against
static const SSphere s_sph;
#pragma optimize("",off)
void GeoPatch::Render(Graphics::Renderer *renderer, const vector3d &campos, const matrix4x4d &modelView, const Graphics::Frustum &frustum)
{
	PROFILE_SCOPED()
	// must update the VBOs to calculate the clipRadius...
	UpdateVBOs(renderer);
	// ...before doing the furstum culling that relies on it.
	if (!frustum.TestPoint(clipCentroid, clipRadius))
		return; // nothing below this patch is visible

	// only want to horizon cull patches that can actually be over the horizon!
	const vector3d camDir(campos - clipCentroid);
	const vector3d camDirNorm(camDir.Normalized());
	const vector3d cenDir(clipCentroid.Normalized());
	const double dotProd = camDirNorm.Dot(cenDir);

	if( dotProd < 0.25 && (camDir.LengthSqr() > (clipRadius*clipRadius)) ) {
		SSphere obj;
		obj.m_centre = clipCentroid;
		obj.m_radius = clipRadius;

		if( !s_sph.HorizonCulling(campos, obj) ) {
			return; // nothing below this patch is visible
		}
	}

	if (kids[0]) {
		for (int i=0; i<NUM_KIDS; i++) kids[i]->Render(renderer, campos, modelView, frustum);
	} else if (heights) {
		RefCountedPtr<Graphics::Material> mat = geosphere->GetSurfaceMaterial();
		Graphics::RenderState *rs = geosphere->GetSurfRenderState();

		const vector3d relpos = clipCentroid - campos;
		renderer->SetTransform(modelView * matrix4x4d::Translation(relpos));

		Pi::statSceneTris += (ctx->GetNumTris());
		++Pi::statNumPatches;

		// per-patch detail texture scaling value
		geosphere->GetMaterialParameters().patchDepth = m_depth;

		renderer->DrawBufferIndexed(m_vertexBuffer.get(), ctx->GetIndexBuffer(), rs, mat.Get());
		SceneGraph::Model *pModel = ctx->GetModelLibrary();
		if(m_numInstances>0 && pModel)
		{
			renderer->SetTransform(matrix4x4f::Identity());
			matrix4x4f mv;
			matrix4x4dtof(modelView * matrix4x4d::Translation(relpos), mv);

			std::vector<matrix4x4f> transforms(m_numInstances);
			for(Uint32 in=0; in<m_numInstances; in++) {
				transforms[in] = mv * matrix4x4f::Translation(instances[in]) * matrix4x4f::ScaleMatrix(0.00000006537);
			}
			pModel->Render(transforms);
		}

#ifdef DEBUG_BOUNDING_SPHERES
		if(m_boundsphere.get()) {
			renderer->SetWireFrameMode(true);
			m_boundsphere->Draw(renderer);
			renderer->SetWireFrameMode(false);
		}
#endif
		renderer->GetStats().AddToStatCount(Graphics::Stats::STAT_PATCHES, 1);
	}
}

void GeoPatch::LODUpdate(const vector3d &campos, const Graphics::Frustum &frustum)
{
	// there should be no LOD update when we have active split requests
	if(mHasJobRequest)
		return;

	bool canSplit = true;
	bool canMerge = bool(kids[0]);

	// always split at first level
	double centroidDist = DBL_MAX;
	if (parent) {
		centroidDist = (campos - centroid).Length();
		const bool errorSplit = (centroidDist < m_roughLength);
		if( !(canSplit && (m_depth < std::min(GEOPATCH_MAX_DEPTH, geosphere->GetMaxDepth())) && errorSplit) ) {
			canSplit = false;
		}
	}

	if (canSplit) {
		if (!kids[0]) {
			// Test if this patch is visible
			if (!frustum.TestPoint(clipCentroid, clipRadius))
				return; // nothing below this patch is visible

			// only want to horizon cull patches that can actually be over the horizon!
			const vector3d camDir(campos - clipCentroid);
			const vector3d camDirNorm(camDir.Normalized());
			const vector3d cenDir(clipCentroid.Normalized());
			const double dotProd = camDirNorm.Dot(cenDir);

			if (dotProd < 0.25 && (camDir.LengthSqr() >(clipRadius*clipRadius))) {
				SSphere obj;
				obj.m_centre = clipCentroid;
				obj.m_radius = clipRadius;

				if (!s_sph.HorizonCulling(campos, obj)) {
					return; // nothing below this patch is visible
				}
			}

			// we can see this patch so submit the jobs!
			assert(!mHasJobRequest);
			mHasJobRequest = true;

			SQuadSplitRequest *ssrd = new SQuadSplitRequest(v0, v1, v2, v3, centroid.Normalized(), m_depth,
						geosphere->GetSystemBody()->GetPath(), mPatchID, ctx->GetEdgeLen()-2,
						ctx->GetFrac(), geosphere->GetTerrain());

			// add to the GeoSphere to be processed at end of all LODUpdate requests
			geosphere->AddQuadSplitRequest(centroidDist, ssrd, this);
		} else {
			for (int i=0; i<NUM_KIDS; i++) {
				kids[i]->LODUpdate(campos, frustum);
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
		mHasJobRequest = true;
		SSingleSplitRequest *ssrd = new SSingleSplitRequest(v0, v1, v2, v3, centroid.Normalized(), m_depth,
					geosphere->GetSystemBody()->GetPath(), mPatchID, ctx->GetEdgeLen()-2, ctx->GetFrac(), geosphere->GetTerrain());
		m_job = Pi::GetAsyncJobQueue()->Queue(new SinglePatchJob(ssrd));
	}
}

void GeoPatch::ReceiveHeightmaps(SQuadSplitResult *psr)
{
	PROFILE_SCOPED()
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
		kids[0]->parent = kids[1]->parent = kids[2]->parent = kids[3]->parent = this;

		for (int i=0; i<NUM_KIDS; i++)
		{
			const SQuadSplitResult::SSplitResultData& data = psr->data(i);
			kids[i]->heights.reset(data.heights);
			kids[i]->normals.reset(data.normals);
			kids[i]->colors.reset(data.colors);
		}
		for (int i=0; i<NUM_KIDS; i++) {
			kids[i]->NeedToUpdateVBOs();
		}
		mHasJobRequest = false;
	}
}

void GeoPatch::ReceiveHeightmap(const SSingleSplitResult *psr)
{
	PROFILE_SCOPED()
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

void GeoPatch::ReceiveJobHandle(Job::Handle job)
{
	assert(!m_job.HasJob());
	m_job = static_cast<Job::Handle&&>(job);
}
