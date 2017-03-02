#include "libs.h"
#include "GeoPatchContext.h"
#include "GeoPatchID.h"
#include "GeoRing.h"
#include "GeoRingPlate.h"
#include "GeoRingJobs.h"
#include "perlin.h"
#include "Pi.h"
#include "galaxy/StarSystem.h"
#include "graphics/Frustum.h"
#include "graphics/Renderer.h"
#include "graphics/Material.h"
#include "graphics/Texture.h"
#include "graphics/TextureBuilder.h"

#include "profiler/Profiler.h"

#include "vcacheopt/vcacheopt.h"

#define DUMP_TO_TEXTURE 0

#if DUMP_TO_TEXTURE
#include "FileSystem.h"
#include "PngWriter.h"
#include "graphics/opengl/TextureGL.h"
void textureDump(const char* destFile, const int width, const int height, const int bpp, const Color3ub* buf)
{
	const std::string dir = "orbital_textures";
	FileSystem::userFiles.MakeDirectory(dir);
	const std::string fname = FileSystem::JoinPathBelow(dir, destFile);

	// pad rows to 4 bytes, which is the default row alignment for OpenGL
	//const int stride = (3*width + 3) & ~3;
	const int stride = width * bpp;

	write_png(FileSystem::userFiles, fname, &buf[0].r, width, height, stride, bpp);

	printf("texture %s saved\n", fname.c_str());
}
#endif

// tri edge lengths
#define GEOPLATE_SUBDIVIDE_AT_CAMDIST	1.5 //5.0
#define GEOPLATE_MAX_DEPTH	15

static const int GEOPLATE_MAX_EDGELEN = 55;

#define lerp(t, a, b) ( a + t * (b - a) )

// must be odd numbers
static const int detail_edgeLen[5] = {
	7, 15, 25, 35, 55
};

inline void SetColour(Color3ub *r, const vector3d &v) { 
	r->r=static_cast<unsigned char>(Clamp(v.x*255.0, 0.0, 255.0)); 
	r->g=static_cast<unsigned char>(Clamp(v.y*255.0, 0.0, 255.0)); 
	r->b=static_cast<unsigned char>(Clamp(v.z*255.0, 0.0, 255.0));
}

// ********************************************************************************
//
// ********************************************************************************

// --------------------------------------
//static 
// --------------------------------------
static bool s_bUseWireframe = false;
static bool s_bCanUpdate = true;

// params
// v0, v1 - define points on the line describing the loop of the ring/orbital.
// depth - 0 is the topmost plate with each depth+1 describing it's depth within the tree.
GeoPlate::GeoPlate(const RefCountedPtr<GeoPatchContext> &ctx_, GeoRing *geoRingPtr, const double halfLength, const vector3d &startVBE, const vector3d &endVBE, 
	const double startAng, const double endAng, const double yoffset, const int depth, const GeoPlateID ID) 
	: ctx(ctx_), normals(nullptr), colors(nullptr), geoRing(geoRingPtr), m_needUpdateVBOs(false), mPatchID(ID), mHasJobRequest(false)
{
	ang[0] = startAng; 
	ang[1] = endAng;
	m_halfLen = halfLength;
	m_yoffset = yoffset;
	clipCentroid = GetSurfacePointCyl(0.5, 0.5, m_halfLen);
	clipRadius = 0;
	vector3d vcorners[4] = {
		GetSurfacePointCyl(0.0, 0.0, m_halfLen),
		GetSurfacePointCyl(1.0, 0.0, m_halfLen),
		GetSurfacePointCyl(0.0, 1.0, m_halfLen),
		GetSurfacePointCyl(1.0, 1.0, m_halfLen)
	};
	for (int i=0; i<4; i++) {
		clipRadius = std::max(clipRadius, (vcorners[i]-clipCentroid).Length());
	}
	m_roughLength = GEOPLATE_SUBDIVIDE_AT_CAMDIST / pow(2.0, depth);

	vbe[0] = startVBE;
	vbe[1] = endVBE;
	m_depth = depth;
	m_needUpdateVBOs = false;
}

void GeoPlate::UpdateVBOs() 
{
	if (m_needUpdateVBOs) 
	{
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
		m_vertexBuffer.reset(Pi::renderer->CreateVertexBuffer(vbd));

		GeoPatchContext::VBOVertex* VBOVtxPtr = m_vertexBuffer->Map<GeoPatchContext::VBOVertex>(Graphics::BUFFER_MAP_WRITE);
		assert(m_vertexBuffer->GetDesc().stride == sizeof(GeoPatchContext::VBOVertex));

		const Sint32 edgeLen = ctx->GetEdgeLen();
		const double frac = ctx->GetFrac();
		const double *hgt = heights.get();
		const vector3f *pNorm = normals.get();
		const Color3ub *pColr = colors.get();

		const vector3d topEnd(0.0, m_yoffset + m_halfLen, 0.0);		// vertices at top edge of circle
		const vector3d btmEnd(0.0, m_yoffset - m_halfLen, 0.0);		// vertices at bottom edge of circle

		const Sint32 innerTop = 1;
		const Sint32 innerBottom = edgeLen - 2;
		const Sint32 innerLeft = 1;
		const Sint32 innerRight = edgeLen - 2;

		const Sint32 outerTop = 0;
		const Sint32 outerBottom = edgeLen - 1;
		const Sint32 outerLeft = 0;
		const Sint32 outerRight = edgeLen - 1;

		double maxh = DBL_MIN;

		// ----------------------------------------------------
		// inner loops
		for (Sint32 y = 1; y<edgeLen-1; y++) 
		{
			const double yFrac = double(y - 1) * frac;
			// point along z-axis by zfrac amount
			const vector3d axisPt = MathUtil::mix( btmEnd, topEnd, yFrac );
			for (Sint32 x = 1; x<edgeLen-1; x++) 
			{
				const double height = *hgt;
				maxh = std::max(height, maxh);
				const double xFrac = double(x - 1) * frac;

				const vector3d pCyl = GetSurfacePointCyl(xFrac, yFrac, m_halfLen);// find point on _surface_ of the cylinder
				const vector3d cDir = (pCyl - axisPt).Normalized();// vector from axis to point-on-surface
				vector3f p = vector3f((pCyl + (cDir * height)) - clipCentroid);// vertex is moved in direction of point-in-axis FROM point-on-surface by height.
				clipRadius = std::max(clipRadius, pCyl.Length());

				GeoPatchContext::VBOVertex* vtxPtr = &VBOVtxPtr[x + (y*edgeLen)];
				vtxPtr->pos = vector3f(p);
				++hgt;	// next height

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
			}
		}
		const double maxhScale = (maxh + 1.0) * 1.00001;//0.99999;
													
		// ----------------------------------------------------
		// vertical edges
		// left-edge
		for (Sint32 y = 1; y < edgeLen - 1; y++) {
			const Sint32 x = innerLeft;
			const double xFrac = double(x - 1) * frac;
			const double yFrac = double(y - 1) * frac;
			const vector3d p((GetSurfacePointCyl(xFrac, yFrac, m_halfLen) * maxhScale) - clipCentroid);

			GeoPatchContext::VBOVertex* vtxPtr = &VBOVtxPtr[outerLeft + (y*edgeLen)];
			GeoPatchContext::VBOVertex* vtxInr = &VBOVtxPtr[innerLeft + (y*edgeLen)];
			vtxPtr->pos = vector3f(p);
			vtxPtr->norm = vtxInr->norm;
			vtxPtr->col = vtxInr->col;
			vtxPtr->uv = vtxInr->uv;
		}
		// right-edge
		for (Sint32 y = 1; y < edgeLen - 1; y++) {
			const Sint32 x = innerRight;
			const double xFrac = double(x - 1) * frac;
			const double yFrac = double(y - 1) * frac;
			const vector3d p((GetSurfacePointCyl(xFrac, yFrac, m_halfLen) * maxhScale) - clipCentroid);

			GeoPatchContext::VBOVertex* vtxPtr = &VBOVtxPtr[outerRight + (y*edgeLen)];
			GeoPatchContext::VBOVertex* vtxInr = &VBOVtxPtr[innerRight + (y*edgeLen)];
			vtxPtr->pos = vector3f(p);
			vtxPtr->norm = vtxInr->norm;
			vtxPtr->col = vtxInr->col;
			vtxPtr->uv = vtxInr->uv;
		}
		// ----------------------------------------------------
		// horizontal edges
		// top-edge
		for (Sint32 x = 1; x < edgeLen - 1; x++) 
		{
			const Sint32 y = innerTop;
			const double xFrac = double(x - 1) * frac;
			const double yFrac = double(y - 1) * frac;
			const vector3d p((GetSurfacePointCyl(xFrac, yFrac, m_halfLen) * maxhScale) - clipCentroid);

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
			const Sint32 y = innerBottom;
			const double xFrac = double(x - 1) * frac;
			const double yFrac = double(y - 1) * frac;
			const vector3d p((GetSurfacePointCyl(xFrac, yFrac, m_halfLen) * maxhScale) - clipCentroid);

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

		normals.reset();
		//heights.reset();
		colors.reset();
	}
}

void GeoPlate::ReceiveHeightmaps(SQuadPlateResult *psr)
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
		for (int i=0; i<4; i++)
		{
			assert(!kids[i]);
			const SQuadPlateResult::SSplitResultData& data = psr->data(i);
			assert(i==data.patchID.GetPatchIdx(nD));
			assert(0==data.patchID.GetPatchIdx(nD+1));
			kids[i].reset(new GeoPlate(ctx, geoRing, data.halfLen, data.vbe0, data.vbe1, data.ang0, data.ang1, data.yoffset, nD, data.patchID));
		}

		for (int i=0; i<4; i++)
		{
			const SQuadPlateResult::SSplitResultData& data = psr->data(i);
			kids[i]->heights.reset(data.heights);
			kids[i]->normals.reset(data.normals);
			kids[i]->colors.reset(data.colors);
		}
		for (int i=0; i<4; i++) {
			kids[i]->NeedToUpdateVBOs();
		}
		mHasJobRequest = false;
	}
}

void GeoPlate::ReceiveHeightmap(const SSinglePlateResult *psr)
{
	PROFILE_SCOPED()
	assert(nullptr != psr);
	assert(mHasJobRequest);
	{
		const SSinglePlateResult::SSplitResultData& data = psr->data();
		heights.reset(data.heights);
		normals.reset(data.normals);
		colors.reset(data.colors);
#if DUMP_TO_TEXTURE
		char filename[1024];
		snprintf(filename, 1024, "%s%lld.png", "Vavatch", data.patchID.GetPlateIdx());
		textureDump(filename, ctx->GetEdgeLen(), ctx->GetEdgeLen(), 3, data.colors);
#endif
	}
	mHasJobRequest = false;
}

void GeoPlate::RequestSinglePatch()
{
	if( !heights ) {
		assert(!mHasJobRequest);
		mHasJobRequest = true;

		SSinglePlateRequest *ssrd = new SSinglePlateRequest(ang[0], ang[1], m_yoffset, m_halfLen, vbe[0], vbe[1], 
			m_depth, geoRing->GetSystemBody()->GetPath(), mPatchID, 
			ctx->GetEdgeLen()-2, ctx->GetFrac(), geoRing->GetTerrain());
		m_job = Pi::GetAsyncJobQueue()->Queue(new SinglePlateJob(ssrd));
	}
}
	
void GeoPlate::Render(Graphics::Renderer *renderer, const vector3d &campos, const matrix4x4d &modelView, const Graphics::Frustum &frustum) 
{
	PROFILE_SCOPED()
	// must update the VBOs to calculate the clipRadius...
	UpdateVBOs();

	// Do frustum culling
	if (!frustum.TestPoint(clipCentroid, clipRadius))
		return; // nothing below this patch is visible

	if (kids[0]) 
	{
		for (int i=0; i<4; i++) {
			kids[i]->Render(renderer, campos, modelView, frustum);
		}
	} 
	else if(m_vertexBuffer.get()) 
	{
		RefCountedPtr<Graphics::Material> mat = geoRing->GetSurfaceMaterial();
		Graphics::RenderState *rs = geoRing->GetSurfRenderState();

		const vector3d relpos = clipCentroid - campos;
		renderer->SetTransform(modelView * matrix4x4d::Translation(relpos));

		Pi::statSceneTris += ctx->GetNumTris();
		++Pi::statNumPatches;

		// per-patch detail texture scaling value
		geoRing->GetMaterialParameters().patchDepth = m_depth;
		renderer->SetWireFrameMode(s_bUseWireframe);
		renderer->DrawBufferIndexed(m_vertexBuffer.get(), ctx->GetIndexBuffer(), rs, mat.Get());
		renderer->SetWireFrameMode(false);
		renderer->GetStats().AddToStatCount(Graphics::Stats::STAT_PLATES, 1);
	}
}

void GeoPlate::LODUpdate(vector3d &campos) 
{
	// there should be no LOD update when we have active split requests
	if(mHasJobRequest || !s_bCanUpdate)
		return;

	bool canSplit = true;
	bool canMerge = bool(kids[0].get()!=nullptr);

	// always split at first level
	const double centroidDist = (campos - clipCentroid).Length();
	const bool errorSplit = (centroidDist < m_roughLength);
	if( !(canSplit && (m_depth < std::min(GEOPLATE_MAX_DEPTH, geoRing->GetMaxDepth())) && errorSplit) ) {
		canSplit = false;
	}

	if (canSplit) 
	{
		if (!kids[0]) 
		{
			// we can see this patch so submit the jobs!
			assert(!mHasJobRequest);
			mHasJobRequest = true;

			SQuadPlateRequest *ssrd = new SQuadPlateRequest(ang[0], ang[1], m_yoffset, m_halfLen, vbe[0], vbe[1], 
				m_depth, geoRing->GetSystemBody()->GetPath(), mPatchID, 
				ctx->GetEdgeLen()-2, ctx->GetFrac(), geoRing->GetTerrain());

			// add to the GeoSphere to be processed at end of all LODUpdate requests
			geoRing->AddQuadPlateRequest(centroidDist, ssrd, this);
		}
		else 
		{
			for (int i=0; i<4; ++i) 
				kids[i]->LODUpdate(campos);
		}
	} 
	else if (canMerge) 
	{
		for (int i=0; i<4; i++) 
		{
			canMerge &= kids[i]->canBeMerged();
		}
		if (canMerge) 
		{
			for (int i=0; i<4; ++i) 
				kids[i].reset();
		}
	}
}
