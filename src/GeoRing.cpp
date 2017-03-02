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

// tri edge lengths
#define GEOPLATE_SUBDIVIDE_AT_CAMDIST	1.5 //5.0
#define GEOPLATE_MAX_DEPTH	15

static const int GEOPLATE_MAX_EDGELEN = 55;

#define lerp(t, a, b) ( a + t * (b - a) )

// must be odd numbers
static const int detail_edgeLen[5] = {
	7, 15, 25, 35, 55
};

class BaseGeo {
public:
	BaseGeo(const RefCountedPtr<GeoPatchContext> &ctx_, GeoRing *geoRingPtr, const double halfLength, const double startAng, const double endAng, const double yoffset, const int depth)
		: ctx(ctx_), normals(nullptr), colors(nullptr), geoRing(geoRingPtr)
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
	}

	vector3d GetSurfacePointCyl(const double x, const double y, const double halfLength) const {
		double theta = lerp( x, ang[1], ang[0] );
		
		const vector3d topEndEdge(sin(theta), m_yoffset + (halfLength * 0.5), cos(theta));		// vertices at top edge of circle
		const vector3d bottomEndEdge(sin(theta), m_yoffset - (halfLength * 0.5), cos(theta));	// vertices at bottom edge of circle
		
		const vector3d res = lerp( y, bottomEndEdge, topEndEdge );
		return res;
	}

	static inline Uint32 NUMVERTICES()			{ return edgeLen*edgeLen; }
	static inline Uint32 VBO_COUNT_ALL_IDX()	{ return (((edgeLen-1)*(edgeLen-1))*2*3); }
	static inline Uint32 VBO_COUNT_ALL_TRIS()	{ return (((edgeLen-1)*(edgeLen-1))*2); }

protected:
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

	// --------------------------------------
	static Uint32 edgeLen;
	static int numTris;
	static double frac;
};
// --------------------------------------
Uint32 BaseGeo::edgeLen = 0xFFFFFFFF;
int BaseGeo::numTris = -1;
double BaseGeo::frac = 0.0;

// must be an odd number
#define GEOPLATE_HULL_LEN 7

class GeoPlateHull : public BaseGeo {
public:
	static RefCountedPtr<Graphics::IndexBuffer> indices;
	std::unique_ptr<vector3d[]> vertices;

	// params
	// v0, v1 - define points on the line describing the loop of the ring/orbital.
	// depth - 0 is the topmost plate with each depth+1 describing it's depth within the tree.
	GeoPlateHull(const RefCountedPtr<GeoPatchContext> &ctx_, GeoRing *geoRingPtr, const double halfLength, const double startAng, const double endAng, const double yoffset, const int depth) 
		: BaseGeo(ctx_, geoRingPtr, halfLength, startAng, endAng, yoffset, depth) 
	{
		normals.reset(new vector3f[NUMVERTICES()]);
		vertices.reset(new vector3d[NUMVERTICES()]);
		colors.reset(new Color3ub[NUMVERTICES()]);
	}

	virtual ~GeoPlateHull() {
	}

	static void Init() 
	{
		edgeLen = GEOPLATE_HULL_LEN;
		frac = 1.0 / double(edgeLen-1);
		numTris = 2*(edgeLen-1)*(edgeLen-1);

		std::vector<Uint32> pl_short;
		pl_short.resize(VBO_COUNT_ALL_IDX());
		Uint32 *idx = &pl_short[0];
		for (int x=0; x<edgeLen-1; x++) 
		{
			for (int y=0; y<edgeLen-1; y++) 
			{
				idx[0] = x + edgeLen*y;			assert(idx[0] < NUMVERTICES());
				idx[1] = x+1 + edgeLen*y;		assert(idx[1] < NUMVERTICES());
				idx[2] = x + edgeLen*(y+1);		assert(idx[2] < NUMVERTICES());
				idx+=3;

				idx[0] = x+1 + edgeLen*y;		assert(idx[0] < NUMVERTICES());
				idx[1] = x+1 + edgeLen*(y+1);	assert(idx[1] < NUMVERTICES());
				idx[2] = x + edgeLen*(y+1);		assert(idx[2] < NUMVERTICES());
				//assert(idx < indices+VBO_COUNT_ALL_IDX());
				idx+=3;
			}
		}

		VertexCacheOptimizerUInt vco;
		VertexCacheOptimizerUInt::Result res = vco.Optimize(&pl_short[0], VBO_COUNT_ALL_IDX()/3);
		assert(0 == res);
		//create buffer & copy
		indices.Reset(Pi::renderer->CreateIndexBuffer(pl_short.size(), Graphics::BUFFER_USAGE_STATIC));
		Uint32* idxPtr = indices->Map(Graphics::BUFFER_MAP_WRITE);
		for (Uint32 j = 0; j < pl_short.size(); j++) {
			idxPtr[j] = pl_short[j];
		}
		indices->Unmap();
	}

	void NeedToUpdateVBOs() 
	{
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
		vbd.numVertices = NUMVERTICES();
		vbd.usage = Graphics::BUFFER_USAGE_STATIC;
		m_vertexBuffer.reset(Pi::renderer->CreateVertexBuffer(vbd));

		GeoPatchContext::VBOVertex* VBOVtxPtr = m_vertexBuffer->Map<GeoPatchContext::VBOVertex>(Graphics::BUFFER_MAP_WRITE);
		assert(m_vertexBuffer->GetDesc().stride == sizeof(GeoPatchContext::VBOVertex));

		const vector3d *vtx = vertices.get();
		const vector3f *pNorm = normals.get();
		const Color3ub *pColr = colors.get();

		for (Sint32 y = 0; y<edgeLen; y++) {
			for (Sint32 x = 0; x<edgeLen; x++) {
				const double xFrac = double(x) * frac;
				const double yFrac = double(y) * frac;

				GeoPatchContext::VBOVertex* vtxPtr = &VBOVtxPtr[x + (y*edgeLen)];
				vector3d p(*vtx - clipCentroid);
				vtxPtr->pos = vector3f(p);
				clipRadius = std::max(clipRadius, p.Length());
				++vtx;	// next height

				vtxPtr->norm = pNorm->Normalized();
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

		// ----------------------------------------------------
		// end of mapping
		m_vertexBuffer->Unmap();

		normals.reset();
		vertices.reset();
		colors.reset();
	}

	// Generates full-detail vertices, and also non-edge normals and colors
	void GenerateMesh() {
		vector3d *vts = vertices.get();
		double xfrac = 0;
		double zfrac = 0;

		const vector3d topEnd(0.0, m_yoffset + m_halfLen, 0.0);		// vertices at top edge of circle
		const vector3d btmEnd(0.0, m_yoffset - m_halfLen, 0.0);		// vertices at bottom edge of circle

		for (int y=0; y<edgeLen; ++y) {	// across the width
			xfrac = 0;
			const vector3d axisPt = lerp( zfrac, btmEnd, topEnd );// point along z-axis by zfrac amount
			for (int x=0; x<edgeLen; ++x) {	// along the length (circumference)
				vector3d p;
				const double height = 0.005;
				assert(x!=edgeLen);
				assert(y!=edgeLen);
				if( y==0 || y==edgeLen-1 ) {	// points in trough
					p = GetSurfacePointCyl(xfrac, (y==0 ? 0.0 : 1.0), m_halfLen);
				} else { // outer hull edge
					p = GetSurfacePointCyl(xfrac, zfrac, m_halfLen);
				}

				// vector from axis to point-on-surface
				const vector3d cDir = (p - axisPt).Normalized();
				// vertex is moved in direction of point-on-surface FROM point-in-axis by height.
				*(vts++) = p + (cDir * height);

				// remember this -- we will need it later
				xfrac += frac;
				colors[x + y*edgeLen] = Color3ub(128, 128, 255);
			}
			zfrac += frac;
		}
		assert(vts == &vertices[NUMVERTICES()]);
		// Generate normals
		for (int y=0; y<edgeLen-1; ++y) {
			for (int x=0; x<edgeLen-1; ++x) {
				// normal
				vector3d xy = vertices[x + y*edgeLen];
				vector3d x1 = vertices[x+1 + y*edgeLen];
				vector3d y1 = vertices[x + (y+1)*edgeLen];

				vector3d n = (x1-xy).Cross(y1-xy);
				normals[x + y*edgeLen] = vector3f(n.Normalized());
			}
		}
	}

	void Render(Graphics::Renderer *renderer, const vector3d &campos, const matrix4x4d &modelView, const Graphics::Frustum &frustum)
	{
		// Do frustum culling
		if (!frustum.TestPoint(clipCentroid, clipRadius))
			return; // nothing below this patch is visible

		RefCountedPtr<Graphics::Material> mat = geoRing->GetSurfaceMaterial();
		Graphics::RenderState *rs = geoRing->GetSurfRenderState();

		const vector3d relpos = clipCentroid - campos;
		renderer->SetTransform(modelView * matrix4x4d::Translation(relpos));

		Pi::statSceneTris += VBO_COUNT_ALL_TRIS();
		++Pi::statNumPatches;

		// per-patch detail texture scaling value
		geoRing->GetMaterialParameters().patchDepth = 0;

		renderer->DrawBufferIndexed(m_vertexBuffer.get(), indices.Get(), rs, mat.Get());

		renderer->GetStats().AddToStatCount(Graphics::Stats::STAT_PLATES, 1);
	}
};
//static 
RefCountedPtr<Graphics::IndexBuffer> GeoPlateHull::indices;

// must be an odd number
#define GEOPLATE_WALL_LEN 7

class GeoPlateWall : public BaseGeo {
public:
	static RefCountedPtr<Graphics::IndexBuffer> indices;
	std::unique_ptr<vector3d[]> vertices;
	bool bInnerWall;

	// params
	// v0, v1 - define points on the line describing the loop of the ring/orbital.
	// depth - 0 is the topmost plate with each depth+1 describing it's depth within the tree.
	GeoPlateWall(const RefCountedPtr<GeoPatchContext> &ctx_, GeoRing *geoRingPtr, const bool isInnerWall, const double halfLength, const double startAng, const double endAng, const double yoffset, const int depth) 
		: BaseGeo(ctx_, geoRingPtr, halfLength, startAng, endAng, yoffset, depth), bInnerWall(isInnerWall)
	{
		normals.reset(new vector3f[NUMVERTICES()]);
		vertices.reset(new vector3d[NUMVERTICES()]);
		colors.reset(new Color3ub[NUMVERTICES()]);
	}

	virtual ~GeoPlateWall() {
	}

	static void Init() 
	{
		edgeLen = GEOPLATE_WALL_LEN;
		frac = 1.0 / double(edgeLen-1);
		numTris = 2*(edgeLen-1)*(edgeLen-1);

		std::vector<Uint32> pl_short;
		pl_short.resize(VBO_COUNT_ALL_IDX());
		Uint32 *idx = &pl_short[0];
		for (int x=0; x<edgeLen-1; x++) 
		{
			for (int y=0; y<edgeLen-1; y++) 
			{
				idx[0] = x + edgeLen*y;			assert(idx[0] < NUMVERTICES());
				idx[1] = x+1 + edgeLen*y;		assert(idx[1] < NUMVERTICES());
				idx[2] = x + edgeLen*(y+1);		assert(idx[2] < NUMVERTICES());
				idx+=3;

				idx[0] = x+1 + edgeLen*y;		assert(idx[0] < NUMVERTICES());
				idx[1] = x+1 + edgeLen*(y+1);	assert(idx[1] < NUMVERTICES());
				idx[2] = x + edgeLen*(y+1);		assert(idx[2] < NUMVERTICES());
				//assert(idx < indices+VBO_COUNT_ALL_IDX());
				idx+=3;
			}
		}

		VertexCacheOptimizerUInt vco;
		VertexCacheOptimizerUInt::Result res = vco.Optimize(&pl_short[0], VBO_COUNT_ALL_IDX()/3);
		assert(0 == res);
		//create buffer & copy
		indices.Reset(Pi::renderer->CreateIndexBuffer(pl_short.size(), Graphics::BUFFER_USAGE_STATIC));
		Uint32* idxPtr = indices->Map(Graphics::BUFFER_MAP_WRITE);
		for (Uint32 j = 0; j < pl_short.size(); j++) {
			idxPtr[j] = pl_short[j];
		}
		indices->Unmap();
	}

	void NeedToUpdateVBOs() {
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
		vbd.numVertices = (edgeLen*edgeLen);
		vbd.usage = Graphics::BUFFER_USAGE_STATIC;
		m_vertexBuffer.reset(Pi::renderer->CreateVertexBuffer(vbd));

		GeoPatchContext::VBOVertex* VBOVtxPtr = m_vertexBuffer->Map<GeoPatchContext::VBOVertex>(Graphics::BUFFER_MAP_WRITE);
		assert(m_vertexBuffer->GetDesc().stride == sizeof(GeoPatchContext::VBOVertex));

		const vector3d *vtx = vertices.get();
		const vector3f *pNorm = normals.get();
		const Color3ub *pColr = colors.get();

		for (Sint32 y = 0; y<edgeLen; y++) {
			for (Sint32 x = 0; x<edgeLen; x++) {
				const double xFrac = double(x) * frac;
				const double yFrac = double(y) * frac;

				GeoPatchContext::VBOVertex* vtxPtr = &VBOVtxPtr[x + (y*edgeLen)];
				vector3d p(*vtx - clipCentroid);
				vtxPtr->pos = vector3f(p);
				clipRadius = std::max(clipRadius, p.Length());
				++vtx;	// next height

				vtxPtr->norm = pNorm->Normalized();
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

		// ----------------------------------------------------
		// end of mapping
		m_vertexBuffer->Unmap();

		normals.reset();
		vertices.reset();
		colors.reset();
	}

	const double heights[2][GEOPLATE_WALL_LEN] = {
		{
			0.005,	// start, bottom
			-0.01,	// above start, top inner wall
			-0.011,	// peak, top inner wall
			-0.011,	// peak, top outer wall
			-0.01,	// top edge outer wall
			0.0,	// base, outer wall
			0.005	// end==start, bottom
		}, {
			0.005,	// end==start, bottom
			0.0,	// base, outer wall
			-0.01,	// top edge outer wall
			-0.011,	// peak, top outer wall
			-0.011,	// peak, top inner wall
			-0.01,	// above start, top inner wall
			0.005	// start, bottom
		}
	};

	const double zvals[2][GEOPLATE_WALL_LEN] = {
		{
			-0.0,	// start, bottom
			-0.0,	// above start, top inner wall
			-0.025,	// peak, top inner wall
			-0.075,	// peak, top outer wall
			-0.1,	// top edge outer wall
			-0.1,	// base, outer wall
			-0.0		// end==start, bottom
		}, {
			1.0,	// end==start, bottom
			1.1,	// base, outer wall
			1.1,	// top edge outer wall
			1.075,	// peak, top outer wall
			1.025,	// peak, top inner wall
			1.0,	// above start, top inner wall
			1.0	// start, bottom
		}
	};

	// Generates full-detail vertices, and also non-edge normals and colors
	void GenerateMesh() 
	{
		vector3d *vts = vertices.get();
		double xfrac = 0;
		double zfrac = 0;

		const vector3d topEnd(0.0, m_yoffset + m_halfLen, 0.0);		// vertices at top edge of circle
		const vector3d btmEnd(0.0, m_yoffset - m_halfLen, 0.0);		// vertices at bottom edge of circle

		const int PriIdx = bInnerWall ? 0 : 1;

		vector3d p;
		for (int y=0; y<edgeLen; ++y) {	// across the width
			xfrac = 0;
			const double height = heights[PriIdx][y];
			const vector3d axisPt = lerp( zvals[PriIdx][y], btmEnd, topEnd );// point along z-axis by zfrac amount
			for (int x=0; x<edgeLen; ++x) {	// along the length (circumference)
				p = GetSurfacePointCyl(xfrac, zvals[PriIdx][y], m_halfLen);
				
				// vector from axis to point-on-surface
				const vector3d cDir = (p - axisPt).Normalized();
				// vertex is moved in direction of point-on-surface FROM point-in-axis by height.
				*(vts++) = p + (cDir * height);

				xfrac += frac;
				colors[x + y*edgeLen] = Color3ub(128, 128, 255);
			}
		}
		assert(vts == &vertices[NUMVERTICES()]);
		// Generate normals
		for (int y=0; y<edgeLen-1; ++y) {
			for (int x=0; x<edgeLen-1; ++x) {
				// normal
				vector3d xy = vertices[x + y*edgeLen];
				vector3d x1 = vertices[x+1 + y*edgeLen];
				vector3d y1 = vertices[x + (y+1)*edgeLen];

				vector3d n = (x1-xy).Cross(y1-xy);
				normals[x + y*edgeLen] = vector3f(n.Normalized());
			}
		}

		// Generate bottom row normals
		for (int y=0; y<edgeLen-1; ++y) {
			const int x=edgeLen-1;
			// normal
			vector3d xy = vertices[x + y*edgeLen];
			vector3d x1 = vertices[x-1 + y*edgeLen];
			vector3d y1 = vertices[x + (y+1)*edgeLen];

			vector3d n = (xy-x1).Cross(y1-xy);
			normals[x + y*edgeLen] = vector3f(n.Normalized());
		}

		// Generate last col normals
		{
			const int y=edgeLen-1;
			for (int x=1; x<edgeLen; ++x) {
				// normal
				vector3d xy = vertices[x + y*edgeLen];
				vector3d x1 = vertices[x-1 + y*edgeLen];
				vector3d y1 = vertices[x + (y-1)*edgeLen];

				vector3d n = (xy-x1).Cross(xy-y1);
				normals[x + y*edgeLen] = vector3f(n.Normalized());
			}
		}

		// corner normals
		{
			const int y=edgeLen-1;
			const int x=0;
			{
				// normal
				vector3d xy = vertices[x + y*edgeLen];
				vector3d x1 = vertices[x+1 + y*edgeLen];
				vector3d y1 = vertices[x + (y-1)*edgeLen];

				vector3d n = (x1-xy).Cross(xy-y1);
				normals[x + y*edgeLen] = vector3f(n.Normalized());
			}
		}
		{
			const int y=0;
			const int x=edgeLen-1;
			{
				// normal
				vector3d xy = vertices[x + y*edgeLen];
				vector3d x1 = vertices[x-1 + y*edgeLen];
				vector3d y1 = vertices[x + (y+1)*edgeLen];

				vector3d n = (xy-x1).Cross(y1-xy);
				normals[x + y*edgeLen] = vector3f(n.Normalized());
			}
		}
	}

	void Render(Graphics::Renderer *renderer, const vector3d &campos, const matrix4x4d &modelView, const Graphics::Frustum &frustum)
	{
		// Do frustum culling
		if (!frustum.TestPoint(clipCentroid, clipRadius))
			return; // nothing below this patch is visible

		RefCountedPtr<Graphics::Material> mat = geoRing->GetSurfaceMaterial();
		Graphics::RenderState *rs = geoRing->GetSurfRenderState();

		const vector3d relpos = clipCentroid - campos;
		renderer->SetTransform(modelView * matrix4x4d::Translation(relpos));

		Pi::statSceneTris += (edgeLen-1)*(edgeLen-1)*2;
		++Pi::statNumPatches;

		// per-patch detail texture scaling value
		geoRing->GetMaterialParameters().patchDepth = 0;

		renderer->DrawBufferIndexed(m_vertexBuffer.get(), indices.Get(), rs, mat.Get());

		renderer->GetStats().AddToStatCount(Graphics::Stats::STAT_PLATES, 1);
	}
};
//static 
RefCountedPtr<Graphics::IndexBuffer> GeoPlateWall::indices;

static std::vector<GeoRing*> s_allGeoRings;

void GeoRing::_UpdateLODs()
{
	for (size_t i=0; i<m_plates.size(); i++) {
		m_plates[i]->LODUpdate(m_tempCampos);
	}
}

void GeoRing::Init()
{
	s_patchContext.Reset(new GeoPatchContext(detail_edgeLen[Pi::detail.planets > 4 ? 4 : Pi::detail.planets]));
	OnChangeDetailLevel();
}

// static
void GeoRing::UpdateAllGeoRings()
{
	PROFILE_SCOPED()
	for(std::vector<GeoRing*>::iterator i = s_allGeoRings.begin(); i != s_allGeoRings.end(); ++i)
	{
		(*i)->Update();
	}
}

void GeoRing::OnChangeDetailLevel()
{
	for(std::vector<GeoRing*>::iterator i = s_allGeoRings.begin(); i != s_allGeoRings.end(); ++i) {
		// remove the plates
		for (size_t p=0; p<(*i)->m_plates.size(); p++) {
			if ((*i)->m_plates[p]) {
				delete (*i)->m_plates[p];
			}
		}
		(*i)->m_plates.clear();

		// and strip away the hull
		for (size_t p=0; p<(*i)->m_hull.size(); p++) {
			if ((*i)->m_hull[p]) {
				delete (*i)->m_hull[p];
			}
		}
		(*i)->m_hull.clear();

		// and strip away the walls:
		// inner first...
		for (size_t p=0; p<(*i)->m_wallInner.size(); p++) {
			if ((*i)->m_wallInner[p]) {
				delete (*i)->m_wallInner[p];
			}
		}
		(*i)->m_wallInner.clear();
		// ... then outer
		for (size_t p=0; p<(*i)->m_wallOuter.size(); p++) {
			if ((*i)->m_wallOuter[p]) {
				delete (*i)->m_wallOuter[p];
			}
		}
		(*i)->m_wallOuter.clear();
	}

	GeoPlateHull::Init();
	GeoPlateWall::Init();
	for(std::vector<GeoRing*>::iterator i = s_allGeoRings.begin(); i != s_allGeoRings.end(); ++i) {
		(*i)->BuildFirstPatches();
	}
}

GeoRing::GeoRing(const SystemBody *body): m_sbody(body), m_terrain(Terrain::InstanceTerrain(body)), m_maxDepth(0), 
	m_hasTempCampos(false), m_tempCampos(0.0), m_tempFrustum(800, 600, 0.5, 1.0, 1000.0)
{
	s_allGeoRings.push_back(this);

	// a little aside, this calculates the radius(R) of a torus
	// that generates 1g for a given rpm.
	// the numbers get very big very quickly :(
	/*double rpm = (1.0/60.0)/24.0;//1.0;
	double cal = ((M_PI*rpm)/30.0);
	double g = 9.81;
	double R = g / (cal*cal);
	printf("Radius of torus with 1g at %.lf rpm is R = %.lf\n", rpm, R);*/
}

GeoRing::~GeoRing()
{
	s_allGeoRings.erase(std::find(s_allGeoRings.begin(), s_allGeoRings.end(), this));

	for (size_t i=0; i<m_plates.size(); i++) {
		if (m_plates[i]) {
			delete m_plates[i];
		}
		if (m_hull[i]) {
			delete m_hull[i];
		}
		if (m_wallInner[i]) {
			delete m_wallInner[i];
		}
		if (m_wallOuter[i]) {
			delete m_wallOuter[i];
		}
	}
}

void GeoRing::BuildFirstPatches()
{
	// This is a somewhat hacky way to decide how "wide" an orbital should be
	// for a given comparison to the earth.
	// The NUM_SEGMENTS defines how many "plates" we'd use if the orbital was the radius of
	// the planet Earth, at 8 this makes each plate approximately half the width & length
	// of a patch used on the GeoSphere Earth.
	// An SBoday radius is stored in terms of it's multiples for the Earths radius, which
	// gives us the perfect multiplier for our NUM_SEGMENTS. Orbitals smaller than Earth
	// will use less segments whereas larger Orbitals will use more.
	#define NUM_SEGMENTS 8
	const int numSegments = (NUM_SEGMENTS * (m_sbody->GetRadius() / EARTH_RADIUS));

	CalculateMaxPatchDepth();

	std::vector<vector3d>	points;
	std::vector<double>		angles;
    for (int i = 0; i <= numSegments; ++i) {
		double angle = ((M_PI * 360.0) / 180.0) * i / numSegments;
		double sinval = sin(angle);
		double cosval = cos(angle);
		vector3d vp(cosval, sinval, 0.0 );
		points.push_back( vp.Normalized() );
		angles.push_back( angle );
    }

	// calculate the width of the orbital
	mRingWidth = (points[1] - points[0]).Length();

	// build the terrain plates
	m_plates.clear();
	for( size_t i=0 ; i<points.size()-1 ; ++i ) {
		m_plates.push_back( new GeoPlate(s_patchContext, this, mRingWidth, points[i], points[i+1], angles[i], angles[i+1], 0.0, 0, i) );
	}

	// request meshes
	for (size_t i=0; i<m_plates.size(); i++) {
		m_plates[i]->RequestSinglePatch();
	}

	// create the outer hull
	m_hull.clear();
	for( size_t i=0 ; i<points.size()-1 ; ++i ) {
		double len = (points[i+1] - points[i]).Length();
		m_hull.push_back( new GeoPlateHull(s_patchContext, this, len, angles[i+1], angles[i], 0.0, 0) );
	}
	for (size_t i=0; i<m_hull.size(); i++) m_hull[i]->GenerateMesh();
	for (size_t i=0; i<m_hull.size(); i++) m_hull[i]->NeedToUpdateVBOs();

	// create the inner wall
	m_wallInner.clear();
	for( size_t i=0 ; i<points.size()-1 ; ++i ) {
		double len = (points[i+1] - points[i]).Length();
		m_wallInner.push_back( new GeoPlateWall(s_patchContext, this, true, len, angles[i+1], angles[i], 0.0, 0) );
	}
	for (size_t i=0; i<m_wallInner.size(); i++) m_wallInner[i]->GenerateMesh();
	for (size_t i=0; i<m_wallInner.size(); i++) m_wallInner[i]->NeedToUpdateVBOs();

	// create the outer wall
	m_wallOuter.clear();
	for( size_t i=0 ; i<points.size()-1 ; ++i ) {
		double len = (points[i+1] - points[i]).Length();
		m_wallOuter.push_back( new GeoPlateWall(s_patchContext, this, false, len, angles[i+1], angles[i], 0.0, 0) );
	}
	for (size_t i=0; i<m_wallOuter.size(); i++) m_wallOuter[i]->GenerateMesh();
	for (size_t i=0; i<m_wallOuter.size(); i++) m_wallOuter[i]->NeedToUpdateVBOs();

	m_initStage = eRequestedFirstPatches;
}

static const double gs_targetPatchTriLength(100.0);
void GeoRing::CalculateMaxPatchDepth()
{
	m_maxDepth = 0;
	const double circumference = 2.0 * M_PI * m_sbody->GetRadius();
	// calculate length of each edge segment (quad) times 4 due to that being the number around the sphere (1 per side, 4 sides for Root).
	double edgeMetres = circumference / double((s_patchContext->GetEdgeLen()-2) * 8);
	// find out what depth we reach the desired resolution
	while (edgeMetres>gs_targetPatchTriLength && m_maxDepth<GEOPLATE_MAX_DEPTH) {
		edgeMetres *= 0.5;
		++m_maxDepth;
	}
	if(m_maxDepth==0)
		m_maxDepth = GEOPLATE_MAX_DEPTH;
}

void GeoRing::Update()
{
	if (m_plates.empty()) {
		BuildFirstPatches();
		return;
	}
	_UpdateLODs();
	switch(m_initStage)
	{
	case eBuildFirstPatches:
		BuildFirstPatches();
		break;
	case eRequestedFirstPatches:
		{
			ProcessSplitResults();
			uint8_t numValidPatches = 0;
			for (size_t i=0; i<m_plates.size(); i++) {
				if(m_plates[i]->HasData()) {
					++numValidPatches;
				}
			}
			m_initStage = (m_plates.size()==numValidPatches) ? eReceivedFirstPatches : eRequestedFirstPatches;
		} break;
	case eReceivedFirstPatches:
		{
			for (size_t i=0; i<m_plates.size(); i++) {
				m_plates[i]->NeedToUpdateVBOs();
			}
			m_initStage = eDefaultUpdateState;
		} break;
	case eDefaultUpdateState:
		if(m_hasTempCampos) {
			ProcessSplitResults();
			for (size_t i=0; i<m_plates.size(); i++) {
				m_plates[i]->LODUpdate(m_tempCampos);//, m_tempFrustum);
			}
			ProcessQuadPlateRequests();
		}
		break;
	}
}

static void DrawAtmosphereSurface(const vector3d &campos, float rad)
{
	/*const int LAT_SEGS = 20;
	const int LONG_SEGS = 20;
	vector3d yaxis = campos.Normalized();
	vector3d zaxis = vector3d(1.0,0.0,0.0).Cross(yaxis).Normalized();
	vector3d xaxis = yaxis.Cross(zaxis);
	const matrix4x4d m = matrix4x4d::MakeRotMatrix(xaxis, yaxis, zaxis).InverseOf();

	glPushMatrix();
	glScalef(rad, rad, rad);
	glMultMatrixd(&m[0]);

	// what is this? Well, angle to the horizon is:
	// acos(planetRadius/viewerDistFromSphereCentre)
	// and angle from this tangent on to atmosphere is:
	// acos(planetRadius/atmosphereRadius) ie acos(1.0/1.01244blah)
	double endAng = acos(1.0/campos.Length())+acos(1.0/rad);
	double latDiff = endAng / double(LAT_SEGS);

	double rot = 0.0;
	float sinCosTable[LONG_SEGS+1][2];
	for (int i=0; i<=LONG_SEGS; i++, rot += 2.0*M_PI/double(LONG_SEGS)) {
		sinCosTable[i][0] = float(sin(rot));
		sinCosTable[i][1] = float(cos(rot));
	}

	// Tri-fan above viewer
	glBegin(GL_TRIANGLE_FAN);
	glVertex3f(0.0f, 1.0f, 0.0f);
	for (int i=0; i<=LONG_SEGS; i++) {
		glVertex3f(sin(latDiff)*sinCosTable[i][0], cos(latDiff), -sin(latDiff)*sinCosTable[i][1]);
	}
	glEnd();

	// and wound latitudinal strips
	double lat = latDiff;
	for (int j=1; j<LAT_SEGS; j++, lat += latDiff) {
		glBegin(GL_TRIANGLE_STRIP);
		float cosLat = cos(lat);
		float sinLat = sin(lat);
		float cosLat2 = cos(lat+latDiff);
		float sinLat2 = sin(lat+latDiff);
		for (int i=0; i<=LONG_SEGS; i++) {
			glVertex3f(sinLat*sinCosTable[i][0], cosLat, -sinLat*sinCosTable[i][1]);
			glVertex3f(sinLat2*sinCosTable[i][0], cosLat2, -sinLat2*sinCosTable[i][1]);
		}
		glEnd();
	}

	glPopMatrix();*/
}

void GeoRing::Render(Graphics::Renderer *renderer, const matrix4x4d &modelView, vector3d campos, const float radius, const std::vector<Camera::Shadow> &shadows) 
{
	PROFILE_SCOPED()
	// store this for later usage in the update method.
	m_tempCampos = campos;
	m_hasTempCampos = true;

	matrix4x4d trans = modelView;
	trans.Translate(-campos.x, -campos.y, -campos.z);
	renderer->SetTransform(trans); //need to set this for the following line to work
	matrix4x4d modv;
	matrix4x4d proj;
	matrix4x4ftod(renderer->GetCurrentModelView(), modv);
	matrix4x4ftod(renderer->GetCurrentProjection(), proj);
	Graphics::Frustum frustum( modv, proj );
	m_tempFrustum = frustum;
	const float atmosRadius = ATMOSPHERE_RADIUS;

	//First draw - create materials (they do not change afterwards)
	if (!m_surfaceMaterial)
		SetUpMaterials();
	
	// no frustum test of entire geoRing, since Space::Render does this
	// for each body using its GetBoundingRadius() value
	{
		//Update material parameters
		//XXX no need to calculate AP every frame
		m_materialParameters.atmosphere = m_sbody->CalcAtmosphereParams();
		m_materialParameters.atmosphere.center = trans * vector3d(0.0, 0.0, 0.0);
		m_materialParameters.atmosphere.planetRadius = radius;

		m_materialParameters.shadows = shadows;

		m_materialParameters.maxPatchDepth = GetMaxDepth();

		m_surfaceMaterial->specialParameter0 = &m_materialParameters;

		//if (m_materialParameters.atmosphere.atmosDensity > 0.0) 
		//{
		//	m_atmosphereMaterial->specialParameter0 = &m_materialParameters;
		//
		//	// make atmosphere sphere slightly bigger than required so
		//	// that the edges of the pixel shader atmosphere jizz doesn't
		//	// show ugly polygonal angles
		//	DrawAtmosphereSurface(renderer, trans, campos,
		//		m_materialParameters.atmosphere.atmosRadius*1.01,
		//		m_atmosRenderState, m_atmosphereMaterial);
		//}
	}

	Color ambient;
	Color &emission = m_surfaceMaterial->emissive;

	// save old global ambient
	const Color oldAmbient = renderer->GetAmbientColor();

	// give planet some ambient lighting if the viewer is close to it
	double camdist = campos.Length();
	camdist = 0.1 / (camdist*camdist);
	// why the fuck is this returning 0.1 when we are sat on the planet??
	// JJ: Because campos is relative to a unit-radius planet - 1.0 at the surface
	// XXX oh well, it is the value we want anyway...
	ambient.r = ambient.g = ambient.b = camdist * 255;
	ambient.a = 255;
	
	renderer->SetAmbientColor(ambient);

	renderer->SetTransform(modelView);
	
	for (size_t i=0; i<m_hull.size(); ++i) {
		m_hull[i]->Render(renderer, campos, modelView, frustum);
	}
	for (size_t i=0; i<m_wallInner.size(); ++i) {
		m_wallInner[i]->Render(renderer, campos, modelView, frustum);
	}
	for (size_t i=0; i<m_wallOuter.size(); ++i) {
		m_wallOuter[i]->Render(renderer, campos, modelView, frustum);
	}

#if 1
	for (size_t i=0; i<m_plates.size(); ++i) {
		m_plates[i]->Render(renderer, campos, modelView, frustum);
	}
#else
	for (size_t i=0; i<m_plates.size(); i+=2) {
		m_plates[i]->Render(renderer, campos, modelView, frustum);
	}
#endif
	renderer->SetAmbientColor(oldAmbient);

	renderer->GetStats().AddToStatCount(Graphics::Stats::STAT_ORBITALS, 1);
}

void GeoRing::SetUpMaterials()
{
	//solid
	Graphics::RenderStateDesc rsd;
	m_surfRenderState = Pi::renderer->CreateRenderState(rsd);

	// Request material for this star or planet, with or without
	// atmosphere. Separate material for surface and sky.
	Graphics::MaterialDescriptor surfDesc;
	//const Uint32 effect_flags = m_terrain->GetSurfaceEffects();
	//if (effect_flags & Terrain::EFFECT_LAVA)
	//	surfDesc.effect = Graphics::EFFECT_GEOSPHERE_TERRAIN_WITH_LAVA;
	//else if (effect_flags & Terrain::EFFECT_WATER)
	//	surfDesc.effect = Graphics::EFFECT_GEOSPHERE_TERRAIN_WITH_WATER;
	//else
		surfDesc.effect = Graphics::EFFECT_GEORING_TERRAIN;

	{
		//planetoid with or without atmosphere
		const SystemBody::AtmosphereParameters ap(m_sbody->CalcAtmosphereParams());
		surfDesc.lighting = true;
		if(ap.atmosDensity > 0.0) {
			surfDesc.quality |= Graphics::HAS_ATMOSPHERE;
		} else {
			surfDesc.quality &= ~Graphics::HAS_ATMOSPHERE;
		}
	}

	surfDesc.quality |= Graphics::HAS_ECLIPSES;
	const bool bEnableDetailMaps = (Pi::config->Int("DisableDetailMaps") == 0);
	if (bEnableDetailMaps) {
		surfDesc.quality |= Graphics::HAS_DETAIL_MAPS;
	}
	m_surfaceMaterial.Reset(Pi::renderer->CreateMaterial(surfDesc));

	m_texHi.Reset( Graphics::TextureBuilder::Model("textures/high.dds").GetOrCreateTexture(Pi::renderer, "model") );
	m_texLo.Reset( Graphics::TextureBuilder::Model("textures/low.dds").GetOrCreateTexture(Pi::renderer, "model") );
	m_surfaceMaterial->texture0 = m_texHi.Get();
	m_surfaceMaterial->texture1 = m_texLo.Get();
}

//static
bool GeoRing::OnAddQuadPlateResult(const SystemPath &path, SQuadPlateResult *res)
{
	// Find the correct GeoSphere via it's system path, and give it the split result
	for(std::vector<GeoRing*>::iterator i=s_allGeoRings.begin(), iEnd=s_allGeoRings.end(); i!=iEnd; ++i) {
		if( path == (*i)->GetSystemBody()->GetPath() ) {
			(*i)->AddQuadPlateResult(res);
			return true;
		}
	}
	// GeoSphere not found to return the data to, cancel and delete it instead
	if( res ) {
		res->OnCancel();
		delete res;
	}
	return false;
}

//static
bool GeoRing::OnAddSinglePlateResult(const SystemPath &path, SSinglePlateResult *res)
{
	// Find the correct GeoSphere via it's system path, and give it the split result
	for(std::vector<GeoRing*>::iterator i=s_allGeoRings.begin(), iEnd=s_allGeoRings.end(); i!=iEnd; ++i) {
		if( path == (*i)->GetSystemBody()->GetPath() ) {
			(*i)->AddSinglePlateResult(res);
			return true;
		}
	}
	// GeoSphere not found to return the data to, cancel and delete it instead
	if( res ) {
		res->OnCancel();
		delete res;
	}
	return false;
}

bool GeoRing::AddQuadPlateResult(SQuadPlateResult *res)
{
	bool result = false;
	assert(res);
	assert(mQuadPlateResults.size()<MAX_SPLIT_OPERATIONS);
	if(mQuadPlateResults.size()<MAX_SPLIT_OPERATIONS) {
		mQuadPlateResults.push_back(res);
		result = true;
	}
	return result;
}

bool GeoRing::AddSinglePlateResult(SSinglePlateResult *res)
{
	bool result = false;
	assert(res);
	assert(mSinglePlateResults.size()<MAX_SPLIT_OPERATIONS);
	if(mSinglePlateResults.size()<MAX_SPLIT_OPERATIONS) {
		mSinglePlateResults.push_back(res);
		result = true;
	}
	return result;
}

void GeoRing::ProcessSplitResults()
{
	// now handle the single split results that define the base level of the quad tree
	{
		std::deque<SSinglePlateResult*>::iterator iter = mSinglePlateResults.begin();
		while(iter!=mSinglePlateResults.end())
		{
			// finally pass SplitResults
			SSinglePlateResult *psr = (*iter);
			assert(psr);

			const uint64_t plateIdx = psr->plate();
			if( m_plates[plateIdx] ) {
				m_plates[plateIdx]->ReceiveHeightmap(psr);
			} else {
				psr->OnCancel();
			}

			// tidyup
			delete psr;

			// Next!
			++iter;
		}
		mSinglePlateResults.clear();
	}

	// now handle the quad split results
	{
		std::deque<SQuadPlateResult*>::iterator iter = mQuadPlateResults.begin();
		while(iter!=mQuadPlateResults.end())
		{
			// finally pass SplitResults
			SQuadPlateResult *psr = (*iter);
			assert(psr);

			const uint64_t plateIdx = psr->plate();
			if( m_plates[plateIdx] ) {
				m_plates[plateIdx]->ReceiveHeightmaps(psr);
			} else {
				psr->OnCancel();
			}

			// tidyup
			delete psr;

			// Next!
			++iter;
		}
		mQuadPlateResults.clear();
	}
}

void GeoRing::AddQuadPlateRequest(double dist, SQuadPlateRequest *pReq, GeoPlate *pPlate)
{
	mQuadPlateRequests.push_back(TDistanceRequest(dist, pReq, pPlate));
}

void GeoRing::ProcessQuadPlateRequests()
{
	class RequestDistanceSort {
	public:
		bool operator()(const TDistanceRequest &a, const TDistanceRequest &b)
		{
			return a.mDistance < b.mDistance;
		}
	};
	std::sort(mQuadPlateRequests.begin(), mQuadPlateRequests.end(), RequestDistanceSort());

	for(auto iter : mQuadPlateRequests) {
		SQuadPlateRequest *ssrd = iter.mpRequest;
		iter.mpRequester->ReceiveJobHandle(Pi::GetAsyncJobQueue()->Queue(new QuadPlateJob(ssrd)));
	}
	mQuadPlateRequests.clear();
}

RefCountedPtr<GeoPatchContext> GeoRing::s_patchContext;
