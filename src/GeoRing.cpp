#include "libs.h"
#include "GeoRing.h"
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

//#define GEOPLATE_USE_THREADING

// tri edge lengths
#define GEOPLATE_SUBDIVIDE_AT_CAMDIST	1.5 //5.0
#define GEOPLATE_MAX_DEPTH	15
// must be an odd number
#define GEOPLATE_EDGELEN_DEFAULT	7//15
#define GEOPLATE_NUMVERTICES	(GEOPLATE_EDGELEN*GEOPLATE_EDGELEN)
#define VBO_COUNT_ALL_IDX  (((GEOPLATE_EDGELEN-1)*(GEOPLATE_EDGELEN-1))*2*3)
#define VBO_COUNT_ALL_TRIS  (((GEOPLATE_EDGELEN-1)*(GEOPLATE_EDGELEN-1))*2)

int GEOPLATE_EDGELEN = GEOPLATE_EDGELEN_DEFAULT;	// 15;
static const int GEOPLATE_MAX_EDGELEN = 55;
static double GEOPLATE_FRAC;
static double GEOPLATEHULL_FRAC;
static double GEOPLATEWALL_FRAC;

#define lerp(t, a, b) ( a + t * (b - a) )

#define PRINT_VECTOR(_v) printf("%f,%f,%f\n", (_v).x, (_v).y, (_v).z);

#define SAFE_DELETE(d) delete d; d=0;
#define SAFE_DELETE_ARRAY(d) delete[] d; d=0;

#pragma pack(4)
struct VBOVertex
{
	vector3f pos;
	vector3f norm;
	Color4ub col;
	vector2f uv;
};
#pragma pack()

#define VBO_COUNT_HI_EDGE  (3*(GEOPLATE_EDGELEN-1))
#define VBO_COUNT_MID_IDX  (4*3*(GEOPLATE_EDGELEN-3) + 2*(GEOPLATE_EDGELEN-3)*(GEOPLATE_EDGELEN-3)*3)
//                          ^^ serrated teeth bit      ^^^ square inner bit

// hold the 16 possible terrain edge connections
typedef struct {
	std::vector<unsigned short> v;
} VecShortHolder;

#define NUM_VBE 2

class BaseGeo {
public:
	BaseGeo(const double halfLength, const double startAng, const double endAng, const double yoffset, const int depth)
		: vertices(NULL), normals(NULL), colors(NULL), geoRing(NULL)
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

	void SetGeoRingPtr( GeoRing *ptr ) { assert(ptr!=NULL); geoRing = ptr; }

protected:
	double ang[NUM_VBE];
	double m_halfLen;
	double m_yoffset;		// offset from the centre i.e. newyoffset = (m_yoffset + m_halfLen);
	double *heights;
	vector3d *vertices;
	vector3d *normals;
	Color3ub *colors;
	GeoRing *geoRing;
	std::unique_ptr<Graphics::VertexBuffer> m_vertexBuffer;
	double m_roughLength;
	vector3d clipCentroid;
	double clipRadius;
};

class GeoPlateHull : public BaseGeo {
public:
	static RefCountedPtr<Graphics::IndexBuffer> indices;

	// params
	// v0, v1 - define points on the line describing the loop of the ring/orbital.
	// depth - 0 is the topmost plate with each depth+1 describing it's depth within the tree.
	GeoPlateHull(const double halfLength, const double startAng, const double endAng, const double yoffset, const int depth) 
		: BaseGeo(halfLength, startAng, endAng, yoffset, depth) 
	{
		normals = new vector3d[GEOPLATE_NUMVERTICES];
		vertices = new vector3d[GEOPLATE_NUMVERTICES];
		colors = new Color3ub[GEOPLATE_NUMVERTICES];
	}

	~GeoPlateHull() {
		delete vertices;
		delete normals;
		delete colors;
	}

	static void Init() 
	{
		GEOPLATEHULL_FRAC = 1.0 / double(GEOPLATE_EDGELEN-1);

		std::vector<Uint32> pl_short;
		pl_short.resize(VBO_COUNT_ALL_IDX);
		Uint32 *idx = &pl_short[0];
		for (int x=0; x<GEOPLATE_EDGELEN-1; x++) 
		{
			for (int y=0; y<GEOPLATE_EDGELEN-1; y++) 
			{
				idx[0] = x + GEOPLATE_EDGELEN*y;		assert(idx[0] < GEOPLATE_NUMVERTICES);
				idx[1] = x+1 + GEOPLATE_EDGELEN*y;		assert(idx[1] < GEOPLATE_NUMVERTICES);
				idx[2] = x + GEOPLATE_EDGELEN*(y+1);	assert(idx[2] < GEOPLATE_NUMVERTICES);
				idx+=3;

				idx[0] = x+1 + GEOPLATE_EDGELEN*y;		assert(idx[0] < GEOPLATE_NUMVERTICES);
				idx[1] = x+1 + GEOPLATE_EDGELEN*(y+1);	assert(idx[1] < GEOPLATE_NUMVERTICES);
				idx[2] = x + GEOPLATE_EDGELEN*(y+1);	assert(idx[2] < GEOPLATE_NUMVERTICES);
				//assert(idx < indices+VBO_COUNT_ALL_IDX);
				idx+=3;
			}
		}

		VertexCacheOptimizerUInt vco;
		VertexCacheOptimizerUInt::Result res = vco.Optimize(&pl_short[0], VBO_COUNT_ALL_IDX/3);
		assert(0 == res);
		//create buffer & copy
		indices.Reset(Pi::renderer->CreateIndexBuffer(pl_short.size(), Graphics::BUFFER_USAGE_STATIC));
		Uint32* idxPtr = indices->Map(Graphics::BUFFER_MAP_WRITE);
		for (Uint32 j = 0; j < pl_short.size(); j++) {
			idxPtr[j] = pl_short[j];
		}
		indices->Unmap();
	}

	void UpdateVBOs() 
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
		vbd.numVertices = GEOPLATE_NUMVERTICES;
		vbd.usage = Graphics::BUFFER_USAGE_STATIC;
		m_vertexBuffer.reset(Pi::renderer->CreateVertexBuffer(vbd));

		VBOVertex* VBOVtxPtr = m_vertexBuffer->Map<VBOVertex>(Graphics::BUFFER_MAP_WRITE);
		assert(m_vertexBuffer->GetDesc().stride == sizeof(VBOVertex));

		const vector3d *vtx = vertices;
		const vector3d *pNorm = normals;
		const Color3ub *pColr = colors;

		for (Sint32 y = 0; y<GEOPLATE_EDGELEN; y++) {
			for (Sint32 x = 0; x<GEOPLATE_EDGELEN; x++) {
				const double xFrac = double(x) * GEOPLATEHULL_FRAC;
				const double yFrac = double(y) * GEOPLATEHULL_FRAC;

				VBOVertex* vtxPtr = &VBOVtxPtr[x + (y*GEOPLATE_EDGELEN)];
				vector3d p(*vtx - clipCentroid);
				vtxPtr->pos = vector3f(p);
				clipRadius = std::max(clipRadius, p.Length());
				++vtx;	// next height

				vtxPtr->norm = vector3f(pNorm->Normalized());
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

		// ----------------------------------------------------
		// end of mapping
		m_vertexBuffer->Unmap();
	}

	// Generates full-detail vertices, and also non-edge normals and colors
	void GenerateMesh() {
		vector3d *vts = vertices;
		double xfrac = 0;
		double zfrac = 0;
#ifdef _DEBUG
		memset(normals, 0, sizeof(vector3d)*GEOPLATE_NUMVERTICES);
#endif
		const vector3d topEnd(0.0, m_yoffset + m_halfLen, 0.0);		// vertices at top edge of circle
		const vector3d btmEnd(0.0, m_yoffset - m_halfLen, 0.0);		// vertices at bottom edge of circle

		for (int y=0; y<GEOPLATE_EDGELEN; ++y) {	// across the width
			xfrac = 0;
			const vector3d axisPt = lerp( zfrac, btmEnd, topEnd );// point along z-axis by zfrac amount
			for (int x=0; x<GEOPLATE_EDGELEN; ++x) {	// along the length (circumference)
				vector3d p;
				const double height = 0.005;
				assert(x!=GEOPLATE_EDGELEN);
				assert(y!=GEOPLATE_EDGELEN);
				if( y==0 || y==GEOPLATE_EDGELEN-1 ) {	// points in trough
					p = GetSurfacePointCyl(xfrac, (y==0 ? 0.0 : 1.0), m_halfLen);
				} else { // outer hull edge
					p = GetSurfacePointCyl(xfrac, zfrac, m_halfLen);
				}

				// vector from axis to point-on-surface
				const vector3d cDir = (p - axisPt).Normalized();
				// vertex is moved in direction of point-on-surface FROM point-in-axis by height.
				*(vts++) = p + (cDir * height);

				// remember this -- we will need it later
				xfrac += GEOPLATEHULL_FRAC;
				colors[x + y*GEOPLATE_EDGELEN] = Color3ub(128, 128, 255);
			}
			zfrac += GEOPLATEHULL_FRAC;
		}
		assert(vts == &vertices[GEOPLATE_NUMVERTICES]);
		// Generate normals
		for (int y=0; y<GEOPLATE_EDGELEN-1; ++y) {
			for (int x=0; x<GEOPLATE_EDGELEN-1; ++x) {
				// normal
				vector3d xy = vertices[x + y*GEOPLATE_EDGELEN];
				vector3d x1 = vertices[x+1 + y*GEOPLATE_EDGELEN];
				vector3d y1 = vertices[x + (y+1)*GEOPLATE_EDGELEN];

				vector3d n = (x1-xy).Cross(y1-xy);
				normals[x + y*GEOPLATE_EDGELEN] = n.Normalized();
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

		Pi::statSceneTris += (VBO_COUNT_ALL_TRIS);
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
#define GEOPLATE_NUM_WALL_VERTICES	(GEOPLATE_WALL_LEN*GEOPLATE_WALL_LEN)

class GeoPlateWall : public BaseGeo {
public:
	static RefCountedPtr<Graphics::IndexBuffer> indices;
	bool bInnerWall;

	// params
	// v0, v1 - define points on the line describing the loop of the ring/orbital.
	// depth - 0 is the topmost plate with each depth+1 describing it's depth within the tree.
	GeoPlateWall(const bool isInnerWall, const double halfLength, const double startAng, const double endAng, const double yoffset, const int depth) 
		: BaseGeo(halfLength, startAng, endAng, yoffset, depth), bInnerWall(isInnerWall)
	{
		normals = new vector3d[GEOPLATE_NUM_WALL_VERTICES];
		vertices = new vector3d[GEOPLATE_NUM_WALL_VERTICES];
		colors = new Color3ub[GEOPLATE_NUM_WALL_VERTICES];
	}

	~GeoPlateWall() {
		delete vertices;
		delete normals;
		delete colors;
	}

	static void Init() {
		GEOPLATEWALL_FRAC = 1.0 / double(GEOPLATE_WALL_LEN-1);

		std::vector<Uint32> pl_short;
		pl_short.resize(VBO_COUNT_ALL_IDX);
		Uint32 *idx = &pl_short[0];
		for (int x=0; x<GEOPLATE_WALL_LEN-1; x++) 
		{
			for (int y=0; y<GEOPLATE_WALL_LEN-1; y++) 
			{
				idx[0] = x + GEOPLATE_WALL_LEN*y;		assert(idx[0] < GEOPLATE_NUMVERTICES);
				idx[1] = x+1 + GEOPLATE_WALL_LEN*y;		assert(idx[1] < GEOPLATE_NUMVERTICES);
				idx[2] = x + GEOPLATE_WALL_LEN*(y+1);	assert(idx[2] < GEOPLATE_NUMVERTICES);
				idx+=3;

				idx[0] = x+1 + GEOPLATE_WALL_LEN*y;		assert(idx[0] < GEOPLATE_NUMVERTICES);
				idx[1] = x+1 + GEOPLATE_WALL_LEN*(y+1);	assert(idx[1] < GEOPLATE_NUMVERTICES);
				idx[2] = x + GEOPLATE_WALL_LEN*(y+1);	assert(idx[2] < GEOPLATE_NUMVERTICES);
				//assert(idx < indices+VBO_COUNT_ALL_IDX);
				idx+=3;
			}
		}

		VertexCacheOptimizerUInt vco;
		VertexCacheOptimizerUInt::Result res = vco.Optimize(&pl_short[0], VBO_COUNT_ALL_IDX/3);
		assert(0 == res);
		//create buffer & copy
		indices.Reset(Pi::renderer->CreateIndexBuffer(pl_short.size(), Graphics::BUFFER_USAGE_STATIC));
		Uint32* idxPtr = indices->Map(Graphics::BUFFER_MAP_WRITE);
		for (Uint32 j = 0; j < pl_short.size(); j++) {
			idxPtr[j] = pl_short[j];
		}
		indices->Unmap();
	}

	void UpdateVBOs() {
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
		vbd.numVertices = (GEOPLATE_WALL_LEN*GEOPLATE_WALL_LEN);
		vbd.usage = Graphics::BUFFER_USAGE_STATIC;
		m_vertexBuffer.reset(Pi::renderer->CreateVertexBuffer(vbd));

		VBOVertex* VBOVtxPtr = m_vertexBuffer->Map<VBOVertex>(Graphics::BUFFER_MAP_WRITE);
		assert(m_vertexBuffer->GetDesc().stride == sizeof(VBOVertex));

		const vector3d *vtx = vertices;
		const vector3d *pNorm = normals;
		const Color3ub *pColr = colors;

		for (Sint32 y = 0; y<GEOPLATE_WALL_LEN; y++) {
			for (Sint32 x = 0; x<GEOPLATE_WALL_LEN; x++) {
				const double xFrac = double(x) * GEOPLATEHULL_FRAC;
				const double yFrac = double(y) * GEOPLATEHULL_FRAC;

				VBOVertex* vtxPtr = &VBOVtxPtr[x + (y*GEOPLATE_WALL_LEN)];
				vector3d p(*vtx - clipCentroid);
				vtxPtr->pos = vector3f(p);
				clipRadius = std::max(clipRadius, p.Length());
				++vtx;	// next height

				vtxPtr->norm = vector3f(pNorm->Normalized());
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

		// ----------------------------------------------------
		// end of mapping
		m_vertexBuffer->Unmap();
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
	void GenerateMesh() {
		vector3d *vts = vertices;
		double xfrac = 0;
		double zfrac = 0;
#ifdef _DEBUG
		memset(normals, 0, sizeof(vector3d)*GEOPLATE_NUM_WALL_VERTICES);
#endif
		const vector3d topEnd(0.0, m_yoffset + m_halfLen, 0.0);		// vertices at top edge of circle
		const vector3d btmEnd(0.0, m_yoffset - m_halfLen, 0.0);		// vertices at bottom edge of circle

		const int PriIdx = bInnerWall ? 0 : 1;

		vector3d p;
		for (int y=0; y<GEOPLATE_WALL_LEN; ++y) {	// across the width
			xfrac = 0;
			const double height = heights[PriIdx][y];
			const vector3d axisPt = lerp( zvals[PriIdx][y], btmEnd, topEnd );// point along z-axis by zfrac amount
			for (int x=0; x<GEOPLATE_WALL_LEN; ++x) {	// along the length (circumference)
				p = GetSurfacePointCyl(xfrac, zvals[PriIdx][y], m_halfLen);
				
				// vector from axis to point-on-surface
				const vector3d cDir = (p - axisPt).Normalized();
				// vertex is moved in direction of point-on-surface FROM point-in-axis by height.
				*(vts++) = p + (cDir * height);

				xfrac += GEOPLATEWALL_FRAC;
				colors[x + y*GEOPLATE_WALL_LEN] = Color3ub(128, 128, 255);
			}
		}
		assert(vts == &vertices[GEOPLATE_NUM_WALL_VERTICES]);
		// Generate normals
		for (int y=0; y<GEOPLATE_WALL_LEN-1; ++y) {
			for (int x=0; x<GEOPLATE_WALL_LEN-1; ++x) {
				// normal
				vector3d xy = vertices[x + y*GEOPLATE_WALL_LEN];
				vector3d x1 = vertices[x+1 + y*GEOPLATE_WALL_LEN];
				vector3d y1 = vertices[x + (y+1)*GEOPLATE_WALL_LEN];

				vector3d n = (x1-xy).Cross(y1-xy);
				normals[x + y*GEOPLATE_WALL_LEN] = (n.Normalized());
			}
		}

		// Generate bottom row normals
		for (int y=0; y<GEOPLATE_WALL_LEN-1; ++y) {
			const int x=GEOPLATE_WALL_LEN-1;
			// normal
			vector3d xy = vertices[x + y*GEOPLATE_WALL_LEN];
			vector3d x1 = vertices[x-1 + y*GEOPLATE_WALL_LEN];
			vector3d y1 = vertices[x + (y+1)*GEOPLATE_WALL_LEN];

			vector3d n = (xy-x1).Cross(y1-xy);
			normals[x + y*GEOPLATE_WALL_LEN] = (n.Normalized());
		}

		// Generate last col normals
		{
			const int y=GEOPLATE_WALL_LEN-1;
			for (int x=1; x<GEOPLATE_WALL_LEN; ++x) {
				// normal
				vector3d xy = vertices[x + y*GEOPLATE_WALL_LEN];
				vector3d x1 = vertices[x-1 + y*GEOPLATE_WALL_LEN];
				vector3d y1 = vertices[x + (y-1)*GEOPLATE_WALL_LEN];

				vector3d n = (xy-x1).Cross(xy-y1);
				normals[x + y*GEOPLATE_WALL_LEN] = (n.Normalized());
			}
		}

		// corner normals
		{
			const int y=GEOPLATE_WALL_LEN-1;
			const int x=0;
			{
				// normal
				vector3d xy = vertices[x + y*GEOPLATE_WALL_LEN];
				vector3d x1 = vertices[x+1 + y*GEOPLATE_WALL_LEN];
				vector3d y1 = vertices[x + (y-1)*GEOPLATE_WALL_LEN];

				vector3d n = (x1-xy).Cross(xy-y1);
				normals[x + y*GEOPLATE_WALL_LEN] = (n.Normalized());
			}
		}
		{
			const int y=0;
			const int x=GEOPLATE_WALL_LEN-1;
			{
				// normal
				vector3d xy = vertices[x + y*GEOPLATE_WALL_LEN];
				vector3d x1 = vertices[x-1 + y*GEOPLATE_WALL_LEN];
				vector3d y1 = vertices[x + (y+1)*GEOPLATE_WALL_LEN];

				vector3d n = (xy-x1).Cross(y1-xy);
				normals[x + y*GEOPLATE_WALL_LEN] = (n.Normalized());
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

		Pi::statSceneTris += (GEOPLATE_WALL_LEN-1)*(GEOPLATE_WALL_LEN-1)*2;
		++Pi::statNumPatches;

		// per-patch detail texture scaling value
		geoRing->GetMaterialParameters().patchDepth = 0;

		renderer->DrawBufferIndexed(m_vertexBuffer.get(), indices.Get(), rs, mat.Get());

		renderer->GetStats().AddToStatCount(Graphics::Stats::STAT_PLATES, 1);
	}
};
//static 
RefCountedPtr<Graphics::IndexBuffer> GeoPlateWall::indices;

class GeoPlate : public BaseGeo {
public:
	vector3d vbe[NUM_VBE];
	GeoPlate *kids[4];
	int m_depth;
	bool m_needUpdateVBOs;

	static RefCountedPtr<Graphics::IndexBuffer> indices;
	static int prevEdgeLen;
	static Uint32 indices_vbo;
	static VBOVertex *vbotemp;

	// params
	// v0, v1 - define points on the line describing the loop of the ring/orbital.
	// depth - 0 is the topmost plate with each depth+1 describing it's depth within the tree.
	GeoPlate(const double halfLength, const vector3d &startVBE, const vector3d &endVBE, 
		const double startAng, const double endAng, 
		const double yoffset, const int depth) 
		: BaseGeo(halfLength, startAng, endAng, yoffset, depth), m_needUpdateVBOs(false)
	{
		for( int i=0; i<4; ++i )
		{
			kids[i] = nullptr;
		}

		vbe[0] = startVBE;
		vbe[1] = endVBE;
		m_depth = depth;
		m_needUpdateVBOs = false;
		heights = new double[GEOPLATE_NUMVERTICES];
		normals = new vector3d[GEOPLATE_NUMVERTICES];
		vertices = new vector3d[GEOPLATE_NUMVERTICES];
		colors = new Color3ub[GEOPLATE_NUMVERTICES];
	}

	~GeoPlate() {
		for (int i=0; i<4; i++) if (kids[i]) delete kids[i];
		SAFE_DELETE_ARRAY(heights);
		SAFE_DELETE_ARRAY(vertices);
		SAFE_DELETE_ARRAY(normals);
		SAFE_DELETE_ARRAY(colors);
	}

	static void Init() 
	{
		GEOPLATE_FRAC = 1.0 / double(GEOPLATE_EDGELEN-1);

		//
		Uint32 *idx;
		std::vector<Uint32> pl_short;

		int tri_count = 0;
		{
			// calculate how many tri's there are
			tri_count = (VBO_COUNT_MID_IDX / 3);
			for (int i = 0; i<4; ++i) {
				tri_count += (VBO_COUNT_HI_EDGE / 3);
			}

			// pre-allocate enough space
			pl_short.reserve(tri_count);

			// add all of the middle indices
			for (int i = 0; i<VBO_COUNT_MID_IDX; ++i) {
				pl_short.push_back(0);
			}
			// add the HI detail indices
			for (int i = 0; i<4; i++) {
				for (int j = 0; j<VBO_COUNT_HI_EDGE; ++j) {
					pl_short.push_back(0);
				}
			}
		}
		// want vtx indices for tris
		idx = &pl_short[0];
		for (int x = 0; x<GEOPLATE_EDGELEN - 1; x++) {
			for (int y = 0; y<GEOPLATE_EDGELEN - 1; y++) {
				// 1st tri
				idx[0] = x + GEOPLATE_EDGELEN*y;
				idx[1] = x + 1 + GEOPLATE_EDGELEN*y;
				idx[2] = x + GEOPLATE_EDGELEN*(y + 1);
				idx += 3;

				// 2nd tri
				idx[0] = x + 1 + GEOPLATE_EDGELEN*y;
				idx[1] = x + 1 + GEOPLATE_EDGELEN*(y + 1);
				idx[2] = x + GEOPLATE_EDGELEN*(y + 1);
				idx += 3;
			}
		}

		// populate the N indices lists from the arrays built during InitTerrainIndices()
		// iterate over each index list and optimize it
		{
			VertexCacheOptimizerUInt vco;
			VertexCacheOptimizerUInt::Result res = vco.Optimize(&pl_short[0], tri_count);
			assert(0 == res);
			//create buffer & copy
			indices.Reset(Pi::renderer->CreateIndexBuffer(pl_short.size(), Graphics::BUFFER_USAGE_STATIC));
			Uint32* idxPtr = indices->Map(Graphics::BUFFER_MAP_WRITE);
			for (Uint32 j = 0; j < pl_short.size(); j++) {
				idxPtr[j] = pl_short[j];
			}
			indices->Unmap();
		}

		prevEdgeLen = GEOPLATE_EDGELEN;
	}

	void UpdateVBOs() 
	{
		m_needUpdateVBOs = true;
	}

	void _UpdateVBOs() 
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
			vbd.numVertices = GEOPLATE_NUMVERTICES;
			vbd.usage = Graphics::BUFFER_USAGE_STATIC;
			m_vertexBuffer.reset(Pi::renderer->CreateVertexBuffer(vbd));

			VBOVertex* VBOVtxPtr = m_vertexBuffer->Map<VBOVertex>(Graphics::BUFFER_MAP_WRITE);
			assert(m_vertexBuffer->GetDesc().stride == sizeof(VBOVertex));

			const vector3d *vtx = vertices;
			const vector3d *pNorm = normals;
			const Color3ub *pColr = colors;

			for (Sint32 y = 0; y<GEOPLATE_EDGELEN; y++) {
				for (Sint32 x = 0; x<GEOPLATE_EDGELEN; x++) {
					const double xFrac = double(x) * GEOPLATEHULL_FRAC;
					const double yFrac = double(y) * GEOPLATEHULL_FRAC;

					VBOVertex* vtxPtr = &VBOVtxPtr[x + (y*GEOPLATE_EDGELEN)];
					vector3d p(*vtx - clipCentroid);
					vtxPtr->pos = vector3f(p);
					clipRadius = std::max(clipRadius, p.Length());
					++vtx;	// next height

					vtxPtr->norm = vector3f(pNorm->Normalized());
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

			// ----------------------------------------------------
			// end of mapping
			m_vertexBuffer->Unmap();
		}
	}
	inline void SetColour(Color3ub *r, const vector3d &v) { 
#if 1
		r->r=static_cast<unsigned char>(Clamp(v.x*255.0, 0.0, 255.0)); 
		r->g=static_cast<unsigned char>(Clamp(v.y*255.0, 0.0, 255.0)); 
		r->b=static_cast<unsigned char>(Clamp(v.z*255.0, 0.0, 255.0));
#else
		r->r = 255;
		r->g = 255;
		r->b = 255;
#endif
	}
	// Generates full-detail vertices, and also non-edge normals and colors
	void GenerateMesh() 
	{
		double *hgt = heights;
		vector3d *vts = vertices;
		vector3d *nrm = normals;
		Color3ub *col = colors;
		double xfrac;
		double zfrac = 0;
		vector3d axisPt(0.0, 0.0, 0.0);

		const vector3d topEnd(0.0, m_yoffset + m_halfLen, 0.0);		// vertices at top edge of circle
		const vector3d btmEnd(0.0, m_yoffset - m_halfLen, 0.0);		// vertices at bottom edge of circle
				
		for (int y=0; y<GEOPLATE_EDGELEN; ++y) 
		{
			xfrac = 0;
			// point along z-axis by zfrac amount
			const vector3d axisPt = lerp( zfrac, btmEnd, topEnd );
			for (int x=0; x<GEOPLATE_EDGELEN; ++x) 
			{
				// find point on _surface_ of the cylinder
				const vector3d p = GetSurfacePointCyl(xfrac, zfrac, m_halfLen);
				// vector from axis to point-on-surface
				const vector3d cDir = (p - axisPt).Normalized();
				// height
				double height = -geoRing->GetHeight(p);
				// vertex is moved in direction of point-in-axis FROM point-on-surface by height.
				*(vts++) = p + (cDir * height);

				// remember this -- we will need it later
				*(hgt++) = -height;

				xfrac += GEOPLATE_FRAC;

				*(nrm++) = (cDir);

				*(col++) = Color3ub::BLACK;
			}
			zfrac += GEOPLATE_FRAC;
		}
		assert(vts == &vertices[GEOPLATE_NUMVERTICES]);
		// Generate normals & colors for non-edge vertices since they never change
		for (int y=1; y<GEOPLATE_EDGELEN-1; y++) 
		{
			for (int x=1; x<GEOPLATE_EDGELEN-1; x++) 
			{
				// normal
				vector3d x1 = vertices[x-1 + y*GEOPLATE_EDGELEN];
				vector3d x2 = vertices[x+1 + y*GEOPLATE_EDGELEN];
				vector3d y1 = vertices[x + (y-1)*GEOPLATE_EDGELEN];
				vector3d y2 = vertices[x + (y+1)*GEOPLATE_EDGELEN];

				vector3d n = (x2-x1).Cross(y2-y1);
				const vector3d &norm = normals[x + y*GEOPLATE_EDGELEN] = (n.Normalized());
				// color
				vector3d p = GetSurfacePointCyl(x*GEOPLATE_FRAC, y*GEOPLATE_FRAC, m_halfLen);
				const double height = heights[x + y*GEOPLATE_EDGELEN];
				SetColour(&colors[x + y*GEOPLATE_EDGELEN], geoRing->GetColor(p, height, norm));
			}
		}
	}
	
	void Render(Graphics::Renderer *renderer, const vector3d &campos, const matrix4x4d &modelView, const Graphics::Frustum &frustum) 
	{
		// Do frustum culling
		if (!frustum.TestPoint(clipCentroid, clipRadius))
			return; // nothing below this patch is visible

		if (kids[0]) {
			for (int i=0; i<4; i++) {
				kids[i]->Render(renderer, campos, modelView, frustum);
			}
		} else {
			_UpdateVBOs();

			RefCountedPtr<Graphics::Material> mat = geoRing->GetSurfaceMaterial();
			Graphics::RenderState *rs = geoRing->GetSurfRenderState();

			const vector3d relpos = clipCentroid - campos;
			renderer->SetTransform(modelView * matrix4x4d::Translation(relpos));

			Pi::statSceneTris += (VBO_COUNT_ALL_TRIS);
			++Pi::statNumPatches;

			// per-patch detail texture scaling value
			geoRing->GetMaterialParameters().patchDepth = m_depth;

			renderer->DrawBufferIndexed(m_vertexBuffer.get(), indices.Get(), rs, mat.Get());

			renderer->GetStats().AddToStatCount(Graphics::Stats::STAT_PLATES, 1);
		}
	}

	void LODUpdate(vector3d &campos) 
	{
		vector3d centroid = GetSurfacePointCyl(0.5, 0.5, m_halfLen);
		centroid = (1.0 + (-geoRing->GetHeight(centroid))) * centroid;

		bool canSplit = true;
		if (!(canSplit && (m_depth < GEOPLATE_MAX_DEPTH) &&
		    ((campos - centroid).Length() < m_roughLength))) {
			canSplit = false;
		}

		bool canMerge = true;

		if (canSplit) {
			if (!kids[0]) {
				const vector3d halfVBE = vbe[0] + 0.5 * (vbe[1] - vbe[0]);

				const double halfAng = ang[0] + 0.5 * (ang[1] - ang[0]);
				const double newLength = (m_halfLen*0.5);
				const double yoffset0 = m_yoffset + (newLength*0.5);
				const double yoffset1 = m_yoffset - (newLength*0.5);
				GeoPlate *_kids[4];
				_kids[0] = new GeoPlate(newLength, halfVBE, vbe[1], halfAng, ang[1], yoffset1, m_depth+1);
				_kids[1] = new GeoPlate(newLength, vbe[0], halfVBE, ang[0], halfAng, yoffset1, m_depth+1);
				_kids[2] = new GeoPlate(newLength, vbe[0], halfVBE, ang[0], halfAng, yoffset0, m_depth+1);
				_kids[3] = new GeoPlate(newLength, halfVBE, vbe[1], halfAng, ang[1], yoffset0, m_depth+1);

				_kids[0]->geoRing = _kids[1]->geoRing = _kids[2]->geoRing = _kids[3]->geoRing = geoRing;

				for (int i=0; i<4; ++i) _kids[i]->GenerateMesh();
				for (int i=0; i<4; ++i) kids[i] = _kids[i];
				for (int i=0; i<4; ++i) {
					kids[i]->UpdateVBOs();
				}
			}
			for (int i=0; i<4; ++i) kids[i]->LODUpdate(campos);
		} else {
			if (canMerge && kids[0]) {
				for (int i=0; i<4; ++i) { 
					delete kids[i]; kids[i] = 0; 
				}
			}
		}
	}
};
//static 
RefCountedPtr<Graphics::IndexBuffer> GeoPlate::indices;
int GeoPlate::prevEdgeLen = 0;
Uint32 GeoPlate::indices_vbo = 0;
VBOVertex *GeoPlate::vbotemp;

#define PLANET_AMBIENT	0.0
static std::vector<GeoRing*> s_allGeoRings;

void GeoRing::_UpdateLODs()
{
	for (size_t i=0; i<m_plates.size(); i++) {
		m_plates[i]->LODUpdate(m_tempCampos);
	}
}

void GeoRing::Init()
{
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

	switch (Pi::detail.planets) {
		case 0: GEOPLATE_EDGELEN = 7;  break;
		case 1: GEOPLATE_EDGELEN = 15; break;
		case 2: GEOPLATE_EDGELEN = 25; break;
		case 3: GEOPLATE_EDGELEN = 35; break;
		default:
		case 4: GEOPLATE_EDGELEN = GEOPLATE_MAX_EDGELEN; break;
	}
	assert(GEOPLATE_EDGELEN <= GEOPLATE_MAX_EDGELEN);
	GeoPlate::Init();
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
	for( int i=0 ; i<points.size()-1 ; ++i ) {
		m_plates.push_back( new GeoPlate(mRingWidth, points[i], points[i+1], angles[i], angles[i+1], 0.0, 0) );
	}

	// set edge friends
	for (size_t i=0; i<m_plates.size(); i++) {
		m_plates[i]->SetGeoRingPtr( this );
	}

	// create mesh(es)
	for (size_t i=0; i<m_plates.size(); i++) m_plates[i]->GenerateMesh();
	for (size_t i=0; i<m_plates.size(); i++) m_plates[i]->UpdateVBOs();

	// create the outer hull
	m_hull.clear();
	for( int i=0 ; i<points.size()-1 ; ++i ) {
		double len = (points[i+1] - points[i]).Length();
		m_hull.push_back( new GeoPlateHull(len, angles[i+1], angles[i], 0.0, 0) );
	}
	for (size_t i=0; i<m_hull.size(); i++) m_hull[i]->SetGeoRingPtr( this );
	for (size_t i=0; i<m_hull.size(); i++) m_hull[i]->GenerateMesh();
	for (size_t i=0; i<m_hull.size(); i++) m_hull[i]->UpdateVBOs();

	// create the inner wall
	m_wallInner.clear();
	for( int i=0 ; i<points.size()-1 ; ++i ) {
		double len = (points[i+1] - points[i]).Length();
		m_wallInner.push_back( new GeoPlateWall(true, len, angles[i+1], angles[i], 0.0, 0) );
	}
	for (size_t i=0; i<m_wallInner.size(); i++) m_wallInner[i]->SetGeoRingPtr( this );
	for (size_t i=0; i<m_wallInner.size(); i++) m_wallInner[i]->GenerateMesh();
	for (size_t i=0; i<m_wallInner.size(); i++) m_wallInner[i]->UpdateVBOs();

	// create the outer wall
	m_wallOuter.clear();
	for( int i=0 ; i<points.size()-1 ; ++i ) {
		double len = (points[i+1] - points[i]).Length();
		m_wallOuter.push_back( new GeoPlateWall(false, len, angles[i+1], angles[i], 0.0, 0) );
	}
	for (size_t i=0; i<m_wallOuter.size(); i++) m_wallOuter[i]->SetGeoRingPtr( this );
	for (size_t i=0; i<m_wallOuter.size(); i++) m_wallOuter[i]->GenerateMesh();
	for (size_t i=0; i<m_wallOuter.size(); i++) m_wallOuter[i]->UpdateVBOs();
}

static const double gs_targetPatchTriLength(100.0);
void GeoRing::CalculateMaxPatchDepth()
{
	m_maxDepth = 0;
	const double circumference = 2.0 * M_PI * m_sbody->GetRadius();
	// calculate length of each edge segment (quad) times 4 due to that being the number around the sphere (1 per side, 4 sides for Root).
	double edgeMetres = circumference / double(GEOPLATE_EDGELEN_DEFAULT * 8);
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
	/*switch(m_initStage)
	{
	case eBuildFirstPatches:
		BuildFirstPatches();
		break;
	case eRequestedFirstPatches:
		{
			ProcessSplitResults();
			uint8_t numValidPatches = 0;
			for (int i=0; i<NUM_PATCHES; i++) {
				if(m_patches[i]->HasHeightData()) {
					++numValidPatches;
				}
			}
			m_initStage = (NUM_PATCHES==numValidPatches) ? eReceivedFirstPatches : eRequestedFirstPatches;
		} break;
	case eReceivedFirstPatches:
		{
			for (int i=0; i<NUM_PATCHES; i++) {
				m_patches[i]->NeedToUpdateVBOs();
			}
			m_initStage = eDefaultUpdateState;
		} break;
	case eDefaultUpdateState:
		if(m_hasTempCampos) {
			ProcessSplitResults();
			for (int i=0; i<NUM_PATCHES; i++) {
				m_patches[i]->LODUpdate(m_tempCampos, m_tempFrustum);
			}
			ProcessQuadSplitRequests();
		}
		break;
	}*/
}

static const float g_ambient[4] = { 0, 0, 0, 1.0 };

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

