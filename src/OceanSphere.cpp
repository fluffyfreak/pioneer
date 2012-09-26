#include "libs.h"
#include "OceanSphere.h"
#include "perlin.h"
#include "Pi.h"
#include "RefCounted.h"
#include "graphics/Material.h"
#include "graphics/Renderer.h"
#include "graphics/Frustum.h"
#include "graphics/Graphics.h"
#include "graphics/VertexArray.h"
#include "graphics/gl2/GeoSphereMaterial.h"
#include "vcacheopt/vcacheopt.h"
#include <deque>
#include <algorithm>

// tri edge lengths
#define OCEANPATCH_SUBDIVIDE_AT_CAMDIST	5.0
#define OCEANPATCH_MAX_DEPTH  15 + (2*Pi::detail.fracmult) //15

static const int OCEANPATCH_MAX_EDGELEN = 55;
int OceanSphere::s_vtxGenCount = 0;
RefCountedPtr<OceanPatchContext> OceanSphere::s_patchContext;

// must be odd numbers
static const int detail_edgeLen[5] = {
	7, 15, 25, 35, 55
};

#define PRINT_VECTOR(_v) printf("%f,%f,%f\n", (_v).x, (_v).y, (_v).z);

#pragma pack(4)
struct VBOVertex
{
	float x,y,z;
	float nx,ny,nz;
	float padding;
	float padding2;
};
#pragma pack()

// hold the 16 possible terrain edge connections
const int NUM_INDEX_LISTS = 16;

class OceanPatchContext : public RefCounted {
public:
	const int edgeLen;
	const int halfEdgeLen;

	inline int VBO_COUNT_LO_EDGE() const { return 3*(halfEdgeLen); }
	inline int VBO_COUNT_HI_EDGE() const { return 3*(edgeLen-1); }
	inline int VBO_COUNT_MID_IDX() const { return (4*3*(edgeLen-3))    + 2*(edgeLen-3)*(edgeLen-3)*3; }
	//                                            ^^ serrated teeth bit  ^^^ square inner bit

	inline int IDX_VBO_LO_OFFSET(int i) const { return i*sizeof(unsigned short)*3*(halfEdgeLen); }
	inline int IDX_VBO_HI_OFFSET(int i) const { return (i*sizeof(unsigned short)*VBO_COUNT_HI_EDGE())+IDX_VBO_LO_OFFSET(4); }
	inline int IDX_VBO_MAIN_OFFSET()    const { return IDX_VBO_HI_OFFSET(4); }
	inline int IDX_VBO_COUNT_ALL_IDX()	const { return ((edgeLen-1)*(edgeLen-1))*2*3; }

	inline int NUMVERTICES() const { return edgeLen*edgeLen; }

	double frac;

	ScopedArray<unsigned short> midIndices;
	ScopedArray<unsigned short> loEdgeIndices[4];
	ScopedArray<unsigned short> hiEdgeIndices[4];
	GLuint indices_vbo;
	GLuint indices_list[NUM_INDEX_LISTS];
	GLuint indices_tri_count;
	GLuint indices_tri_counts[NUM_INDEX_LISTS];
	VBOVertex *vbotemp;

	OceanPatchContext(const int _edgeLen) : edgeLen(_edgeLen), halfEdgeLen(_edgeLen>>1) {
		Init();
	}

	~OceanPatchContext() {
		Cleanup();
	}

	void Refresh() {
		Cleanup();
		Init();
	}

	void Cleanup() {
		midIndices.Reset();
		for (int i=0; i<4; i++) {
			loEdgeIndices[i].Reset();
			hiEdgeIndices[i].Reset();
		}
		if (indices_vbo) {
			indices_vbo = 0;
		}
		for (int i=0; i<NUM_INDEX_LISTS; i++) {
			if (indices_list[i]) {
				glDeleteBuffersARB(1, &indices_list[i]);
			}
		}
		delete [] vbotemp;
	}

	void updateIndexBufferId(const GLuint edge_hi_flags) {
		assert(edge_hi_flags < GLuint(NUM_INDEX_LISTS));
		indices_vbo = indices_list[edge_hi_flags];
		indices_tri_count = indices_tri_counts[edge_hi_flags];
	}

	int getIndices(std::vector<unsigned short> &pl, const unsigned int edge_hi_flags)
	{
		// calculate how many tri's there are
		int tri_count = (VBO_COUNT_MID_IDX() / 3); 
		for( int i=0; i<4; ++i ) {
			if( edge_hi_flags & (1 << i) ) {
				tri_count += (VBO_COUNT_HI_EDGE() / 3);
			} else {
				tri_count += (VBO_COUNT_LO_EDGE() / 3);
			}
		}

		// pre-allocate enough space
		pl.reserve(tri_count);

		// add all of the middle indices
		for(int i=0; i<VBO_COUNT_MID_IDX(); ++i) {
			pl.push_back(midIndices[i]);
		}
		// selectively add the HI or LO detail indices
		for (int i=0; i<4; i++) {
			if( edge_hi_flags & (1 << i) ) {
				for(int j=0; j<VBO_COUNT_HI_EDGE(); ++j) {
					pl.push_back(hiEdgeIndices[i][j]);
				}
			} else {
				for(int j=0; j<VBO_COUNT_LO_EDGE(); ++j) {
					pl.push_back(loEdgeIndices[i][j]);
				}
			}
		}

		return tri_count;
	}

	void Init() {
		frac = 1.0 / double(edgeLen-1);

		vbotemp = new VBOVertex[NUMVERTICES()];

		unsigned short *idx;
		midIndices.Reset(new unsigned short[VBO_COUNT_MID_IDX()]);
		for (int i=0; i<4; i++) {
			loEdgeIndices[i].Reset(new unsigned short[VBO_COUNT_LO_EDGE()]);
			hiEdgeIndices[i].Reset(new unsigned short[VBO_COUNT_HI_EDGE()]);
		}
		/* also want vtx indices for tris not touching edge of patch */
		idx = midIndices.Get();
		for (int x=1; x<edgeLen-2; x++) {
			for (int y=1; y<edgeLen-2; y++) {
				idx[0] = x + edgeLen*y;
				idx[1] = x+1 + edgeLen*y;
				idx[2] = x + edgeLen*(y+1);
				idx+=3;

				idx[0] = x+1 + edgeLen*y;
				idx[1] = x+1 + edgeLen*(y+1);
				idx[2] = x + edgeLen*(y+1);
				idx+=3;
			}
		}
		{
			for (int x=1; x<edgeLen-3; x+=2) {
				// razor teeth near edge 0
				idx[0] = x + edgeLen;
				idx[1] = x+1;
				idx[2] = x+1 + edgeLen;
				idx+=3;
				idx[0] = x+1;
				idx[1] = x+2 + edgeLen;
				idx[2] = x+1 + edgeLen;
				idx+=3;
			}
			for (int x=1; x<edgeLen-3; x+=2) {
				// near edge 2
				idx[0] = x + edgeLen*(edgeLen-2);
				idx[1] = x+1 + edgeLen*(edgeLen-2);
				idx[2] = x+1 + edgeLen*(edgeLen-1);
				idx+=3;
				idx[0] = x+1 + edgeLen*(edgeLen-2);
				idx[1] = x+2 + edgeLen*(edgeLen-2);
				idx[2] = x+1 + edgeLen*(edgeLen-1);
				idx+=3;
			}
			for (int y=1; y<edgeLen-3; y+=2) {
				// near edge 1
				idx[0] = edgeLen-2 + y*edgeLen;
				idx[1] = edgeLen-1 + (y+1)*edgeLen;
				idx[2] = edgeLen-2 + (y+1)*edgeLen;
				idx+=3;
				idx[0] = edgeLen-2 + (y+1)*edgeLen;
				idx[1] = edgeLen-1 + (y+1)*edgeLen;
				idx[2] = edgeLen-2 + (y+2)*edgeLen;
				idx+=3;
			}
			for (int y=1; y<edgeLen-3; y+=2) {
				// near edge 3
				idx[0] = 1 + y*edgeLen;
				idx[1] = 1 + (y+1)*edgeLen;
				idx[2] = (y+1)*edgeLen;
				idx+=3;
				idx[0] = 1 + (y+1)*edgeLen;
				idx[1] = 1 + (y+2)*edgeLen;
				idx[2] = (y+1)*edgeLen;
				idx+=3;
			}
		}
		// full detail edge triangles
		{
			idx = hiEdgeIndices[0].Get();
			for (int x=0; x<edgeLen-1; x+=2) {
				idx[0] = x; idx[1] = x+1; idx[2] = x+1 + edgeLen;
				idx+=3;
				idx[0] = x+1; idx[1] = x+2; idx[2] = x+1 + edgeLen;
				idx+=3;
			}
			idx = hiEdgeIndices[1].Get();
			for (int y=0; y<edgeLen-1; y+=2) {
				idx[0] = edgeLen-1 + y*edgeLen;
				idx[1] = edgeLen-1 + (y+1)*edgeLen;
				idx[2] = edgeLen-2 + (y+1)*edgeLen;
				idx+=3;
				idx[0] = edgeLen-1 + (y+1)*edgeLen;
				idx[1] = edgeLen-1 + (y+2)*edgeLen;
				idx[2] = edgeLen-2 + (y+1)*edgeLen;
				idx+=3;
			}
			idx = hiEdgeIndices[2].Get();
			for (int x=0; x<edgeLen-1; x+=2) {
				idx[0] = x + (edgeLen-1)*edgeLen;
				idx[1] = x+1 + (edgeLen-2)*edgeLen;
				idx[2] = x+1 + (edgeLen-1)*edgeLen;
				idx+=3;
				idx[0] = x+1 + (edgeLen-2)*edgeLen;
				idx[1] = x+2 + (edgeLen-1)*edgeLen;
				idx[2] = x+1 + (edgeLen-1)*edgeLen;
				idx+=3;
			}
			idx = hiEdgeIndices[3].Get();
			for (int y=0; y<edgeLen-1; y+=2) {
				idx[0] = y*edgeLen;
				idx[1] = 1 + (y+1)*edgeLen;
				idx[2] = (y+1)*edgeLen;
				idx+=3;
				idx[0] = (y+1)*edgeLen;
				idx[1] = 1 + (y+1)*edgeLen;
				idx[2] = (y+2)*edgeLen;
				idx+=3;
			}
		}
		// these edge indices are for patches with no
		// neighbour of equal or greater detail -- they reduce
		// their edge complexity by 1 division
		{
			idx = loEdgeIndices[0].Get();
			for (int x=0; x<edgeLen-2; x+=2) {
				idx[0] = x;
				idx[1] = x+2;
				idx[2] = x+1+edgeLen;
				idx += 3;
			}
			idx = loEdgeIndices[1].Get();
			for (int y=0; y<edgeLen-2; y+=2) {
				idx[0] = (edgeLen-1) + y*edgeLen;
				idx[1] = (edgeLen-1) + (y+2)*edgeLen;
				idx[2] = (edgeLen-2) + (y+1)*edgeLen;
				idx += 3;
			}
			idx = loEdgeIndices[2].Get();
			for (int x=0; x<edgeLen-2; x+=2) {
				idx[0] = x+edgeLen*(edgeLen-1);
				idx[2] = x+2+edgeLen*(edgeLen-1);
				idx[1] = x+1+edgeLen*(edgeLen-2);
				idx += 3;
			}
			idx = loEdgeIndices[3].Get();
			for (int y=0; y<edgeLen-2; y+=2) {
				idx[0] = y*edgeLen;
				idx[2] = (y+2)*edgeLen;
				idx[1] = 1 + (y+1)*edgeLen;
				idx += 3;
			}
		}

		// these will hold the optimised indices
		std::vector<unsigned short> pl_short[NUM_INDEX_LISTS];
		// populate the N indices lists from the arrays built during InitTerrainIndices()
		for( int i=0; i<NUM_INDEX_LISTS; ++i ) {
			const unsigned int edge_hi_flags = i;
			indices_tri_counts[i] = getIndices(pl_short[i], edge_hi_flags);
		}

		// iterate over each index list and optimize it
		for( int i=0; i<NUM_INDEX_LISTS; ++i ) {
			int tri_count = indices_tri_counts[i];
			VertexCacheOptimizerUShort vco;
			VertexCacheOptimizerUShort::Result res = vco.Optimize(&pl_short[i][0], tri_count);
			assert(0 == res);
		}

		// everything should be hunky-dory for setting up as OpenGL index buffers now.
		for( int i=0; i<NUM_INDEX_LISTS; ++i ) {
			glGenBuffersARB(1, &indices_list[i]);
			glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, indices_list[i]);
			glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned short)*indices_tri_counts[i]*3, &(pl_short[i][0]), GL_STATIC_DRAW);
			glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);
		}

		// default it to the last entry which uses the hi-res borders
		indices_vbo			= indices_list[NUM_INDEX_LISTS-1];
		indices_tri_count	= indices_tri_counts[NUM_INDEX_LISTS-1];

		if (midIndices) {
			midIndices.Reset();
			for (int i=0; i<4; i++) {
				loEdgeIndices[i].Reset();
				hiEdgeIndices[i].Reset();
			}
		}
	}

	void GetEdge(vector3d *array, int edge, vector3d *ev) {
		if (edge == 0) {
			for (int x=0; x<edgeLen; x++) ev[x] = array[x];
		} else if (edge == 1) {
			const int x = edgeLen-1;
			for (int y=0; y<edgeLen; y++) ev[y] = array[x + y*edgeLen];
		} else if (edge == 2) {
			const int y = edgeLen-1;
			for (int x=0; x<edgeLen; x++) ev[x] = array[(edgeLen-1)-x + y*edgeLen];
		} else {
			for (int y=0; y<edgeLen; y++) ev[y] = array[0 + ((edgeLen-1)-y)*edgeLen];
		}
	}

	void SetEdge(vector3d *array, int edge, const vector3d *ev) {
		if (edge == 0) {
			for (int x=0; x<edgeLen; x++) array[x] = ev[x];
		} else if (edge == 1) {
			const int x = edgeLen-1;
			for (int y=0; y<edgeLen; y++) array[x + y*edgeLen] = ev[y];
		} else if (edge == 2) {
			const int y = edgeLen-1;
			for (int x=0; x<edgeLen; x++) array[(edgeLen-1)-x + y*edgeLen] = ev[x];
		} else {
			for (int y=0; y<edgeLen; y++) array[0 + ((edgeLen-1)-y)*edgeLen] = ev[y];
		}
	}
};


class OceanPatch {
public:
	RefCountedPtr<OceanPatchContext> ctx;
	vector3d v[4];
	vector3d *vertices;
	vector3d *normals;
	GLuint m_vbo;
	OceanPatch *kids[4];
	OceanPatch *parent;
	OceanPatch *edgeFriend[4]; // [0]=v01, [1]=v12, [2]=v20
	OceanSphere *oceansphere;
	double m_roughLength;
	vector3d clipCentroid, centroid;
	double clipRadius;
	int m_depth;
	SDL_mutex *m_kidsLock;
	bool m_needUpdateVBOs;
	double m_distMult;

	OceanPatch(const RefCountedPtr<OceanPatchContext> &_ctx, OceanSphere *gs, vector3d v0, vector3d v1, vector3d v2, vector3d v3, int depth) {
		memset(this, 0, sizeof(OceanPatch));

		ctx = _ctx;

		oceansphere = gs;

		m_kidsLock = SDL_CreateMutex();
		v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3;
		//depth -= Pi::detail.fracmult;
		m_depth = depth;
		clipCentroid = (v0+v1+v2+v3) * 0.25;
		clipRadius = 0;
		for (int i=0; i<4; i++) {
			clipRadius = std::max(clipRadius, (v[i]-clipCentroid).Length());
		}
		if (oceansphere->m_sbody->type < SystemBody::TYPE_PLANET_ASTEROID) {
 			m_distMult = 10 / Clamp(depth, 1, 10);
 		} else {
 			m_distMult = 5 / Clamp(depth, 1, 5);
 		}
		m_roughLength = OCEANPATCH_SUBDIVIDE_AT_CAMDIST / pow(2.0, depth) * m_distMult;
		m_needUpdateVBOs = false;
		normals = new vector3d[ctx->NUMVERTICES()];
		vertices = new vector3d[ctx->NUMVERTICES()];
	}

	~OceanPatch() {
		SDL_DestroyMutex(m_kidsLock);
		for (int i=0; i<4; i++) {
			if (edgeFriend[i]) edgeFriend[i]->NotifyEdgeFriendDeleted(this);
		}
		for (int i=0; i<4; i++) if (kids[i]) delete kids[i];
		delete[] vertices;
		delete[] normals;
		oceansphere->AddVBOToDestroy(m_vbo);
	}

	void UpdateVBOs() {
		m_needUpdateVBOs = true;
	}

	void _UpdateVBOs() {
		if (m_needUpdateVBOs) {
			if (!m_vbo) glGenBuffersARB(1, &m_vbo);
			m_needUpdateVBOs = false;
			glBindBufferARB(GL_ARRAY_BUFFER, m_vbo);
			glBufferDataARB(GL_ARRAY_BUFFER, sizeof(VBOVertex)*ctx->NUMVERTICES(), 0, GL_DYNAMIC_DRAW);
			for (int i=0; i<ctx->NUMVERTICES(); i++)
			{
				clipRadius = std::max(clipRadius, (vertices[i]-clipCentroid).Length());
				VBOVertex *pData = ctx->vbotemp + i;
				pData->x = float(vertices[i].x - clipCentroid.x);
				pData->y = float(vertices[i].y - clipCentroid.y);
				pData->z = float(vertices[i].z - clipCentroid.z);
				pData->nx = float(normals[i].x);
				pData->ny = float(normals[i].y);
				pData->nz = float(normals[i].z);
			}
			glBufferDataARB(GL_ARRAY_BUFFER, sizeof(VBOVertex)*ctx->NUMVERTICES(), ctx->vbotemp, GL_DYNAMIC_DRAW);
			glBindBufferARB(GL_ARRAY_BUFFER, 0);
		}
	}
	/* not quite edge, since we share edge vertices so that would be
	 * fucking pointless. one position inwards. used to make edge normals
	 * for adjacent tiles */
	void GetEdgeMinusOneVerticesFlipped(int edge, vector3d *ev) const {
		if (edge == 0) {
			for (int x=0; x<ctx->edgeLen; x++) ev[ctx->edgeLen-1-x] = vertices[x + ctx->edgeLen];
		} else if (edge == 1) {
			const int x = ctx->edgeLen-2;
			for (int y=0; y<ctx->edgeLen; y++) ev[ctx->edgeLen-1-y] = vertices[x + y*ctx->edgeLen];
		} else if (edge == 2) {
			const int y = ctx->edgeLen-2;
			for (int x=0; x<ctx->edgeLen; x++) ev[ctx->edgeLen-1-x] = vertices[(ctx->edgeLen-1)-x + y*ctx->edgeLen];
		} else {
			for (int y=0; y<ctx->edgeLen; y++) ev[ctx->edgeLen-1-y] = vertices[1 + ((ctx->edgeLen-1)-y)*ctx->edgeLen];
		}
	}
	int GetEdgeIdxOf(const OceanPatch *e) const {
		for (int i=0; i<4; i++) {
			if (edgeFriend[i] == e) return i;
		}
		abort();
		return -1;
	}


	int GetChildIdx(OceanPatch *child) const {
		for (int i=0; i<4; i++) {
			if (kids[i] == child) return i;
		}
		abort();
		return -1;
	}

	void FixEdgeFromParentInterpolated(int edge) {
		// noticeable artefacts from not doing so...
		vector3d ev[OCEANPATCH_MAX_EDGELEN];
		vector3d en[OCEANPATCH_MAX_EDGELEN];
		vector3d ev2[OCEANPATCH_MAX_EDGELEN];
		vector3d en2[OCEANPATCH_MAX_EDGELEN];
		ctx->GetEdge(parent->vertices, edge, ev);
		ctx->GetEdge(parent->normals, edge, en);

		int kid_idx = parent->GetChildIdx(this);
		if (edge == kid_idx) {
			// use first half of edge
			for (int i=0; i<=ctx->halfEdgeLen; i++) {
				ev2[i<<1] = ev[i];
				en2[i<<1] = en[i];
			}
		} else {
			// use 2nd half of edge
			for (int i=ctx->halfEdgeLen; i<ctx->edgeLen; i++) {
				ev2[(i-(ctx->halfEdgeLen))<<1] = ev[i];
				en2[(i-(ctx->halfEdgeLen))<<1] = en[i];
			}
		}
		// interpolate!!
		for (int i=1; i<ctx->edgeLen; i+=2) {
			ev2[i] = (ev2[i-1]+ev2[i+1]) * 0.5;
			en2[i] = (en2[i-1]+en2[i+1]).Normalized();
		}
		ctx->SetEdge(this->vertices, edge, ev2);
		ctx->SetEdge(this->normals, edge, en2);
	}

	/* in patch surface coords, [0,1] */
	vector3d GetSpherePoint(double x, double y) {
		return (v[0] + x*(1.0-y)*(v[1]-v[0]) +
			    x*y*(v[2]-v[0]) +
			    (1.0-x)*y*(v[3]-v[0])).Normalized();
	}

	/** Generates full-detail vertices, and also non-edge normals */
	void GenerateMesh() {
		centroid = clipCentroid.Normalized();
		centroid = (1.0 + oceansphere->GetHeight(centroid)) * centroid;
		vector3d *vts = vertices;
		vector3d *nor = normals;
		double xfrac;
		double yfrac = 0;
		for (int y=0; y<ctx->edgeLen; y++) {
			xfrac = 0;
			for (int x=0; x<ctx->edgeLen; x++) {
				const vector3d p = GetSpherePoint(xfrac, yfrac);
				const double height = oceansphere->GetHeight(p);
				/*const double waveTime = 1.0;
				const double waveWidth = 0.1;
				const double waveHeight = 3.0;
			//	const double r = pow(p.x*p.x + p.y*p.y + p.z*p.z, 0.5);
				const double theta = atan(p.z/pow(p.x*p.x+p.y*p.y, 0.5));//atan(p.z/sqrt(p.x*p.x+p.y*p.y));
				const double phi = atan(p.y/p.x);
				const double height = (sin(waveWidth * (theta*1000.0) + waveTime) * cos(waveWidth * (phi*1000.0) + waveTime) * waveHeight) * 0.01;*/
				*(vts++) = p * (height + 1.0);
				*(nor++) = p;
				xfrac += ctx->frac;
			}
			yfrac += ctx->frac;
		}
		assert(vts == &vertices[ctx->NUMVERTICES()]);
		assert(nor == &normals[ctx->NUMVERTICES()]);
	}
	void OnEdgeFriendChanged(const int edge, OceanPatch *e) {
		edgeFriend[edge] = e;
		vector3d ev[OCEANPATCH_MAX_EDGELEN];
		const int we_are = e->GetEdgeIdxOf(this);
		e->GetEdgeMinusOneVerticesFlipped(we_are, ev);
		/* now we have a valid edge, fix the edge vertices */
		if (edge == 0) {
			for (int x=0; x<ctx->edgeLen; x++) {
				vector3d p = GetSpherePoint(x * ctx->frac, 0);
				double height = oceansphere->GetHeight(p);
				vertices[x] = p * (height + 1.0);
			}
		} else if (edge == 1) {
			for (int y=0; y<ctx->edgeLen; y++) {
				vector3d p = GetSpherePoint(1.0, y * ctx->frac);
				double height = oceansphere->GetHeight(p);
				int pos = (ctx->edgeLen-1) + y*ctx->edgeLen;
				vertices[pos] = p * (height + 1.0);
			}
		} else if (edge == 2) {
			for (int x=0; x<ctx->edgeLen; x++) {
				vector3d p = GetSpherePoint(x * ctx->frac, 1.0);
				double height = oceansphere->GetHeight(p);
				int pos = x + (ctx->edgeLen-1)*ctx->edgeLen;
				vertices[pos] = p * (height + 1.0);
			}
		} else {
			for (int y=0; y<ctx->edgeLen; y++) {
				vector3d p = GetSpherePoint(0, y * ctx->frac);
				double height = oceansphere->GetHeight(p);
				int pos = y * ctx->edgeLen;
				vertices[pos] = p * (height + 1.0);
			}
		}

		UpdateVBOs();

		if (kids[0]) {
			if (edge == 0) {
				kids[0]->FixEdgeFromParentInterpolated(0);
				kids[0]->UpdateVBOs();
				kids[1]->FixEdgeFromParentInterpolated(0);
				kids[1]->UpdateVBOs();
			} else if (edge == 1) {
				kids[1]->FixEdgeFromParentInterpolated(1);
				kids[1]->UpdateVBOs();
				kids[2]->FixEdgeFromParentInterpolated(1);
				kids[2]->UpdateVBOs();
			} else if (edge == 2) {
				kids[2]->FixEdgeFromParentInterpolated(2);
				kids[2]->UpdateVBOs();
				kids[3]->FixEdgeFromParentInterpolated(2);
				kids[3]->UpdateVBOs();
			} else {
				kids[3]->FixEdgeFromParentInterpolated(3);
				kids[3]->UpdateVBOs();
				kids[0]->FixEdgeFromParentInterpolated(3);
				kids[0]->UpdateVBOs();
			}
		}
	}
	void NotifyEdgeFriendSplit(OceanPatch *e) {
		int idx = GetEdgeIdxOf(e);
		int we_are = e->GetEdgeIdxOf(this);
		if (!kids[0]) return;
		// match e's new kids to our own... :/
		kids[idx]->OnEdgeFriendChanged(idx, e->kids[(we_are+1)%4]);
		kids[(idx+1)%4]->OnEdgeFriendChanged(idx, e->kids[we_are]);
	}

	void NotifyEdgeFriendDeleted(OceanPatch *e) {
		int idx = GetEdgeIdxOf(e);
		edgeFriend[idx] = 0;
		if (!parent) return;
		if (parent->edgeFriend[idx]) {
			FixEdgeFromParentInterpolated(idx);
			UpdateVBOs();
		} else {
			// XXX TODO XXX
			// Bad. not fixing up edges in this case!!!
		}
	}

	OceanPatch *GetEdgeFriendForKid(const int kid, const int edge) const {
		OceanPatch *e = edgeFriend[edge];
		if (!e) return 0;
		//assert (e);// && (e->m_depth >= m_depth));
		const int we_are = e->GetEdgeIdxOf(this);
		// neighbour patch has not split yet (is at depth of this patch), so kids of this patch do
		// not have same detail level neighbours yet
		if (edge == kid) return e->kids[(we_are+1)%4];
		else return e->kids[we_are];
	}

	GLuint determineIndexbuffer() const {
		return // index buffers are ordered by edge resolution flags
			(edgeFriend[0] ? 1u : 0u) |
			(edgeFriend[1] ? 2u : 0u) |
			(edgeFriend[2] ? 4u : 0u) |
			(edgeFriend[3] ? 8u : 0u);
	}

	void Render(const vector3d &campos, const Graphics::Frustum &frustum) {
		PiVerify(SDL_mutexP(m_kidsLock)==0);
		if (kids[0]) {
			for (int i=0; i<4; i++) kids[i]->Render(campos, frustum);
			SDL_mutexV(m_kidsLock);
		} else {
			SDL_mutexV(m_kidsLock);
			_UpdateVBOs();

			if (!frustum.TestPoint(clipCentroid, clipRadius))
				return;

			vector3d relpos = clipCentroid - campos;
			glPushMatrix();
			glTranslated(relpos.x, relpos.y, relpos.z);

			Pi::statSceneTris += 2*(ctx->edgeLen-1)*(ctx->edgeLen-1);
			glEnableClientState(GL_VERTEX_ARRAY);
			glEnableClientState(GL_NORMAL_ARRAY);

			// update the indices used for rendering
			ctx->updateIndexBufferId(determineIndexbuffer());

			glBindBufferARB(GL_ARRAY_BUFFER, m_vbo);
			glVertexPointer(3, GL_FLOAT, sizeof(VBOVertex), 0);
			glNormalPointer(GL_FLOAT, sizeof(VBOVertex), reinterpret_cast<void *>(3*sizeof(float)));
			glColorPointer(4, GL_UNSIGNED_BYTE, sizeof(VBOVertex), reinterpret_cast<void *>(6*sizeof(float)));
			glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, ctx->indices_vbo);
			glDrawElements(GL_TRIANGLES, ctx->indices_tri_count*3, GL_UNSIGNED_SHORT, 0);
			glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
			glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);

			glDisableClientState(GL_VERTEX_ARRAY);
			glDisableClientState(GL_NORMAL_ARRAY);
			glPopMatrix();
		}
	}

	void LODUpdate(const vector3d &campos) {
		bool canSplit = true;
		for (int i=0; i<4; i++) {
			if (!edgeFriend[i]) { canSplit = false; break; }
			if (edgeFriend[i] && (edgeFriend[i]->m_depth < m_depth)) {
				canSplit = false;
				break;
			}
		}
		if (!(canSplit && (m_depth < OCEANPATCH_MAX_DEPTH) &&
		    ((campos - centroid).Length() < m_roughLength)))
			canSplit = false;
		// always split at first level
		if (!parent) canSplit = true;
		//printf(campos.Length());

		bool canMerge = true;

		if (canSplit) {
			if (!kids[0]) {
				vector3d v01, v12, v23, v30, cn;
				cn = centroid.Normalized();
				v01 = (v[0]+v[1]).Normalized();
				v12 = (v[1]+v[2]).Normalized();
				v23 = (v[2]+v[3]).Normalized();
				v30 = (v[3]+v[0]).Normalized();
				OceanPatch *_kids[4];
				_kids[0] = new OceanPatch(ctx, oceansphere, v[0], v01, cn, v30, m_depth+1);
				_kids[1] = new OceanPatch(ctx, oceansphere, v01, v[1], v12, cn, m_depth+1);
				_kids[2] = new OceanPatch(ctx, oceansphere, cn, v12, v[2], v23, m_depth+1);
				_kids[3] = new OceanPatch(ctx, oceansphere, v30, cn, v23, v[3], m_depth+1);
				// hm.. edges. Not right to pass this
				// edgeFriend...
				_kids[0]->edgeFriend[0] = GetEdgeFriendForKid(0, 0);
				_kids[0]->edgeFriend[1] = _kids[1];
				_kids[0]->edgeFriend[2] = _kids[3];
				_kids[0]->edgeFriend[3] = GetEdgeFriendForKid(0, 3);
				_kids[1]->edgeFriend[0] = GetEdgeFriendForKid(1, 0);
				_kids[1]->edgeFriend[1] = GetEdgeFriendForKid(1, 1);
				_kids[1]->edgeFriend[2] = _kids[2];
				_kids[1]->edgeFriend[3] = _kids[0];
				_kids[2]->edgeFriend[0] = _kids[1];
				_kids[2]->edgeFriend[1] = GetEdgeFriendForKid(2, 1);
				_kids[2]->edgeFriend[2] = GetEdgeFriendForKid(2, 2);
				_kids[2]->edgeFriend[3] = _kids[3];
				_kids[3]->edgeFriend[0] = _kids[0];
				_kids[3]->edgeFriend[1] = _kids[2];
				_kids[3]->edgeFriend[2] = GetEdgeFriendForKid(3, 2);
				_kids[3]->edgeFriend[3] = GetEdgeFriendForKid(3, 3);
				_kids[0]->parent = _kids[1]->parent = _kids[2]->parent = _kids[3]->parent = this;
				for (int i=0; i<4; i++) _kids[i]->GenerateMesh();
				PiVerify(SDL_mutexP(m_kidsLock)==0);
				for (int i=0; i<4; i++) kids[i] = _kids[i];
				for (int i=0; i<4; i++) edgeFriend[i]->NotifyEdgeFriendSplit(this);
				for (int i=0; i<4; i++) {
					kids[i]->UpdateVBOs();
				}
				PiVerify(SDL_mutexV(m_kidsLock)!=-1);
			}
			for (int i=0; i<4; i++) kids[i]->LODUpdate(campos);
		} else {
			if (canMerge && kids[0]) {
				PiVerify(SDL_mutexP(m_kidsLock)==0);
				for (int i=0; i<4; i++) { delete kids[i]; kids[i] = 0; }
				PiVerify(SDL_mutexV(m_kidsLock)!=-1);
			}
		}
	}
};

static const int ocean_sphere_edge_friends[6][4] = {
	{ 3, 4, 1, 2 },
	{ 0, 4, 5, 2 },
	{ 0, 1, 5, 3 },
	{ 0, 2, 5, 4 },
	{ 0, 3, 5, 1 },
	{ 1, 4, 3, 2 }
};

/* updates oceansphere level of detail thingies */
int OceanSphere::UpdateLODs(const vector3d campos)
{
	// update the patches
	for (int n=0; n<6; n++) {
		m_patches[n]->LODUpdate(campos);
	}

	return 0;
}

void OceanSphere::Init()
{
	s_patchContext.Reset(new OceanPatchContext(detail_edgeLen[Pi::detail.planets > 4 ? 4 : Pi::detail.planets]));
	assert(s_patchContext->edgeLen <= OCEANPATCH_MAX_EDGELEN);
}

void OceanSphere::Uninit()
{
	assert (s_patchContext.Unique());
	s_patchContext.Reset();
}

void OceanSphere::OnChangeDetailLevel()
{
	s_patchContext.Reset(new OceanPatchContext(detail_edgeLen[Pi::detail.planets > 4 ? 4 : Pi::detail.planets]));
	assert(s_patchContext->edgeLen <= OCEANPATCH_MAX_EDGELEN);
}

#define OCEANSPHERE_TYPE	(m_sbody->type)

OceanSphere::OceanSphere(const SystemBody *body)
{
	//m_terrain = Terrain::InstanceTerrain(body);

	m_vbosToDestroyLock = SDL_CreateMutex();
	m_sbody = body;
	memset(m_patches, 0, 6*sizeof(OceanPatch*));

	//SetUpMaterials is not called until first Render since light count is zero :)
}

OceanSphere::~OceanSphere()
{
	for (int i=0; i<6; i++) {
		if (m_patches[i]) delete m_patches[i];
	}
	DestroyVBOs();
	SDL_DestroyMutex(m_vbosToDestroyLock);

	//delete m_terrain;
}

void OceanSphere::AddVBOToDestroy(GLuint vbo)
{
	SDL_mutexP(m_vbosToDestroyLock);
	m_vbosToDestroy.push_back(vbo);
	SDL_mutexV(m_vbosToDestroyLock);
}

void OceanSphere::DestroyVBOs()
{
	SDL_mutexP(m_vbosToDestroyLock);
	for (std::list<GLuint>::iterator i = m_vbosToDestroy.begin();
			i != m_vbosToDestroy.end(); ++i) {
		glDeleteBuffersARB(1, &(*i));
	}
	m_vbosToDestroy.clear();
	SDL_mutexV(m_vbosToDestroyLock);
}

void OceanSphere::BuildFirstPatches()
{
	// generate initial wank
	vector3d p1(1,1,1);
	vector3d p2(-1,1,1);
	vector3d p3(-1,-1,1);
	vector3d p4(1,-1,1);
	vector3d p5(1,1,-1);
	vector3d p6(-1,1,-1);
	vector3d p7(-1,-1,-1);
	vector3d p8(1,-1,-1);
	p1 = p1.Normalized();
	p2 = p2.Normalized();
	p3 = p3.Normalized();
	p4 = p4.Normalized();
	p5 = p5.Normalized();
	p6 = p6.Normalized();
	p7 = p7.Normalized();
	p8 = p8.Normalized();

	m_patches[0] = new OceanPatch(s_patchContext, this, p1, p2, p3, p4, 0);
	m_patches[1] = new OceanPatch(s_patchContext, this, p4, p3, p7, p8, 0);
	m_patches[2] = new OceanPatch(s_patchContext, this, p1, p4, p8, p5, 0);
	m_patches[3] = new OceanPatch(s_patchContext, this, p2, p1, p5, p6, 0);
	m_patches[4] = new OceanPatch(s_patchContext, this, p3, p2, p6, p7, 0);
	m_patches[5] = new OceanPatch(s_patchContext, this, p8, p7, p6, p5, 0);
	for (int i=0; i<6; i++) {
		for (int j=0; j<4; j++) {
			m_patches[i]->edgeFriend[j] = m_patches[ocean_sphere_edge_friends[i][j]];
		}
	}
	for (int i=0; i<6; i++) m_patches[i]->GenerateMesh();
	for (int i=0; i<6; i++) m_patches[i]->UpdateVBOs();
}

void OceanSphere::DestroyPatches() {
	for (int p=0; p<6; p++) {
		// delete patches
		if (m_patches[p]) {
			delete m_patches[p];
			m_patches[p] = 0;
		}
	}
}

static const float g_ambient[4] = { 0, 0, 0, 1.0 };

void OceanSphere::Render(Graphics::Renderer *renderer, const vector3d campos, const float radius, const float scale) {
	glPushMatrix();
	glTranslated(-campos.x, -campos.y, -campos.z);
	Graphics::Frustum frustum = Graphics::Frustum::FromGLState();

	// no frustum test of entire oceansphere, since Space::Render does this
	// for each body using its GetBoundingRadius() value

	//First draw - create materials (they do not change afterwards)
	if (!m_oceanMaterial.Valid())
		SetUpMaterials();

	if (Graphics::AreShadersEnabled()) {
		matrix4x4d modelMatrix;
		glGetDoublev (GL_MODELVIEW_MATRIX, &modelMatrix[0]);

		//Update material parameters
		//XXX no need to calculate AP every frame
		m_atmosphereParameters = m_sbody->CalcAtmosphereParams();
		m_atmosphereParameters.center = modelMatrix * vector3d(0.0, 0.0, 0.0);
		m_atmosphereParameters.planetRadius = radius;
		m_atmosphereParameters.scale = scale;

		m_oceanMaterial->specialParameter0 = &m_atmosphereParameters;
	}
	glPopMatrix();

	if (!m_patches[0]) BuildFirstPatches();

	Color ambient;
	Color &emission = m_oceanMaterial->emissive;

	// save old global ambient
	const Color oldAmbient = renderer->GetAmbientColor();

	// give planet some ambient lighting if the viewer is close to it
	double camdist = campos.Length();
	camdist = 0.1 / (camdist*camdist);
	// why the fuck is this returning 0.1 when we are sat on the planet??
	// JJ: Because campos is relative to a unit-radius planet - 1.0 at the surface
	// XXX oh well, it is the value we want anyway...
	ambient.r = ambient.g = ambient.b = float(camdist);
	ambient.a = 1.0f;

	renderer->SetAmbientColor(ambient);
	// this is pretty much the only place where a non-renderer is allowed to call Apply()
	// to be removed when someone rewrites terrain
	m_oceanMaterial->Apply();

	for (int i=0; i<6; i++) {
		m_patches[i]->Render(campos, frustum);
	}

	m_oceanMaterial->Unapply();

	renderer->SetAmbientColor(oldAmbient);

	// if the update thread has deleted any oceanpatches, destroy the vbos
	// associated with them
	DestroyVBOs();
}

void OceanSphere::SetUpMaterials()
{
	//SystemBody::TYPE_PLANET_GAS_GIANT
	//SystemBody::TYPE_PLANET_ASTEROID
	assert(m_sbody->type == SystemBody::TYPE_PLANET_TERRESTRIAL);
	// Request material for this star or planet, with or without
	// atmosphere. Separate material for surface and sky.
	Graphics::MaterialDescriptor surfDesc;
	surfDesc.effect = Graphics::EFFECT_OCEANSPHERE_PATCH;
	//planetoid with or without atmosphere
	const SystemBody::AtmosphereParameters ap(m_sbody->CalcAtmosphereParams());
	surfDesc.lighting = true;
	surfDesc.atmosphere = (ap.atmosDensity > 0.0);
	m_oceanMaterial.Reset(Pi::renderer->CreateMaterial(surfDesc));
}
