#include "libs.h"
#include "GeoRing.h"
#include "perlin.h"
#include "Pi.h"
#include "StarSystem.h"
#include "Render.h"

#include "profiler/Profiler.h"

#define GEOPLATE_USE_THREADING

// tri edge lengths
#define GEOPLATE_SUBDIVIDE_AT_CAMDIST	5.0
#define GEOPLATE_MAX_DEPTH	15
// must be an odd number
//#define GEOPLATE_EDGELEN	15
#define GEOPLATE_NUMVERTICES	(GEOPLATE_EDGELEN*GEOPLATE_EDGELEN)
#define GEORING_USE_THREADING

int GEOPLATE_EDGELEN = 15;
static const int GEOPLATE_MAX_EDGELEN = 65;
static double GEOPLATE_FRAC;
static double GEOPLATEHULL_FRAC;

#define lerp(t, a, b) ( a + t * (b - a) )

#define PRINT_VECTOR(_v) printf("%f,%f,%f\n", (_v).x, (_v).y, (_v).z);

SHADER_CLASS_BEGIN(GeoRingShader)
	SHADER_UNIFORM_VEC4(atmosColor)
	SHADER_UNIFORM_FLOAT(geosphereScale)
	SHADER_UNIFORM_FLOAT(geosphereAtmosTopRad)
	SHADER_UNIFORM_VEC3(geosphereCenter)
	SHADER_UNIFORM_FLOAT(geosphereAtmosFogDensity)
SHADER_CLASS_END()

static GeoRingShader *s_geoRingSurfaceShader[4];//, *s_geoRingSkyShader[4];

#pragma pack(4)
struct VBOVertex
{
	float x,y,z;
	float nx,ny,nz;
	unsigned char col[4];
	float padding;
};
#pragma pack()

#define VBO_COUNT_LO_EDGE  (3*(GEOPLATE_EDGELEN/2))
#define VBO_COUNT_HI_EDGE  (3*(GEOPLATE_EDGELEN-1))
#define VBO_COUNT_MID_IDX  (4*3*(GEOPLATE_EDGELEN-3) + 2*(GEOPLATE_EDGELEN-3)*(GEOPLATE_EDGELEN-3)*3)
#define VBO_COUNT_ALL_IDX  (((GEOPLATE_EDGELEN-1)*(GEOPLATE_EDGELEN-1))*2*3)

//                          ^^ serrated teeth bit      ^^^ square inner bit
#define IDX_VBO_LO_OFFSET(_i) ((_i)*sizeof(unsigned short)*3*(GEOPLATE_EDGELEN/2))
#define IDX_VBO_HI_OFFSET(_i) (((_i)*sizeof(unsigned short)*VBO_COUNT_HI_EDGE)+IDX_VBO_LO_OFFSET(4))
#define IDX_VBO_MAIN_OFFSET IDX_VBO_HI_OFFSET(4)

// for glDrawRangeElements
static int s_loMinIdx[4], s_loMaxIdx[4];
static int s_hiMinIdx[4], s_hiMaxIdx[4];

#define NUM_VBE 2

bool rayIntersectsTriangle(
	const vector3d p, 
	const vector3d d,
	const vector3d v0, 
	const vector3d v1, 
	const vector3d v2,
	double &t) 
{
	vector3d e1,e2,h,s,q;
	float a,f,u,v;
	e1 = v1 - v0;
	e2 = v2 - v0;

	h = d.Cross(e2);
	a = e1.Dot(h);

	if (a > -0.00001 && a < 0.00001)
		return false;

	f = 1.0/a;
	s = p - v0;
	u = f * s.Dot(h);

	if (u < 0.0 || u > 1.0)
		return false;

	q = s.Cross(e1);
	v = f * d.Dot(q);

	if (v < 0.0 || u + v > 1.0)
		return false;

	// at this stage we can compute t to find out where
	// the intersection point is on the line
	t = f * e2.Dot(q);

	if (t > 0.00001) {// ray intersection
		return true;
	} else { // this means that there is a line intersection
		return false;// but not a ray intersection
	}
}

class GeoPlate {
public:
	double ang[NUM_VBE];
	double m_halfLen;
	double m_zoffset;		// offset from the centre i.e. newzoffset = (m_zoffset + m_halfLen);
	vector3d *vertices;
	vector3d *normals;
	vector3d *colors;
	GLuint m_vbo;
	static unsigned short *midIndices;
	static unsigned short *loEdgeIndices[4];
	static unsigned short *hiEdgeIndices[4];
	static GLuint indices_vbo;
	static VBOVertex *vbotemp;
	GeoPlate *kids[4];
	GeoPlate *parent;
	GeoPlate *edgeFriend[4]; // [0]=v01, [1]=v12, [2]=v23, [3]=30 - original from GeoPatch...
	GeoRing *geoRing;
	double m_roughLength;
	vector3d clipCentroid;
	double clipRadius;
	int m_depth;
	int m_cIdx;
	SDL_mutex *m_kidsLock;
	bool m_needUpdateVBOs;

	static GeoPlate** s_geoPlates[4];
	static SDL_mutex *s_geoPlateLock[4];
	static SDL_sem* s_geoPlateSem[4];
	static SDL_sem* s_geoRingSem[4];
	static SDL_Thread* s_geoRingThread[4];

	bool RayIntersect(const vector3d& ray, double &minTime, std::vector<double> &timeOI) const
	{
		bool result=false;
		int pointInTriTests=0;
		const vector3d zero(0.0, 0.0, 0.0);
		double time=DBL_MAX;

		for(int y=0; y<GEOPLATE_EDGELEN-1; ++y) 
		{
			for(int x=0; x<GEOPLATE_EDGELEN-1; ++x) 
			{
				// check
				if( rayIntersectsTriangle( zero, ray,
					vertices[x + (y*GEOPLATE_EDGELEN)],
					vertices[(x+1) + (y*GEOPLATE_EDGELEN)],
					vertices[x + ((y+1)*GEOPLATE_EDGELEN)],
					time) ) 
				{
					++pointInTriTests;
					result = true;
					timeOI.push_back( time );
					minTime = std::min( time, minTime );
				}

				// check
				if( rayIntersectsTriangle( zero, ray,
					vertices[(x+1) + (y*GEOPLATE_EDGELEN)],
					vertices[(x+1) + ((y+1)*GEOPLATE_EDGELEN)],
					vertices[x + ((y+1)*GEOPLATE_EDGELEN)],
					time) ) 
				{
					++pointInTriTests;
					result = true;
					timeOI.push_back( time );
					minTime = std::min( time, minTime );
				}
			}
		}
		
		return result;
	}

	double DistFromSurface(const double rad, const double z) const
	{
		PiAssert( IsWithinAng(rad, z) );

		// do we have kids?
		if(kids[0]) {
			// yes, so recurse into them.
			for(int i=0;i<4;++i) {
				if( kids[i]->IsWithinAng(rad, z) ) {
					return kids[i]->DistFromSurface(rad, z);
				}
			}
		} else {
			// find point on _surface_ of the cylinder
			const vector3d p = GetSurfacePointCyl(rad, z, m_halfLen);
			return geoRing->GetHeight(p);
		}
		return DBL_MAX;
	}

	bool IsWithinAng(const double rad, const double z) const
	{
		if( rad>=ang[0] && rad<=ang[1] ) {
			const double zoffsetpos = (m_zoffset + m_halfLen);
			const double zoffsetneg = (m_zoffset - m_halfLen);
			if( z>=zoffsetneg && z<=zoffsetpos ) {
				return true;
			}
		}
		return false;
	}

	/* Thread(s) that generate the mesh data for a geopatch */
	static int UpdateGeoPlateThread(void *data)
	{
		PROFILE_THREAD_SCOPED()

		uint8_t idx = (*(uint8_t*)data);
		PiAssert( idx>=0 && idx<4 );
		PiAssert(s_geoPlateLock[idx]);
		PiAssert(s_geoPlateSem[idx]);
		for(;;) {
			SDL_SemWait(s_geoPlateSem[idx]);
			SDL_mutexP(s_geoPlateLock[idx]);
			if (s_geoPlates[idx]) { 
				(*s_geoPlates[idx])->GenerateMesh();
			}
			s_geoPlates[idx] = NULL;
			SDL_mutexV(s_geoPlateLock[idx]);
			SDL_SemPost(s_geoRingSem[idx]);
		}
		return 0;
	}
	
	// params
	// v0, v1 - define points on the line describing the loop of the ring/orbital.
	// depth - 0 is the topmost plate with each depth+1 describing it's depth within the tree.
	GeoPlate(const double halfLength, const double startAng, const double endAng, 
		const double zoffset, const int depth, const int cIdx) {
		//PROFILE_SCOPED()
		memset(this, 0, sizeof(GeoPlate));
		m_kidsLock = SDL_CreateMutex();
		ang[0] = startAng;
		ang[1] = endAng;
		m_halfLen = halfLength;
		m_zoffset = zoffset;
		m_depth = depth;
		m_cIdx = cIdx;
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
		m_needUpdateVBOs = false;
		normals = new vector3d[GEOPLATE_NUMVERTICES];
		vertices = new vector3d[GEOPLATE_NUMVERTICES];
		colors = new vector3d[GEOPLATE_NUMVERTICES];
#ifdef _DEBUG
		memset(vertices, 0, sizeof(vector3d)*GEOPLATE_NUMVERTICES);
		memset(normals, 0, sizeof(vector3d)*GEOPLATE_NUMVERTICES);
		memset(colors, 0, sizeof(vector3d)*GEOPLATE_NUMVERTICES);
#endif
	}

	static void Init() {
		//PROFILE_SCOPED()
		GEOPLATE_FRAC = 1.0 / double(GEOPLATE_EDGELEN-1);

		if (midIndices) {
			delete [] midIndices;
			for (int i=0; i<4; i++) {
				delete [] loEdgeIndices[i];
				delete [] hiEdgeIndices[i];
			}
			if (indices_vbo) {
				glDeleteBuffersARB(1, &indices_vbo);
			}
			delete [] vbotemp;
		}

		{
			vbotemp = new VBOVertex[GEOPLATE_NUMVERTICES];
				
			unsigned short *idx;
			midIndices = new unsigned short[VBO_COUNT_MID_IDX];
			for (int i=0; i<4; i++) {
				loEdgeIndices[i] = new unsigned short[VBO_COUNT_LO_EDGE];
				hiEdgeIndices[i] = new unsigned short[VBO_COUNT_HI_EDGE];
			}
			/* also want vtx indices for tris not touching edge of patch */
			idx = midIndices;
			for (int x=1; x<GEOPLATE_EDGELEN-2; x++) {
				for (int y=1; y<GEOPLATE_EDGELEN-2; y++) {
					idx[0] = x + GEOPLATE_EDGELEN*y;
					idx[1] = x+1 + GEOPLATE_EDGELEN*y;
					idx[2] = x + GEOPLATE_EDGELEN*(y+1);
					idx+=3;

					idx[0] = x+1 + GEOPLATE_EDGELEN*y;
					idx[1] = x+1 + GEOPLATE_EDGELEN*(y+1);
					idx[2] = x + GEOPLATE_EDGELEN*(y+1);
					idx+=3;
				}
			}
			{
				for (int x=1; x<GEOPLATE_EDGELEN-3; x+=2) {
					// razor teeth near edge 0
					idx[0] = x + GEOPLATE_EDGELEN;
					idx[1] = x+1;
					idx[2] = x+1 + GEOPLATE_EDGELEN;
					idx+=3;
					idx[0] = x+1;
					idx[1] = x+2 + GEOPLATE_EDGELEN;
					idx[2] = x+1 + GEOPLATE_EDGELEN;
					idx+=3;
				}
				for (int x=1; x<GEOPLATE_EDGELEN-3; x+=2) {
					// near edge 2
					idx[0] = x + GEOPLATE_EDGELEN*(GEOPLATE_EDGELEN-2);
					idx[1] = x+1 + GEOPLATE_EDGELEN*(GEOPLATE_EDGELEN-2);
					idx[2] = x+1 + GEOPLATE_EDGELEN*(GEOPLATE_EDGELEN-1);
					idx+=3;
					idx[0] = x+1 + GEOPLATE_EDGELEN*(GEOPLATE_EDGELEN-2);
					idx[1] = x+2 + GEOPLATE_EDGELEN*(GEOPLATE_EDGELEN-2);
					idx[2] = x+1 + GEOPLATE_EDGELEN*(GEOPLATE_EDGELEN-1);
					idx+=3;
				}
				for (int y=1; y<GEOPLATE_EDGELEN-3; y+=2) {
					// near edge 1
					idx[0] = GEOPLATE_EDGELEN-2 + y*GEOPLATE_EDGELEN;
					idx[1] = GEOPLATE_EDGELEN-1 + (y+1)*GEOPLATE_EDGELEN;
					idx[2] = GEOPLATE_EDGELEN-2 + (y+1)*GEOPLATE_EDGELEN;
					idx+=3;
					idx[0] = GEOPLATE_EDGELEN-2 + (y+1)*GEOPLATE_EDGELEN;
					idx[1] = GEOPLATE_EDGELEN-1 + (y+1)*GEOPLATE_EDGELEN;
					idx[2] = GEOPLATE_EDGELEN-2 + (y+2)*GEOPLATE_EDGELEN;
					idx+=3;
				}
				for (int y=1; y<GEOPLATE_EDGELEN-3; y+=2) {
					// near edge 3
					idx[0] = 1 + y*GEOPLATE_EDGELEN;
					idx[1] = 1 + (y+1)*GEOPLATE_EDGELEN;
					idx[2] = (y+1)*GEOPLATE_EDGELEN;
					idx+=3;
					idx[0] = 1 + (y+1)*GEOPLATE_EDGELEN;
					idx[1] = 1 + (y+2)*GEOPLATE_EDGELEN;
					idx[2] = (y+1)*GEOPLATE_EDGELEN;
					idx+=3;
				}
			}
			// full detail edge triangles
			{
				idx = hiEdgeIndices[0];
				for (int x=0; x<GEOPLATE_EDGELEN-1; x+=2) {
					idx[0] = x; idx[1] = x+1; idx[2] = x+1 + GEOPLATE_EDGELEN;
					idx+=3;
					idx[0] = x+1; idx[1] = x+2; idx[2] = x+1 + GEOPLATE_EDGELEN;
					idx+=3;
				}
				idx = hiEdgeIndices[1];
				for (int y=0; y<GEOPLATE_EDGELEN-1; y+=2) {
					idx[0] = GEOPLATE_EDGELEN-1 + y*GEOPLATE_EDGELEN;
					idx[1] = GEOPLATE_EDGELEN-1 + (y+1)*GEOPLATE_EDGELEN;
					idx[2] = GEOPLATE_EDGELEN-2 + (y+1)*GEOPLATE_EDGELEN;
					idx+=3;
					idx[0] = GEOPLATE_EDGELEN-1 + (y+1)*GEOPLATE_EDGELEN;
					idx[1] = GEOPLATE_EDGELEN-1 + (y+2)*GEOPLATE_EDGELEN;
					idx[2] = GEOPLATE_EDGELEN-2 + (y+1)*GEOPLATE_EDGELEN;
					idx+=3;
				}
				idx = hiEdgeIndices[2];
				for (int x=0; x<GEOPLATE_EDGELEN-1; x+=2) {
					idx[0] = x + (GEOPLATE_EDGELEN-1)*GEOPLATE_EDGELEN;
					idx[1] = x+1 + (GEOPLATE_EDGELEN-2)*GEOPLATE_EDGELEN;
					idx[2] = x+1 + (GEOPLATE_EDGELEN-1)*GEOPLATE_EDGELEN;
					idx+=3;
					idx[0] = x+1 + (GEOPLATE_EDGELEN-2)*GEOPLATE_EDGELEN;
					idx[1] = x+2 + (GEOPLATE_EDGELEN-1)*GEOPLATE_EDGELEN;
					idx[2] = x+1 + (GEOPLATE_EDGELEN-1)*GEOPLATE_EDGELEN;
					idx+=3;
				}
				idx = hiEdgeIndices[3];
				for (int y=0; y<GEOPLATE_EDGELEN-1; y+=2) {
					idx[0] = y*GEOPLATE_EDGELEN;
					idx[1] = 1 + (y+1)*GEOPLATE_EDGELEN;
					idx[2] = (y+1)*GEOPLATE_EDGELEN;
					idx+=3;
					idx[0] = (y+1)*GEOPLATE_EDGELEN;
					idx[1] = 1 + (y+1)*GEOPLATE_EDGELEN;
					idx[2] = (y+2)*GEOPLATE_EDGELEN;
					idx+=3;
				}
			}
			// these edge indices are for patches with no
			// neighbour of equal or greater detail -- they reduce
			// their edge complexity by 1 division
			{
				idx = loEdgeIndices[0];
				for (int x=0; x<GEOPLATE_EDGELEN-2; x+=2) {
					idx[0] = x;
					idx[1] = x+2;
					idx[2] = x+1+GEOPLATE_EDGELEN;
					idx += 3;
				}
				idx = loEdgeIndices[1];
				for (int y=0; y<GEOPLATE_EDGELEN-2; y+=2) {
					idx[0] = (GEOPLATE_EDGELEN-1) + y*GEOPLATE_EDGELEN;
					idx[1] = (GEOPLATE_EDGELEN-1) + (y+2)*GEOPLATE_EDGELEN;
					idx[2] = (GEOPLATE_EDGELEN-2) + (y+1)*GEOPLATE_EDGELEN;
					idx += 3;
				}
				idx = loEdgeIndices[2];
				for (int x=0; x<GEOPLATE_EDGELEN-2; x+=2) {
					idx[0] = x+GEOPLATE_EDGELEN*(GEOPLATE_EDGELEN-1);
					idx[2] = x+2+GEOPLATE_EDGELEN*(GEOPLATE_EDGELEN-1);
					idx[1] = x+1+GEOPLATE_EDGELEN*(GEOPLATE_EDGELEN-2);
					idx += 3;
				}
				idx = loEdgeIndices[3];
				for (int y=0; y<GEOPLATE_EDGELEN-2; y+=2) {
					idx[0] = y*GEOPLATE_EDGELEN;
					idx[2] = (y+2)*GEOPLATE_EDGELEN;
					idx[1] = 1 + (y+1)*GEOPLATE_EDGELEN;
					idx += 3;
				}
			}
			// find min/max indices
			for (int i=0; i<4; i++) {
				s_loMinIdx[i] = s_hiMinIdx[i] = 1<<30;
				s_loMaxIdx[i] = s_hiMaxIdx[i] = 0;
				for (int j=0; j<3*(GEOPLATE_EDGELEN/2); j++) {
					if (loEdgeIndices[i][j] < s_loMinIdx[i]) s_loMinIdx[i] = loEdgeIndices[i][j];
					if (loEdgeIndices[i][j] > s_loMaxIdx[i]) s_loMaxIdx[i] = loEdgeIndices[i][j];
				}
				for (int j=0; j<VBO_COUNT_HI_EDGE; j++) {
					if (hiEdgeIndices[i][j] < s_hiMinIdx[i]) s_hiMinIdx[i] = hiEdgeIndices[i][j];
					if (hiEdgeIndices[i][j] > s_hiMaxIdx[i]) s_hiMaxIdx[i] = hiEdgeIndices[i][j];
				}
				//printf("%d:\nLo %d:%d\nHi: %d:%d\n", i, s_loMinIdx[i], s_loMaxIdx[i], s_hiMinIdx[i], s_hiMaxIdx[i]);
			}

			glGenBuffersARB(1, &indices_vbo);
			glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, indices_vbo);
			glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER, IDX_VBO_MAIN_OFFSET + sizeof(unsigned short)*VBO_COUNT_MID_IDX, 0, GL_STATIC_DRAW);
			for (int i=0; i<4; i++) {
				glBufferSubDataARB(GL_ELEMENT_ARRAY_BUFFER, 
					IDX_VBO_LO_OFFSET(i),
					sizeof(unsigned short)*3*(GEOPLATE_EDGELEN/2),
					loEdgeIndices[i]);
			}
			for (int i=0; i<4; i++) {
				glBufferSubDataARB(GL_ELEMENT_ARRAY_BUFFER,
					IDX_VBO_HI_OFFSET(i),
					sizeof(unsigned short)*VBO_COUNT_HI_EDGE,
					hiEdgeIndices[i]);
			}
			glBufferSubDataARB(GL_ELEMENT_ARRAY_BUFFER,
					IDX_VBO_MAIN_OFFSET,
					sizeof(unsigned short)*VBO_COUNT_MID_IDX,
					midIndices);
			glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);

#ifdef GEOPLATE_USE_THREADING
			static const int indexes[4] = {0,1,2,3};
			for( int i = 0; i<4; ++i ) {
				assert(NULL==s_geoPlates[i]);
				s_geoPlates[i] = NULL;

				if(NULL==s_geoPlateLock[i])
					s_geoPlateLock[i] = SDL_CreateMutex();

				if(!s_geoPlateSem[i])
					s_geoPlateSem[i] = SDL_CreateSemaphore(0);

				if(!s_geoRingSem[i])
					s_geoRingSem[i] = SDL_CreateSemaphore(0);

				if(!s_geoRingThread[i])
					s_geoRingThread[i] = SDL_CreateThread(&GeoPlate::UpdateGeoPlateThread, (void*)&indexes[i]);
			}
#endif /* GEOPLATE_USE_THREADING */
		}
	}

	~GeoPlate() {
		SDL_DestroyMutex(m_kidsLock);
		for (int i=0; i<4; i++) {
			if (edgeFriend[i]) edgeFriend[i]->NotifyEdgeFriendDeleted(this);
		}
		for (int i=0; i<4; i++) if (kids[i]) delete kids[i];
		delete vertices;
		delete normals;
		delete colors;
		geoRing->AddVBOToDestroy(m_vbo);
	}
	void UpdateVBOs() {
		m_needUpdateVBOs = true;
	}

	void _UpdateVBOs() {
		//PROFILE_SCOPED()
		if (m_needUpdateVBOs) {
			if (!m_vbo) glGenBuffersARB(1, &m_vbo);
			m_needUpdateVBOs = false;
			glBindBufferARB(GL_ARRAY_BUFFER, m_vbo);
			glBufferDataARB(GL_ARRAY_BUFFER, sizeof(VBOVertex)*GEOPLATE_NUMVERTICES, 0, GL_DYNAMIC_DRAW);
			for (int i=0; i<GEOPLATE_NUMVERTICES; i++)
			{
				clipRadius = std::max(clipRadius, (vertices[i]-clipCentroid).Length());
				VBOVertex *pData = vbotemp + i;
				pData->x = float(vertices[i].x);
				pData->y = float(vertices[i].y);
				pData->z = float(vertices[i].z);
				pData->nx = float(normals[i].x);
				pData->ny = float(normals[i].y);
				pData->nz = float(normals[i].z);
				pData->col[0] = static_cast<unsigned char>(Clamp(colors[i].x*255.0, 0.0, 255.0));
				pData->col[1] = static_cast<unsigned char>(Clamp(colors[i].y*255.0, 0.0, 255.0));
				pData->col[2] = static_cast<unsigned char>(Clamp(colors[i].z*255.0, 0.0, 255.0));
				pData->col[3] = 1.0;
			}
			glBufferDataARB(GL_ARRAY_BUFFER, sizeof(VBOVertex)*GEOPLATE_NUMVERTICES, vbotemp, GL_DYNAMIC_DRAW);
			glBindBufferARB(GL_ARRAY_BUFFER, 0);
		}
	}	
	/* not quite edge, since we share edge vertices so that would be
	 * fucking pointless. one position inwards. used to make edge normals
	 * for adjacent tiles */
	void GetEdgeMinusOneVerticesFlipped(int edge, vector3d *ev) {
		//PROFILE_SCOPED()
		if( edge<0 || edge>4) return;
		if (edge == 0) {
			for (int x=0; x<GEOPLATE_EDGELEN; x++) ev[GEOPLATE_EDGELEN-1-x] = vertices[x + GEOPLATE_EDGELEN];
		} else if (edge == 1) {
			const int x = GEOPLATE_EDGELEN-2;
			for (int y=0; y<GEOPLATE_EDGELEN; y++) ev[GEOPLATE_EDGELEN-1-y] = vertices[x + y*GEOPLATE_EDGELEN];
		} else if (edge == 2) {
			const int y = GEOPLATE_EDGELEN-2;
			for (int x=0; x<GEOPLATE_EDGELEN; x++) ev[GEOPLATE_EDGELEN-1-x] = vertices[(GEOPLATE_EDGELEN-1)-x + y*GEOPLATE_EDGELEN];
		} else {
			for (int y=0; y<GEOPLATE_EDGELEN; y++) ev[GEOPLATE_EDGELEN-1-y] = vertices[1 + ((GEOPLATE_EDGELEN-1)-y)*GEOPLATE_EDGELEN];
		}
	}
	static void GetEdge(vector3d *verts, int edge, vector3d *ev) {
		//PROFILE_SCOPED()
		if( edge<0 || edge>4) return;
		if (edge == 0) {
			for (int x=0; x<GEOPLATE_EDGELEN; x++) ev[x] = verts[x];
		} else if (edge == 1) {
			const int x = GEOPLATE_EDGELEN-1;
			for (int y=0; y<GEOPLATE_EDGELEN; y++) ev[y] = verts[x + y*GEOPLATE_EDGELEN];
		} else if (edge == 2) {
			const int y = GEOPLATE_EDGELEN-1;
			for (int x=0; x<GEOPLATE_EDGELEN; x++) ev[x] = verts[(GEOPLATE_EDGELEN-1)-x + y*GEOPLATE_EDGELEN];
		} else {
			for (int y=0; y<GEOPLATE_EDGELEN; y++) ev[y] = verts[0 + ((GEOPLATE_EDGELEN-1)-y)*GEOPLATE_EDGELEN];
		}
	}
	static void SetEdge(vector3d *verts, int edge, const vector3d *ev) {
		//PROFILE_SCOPED()
		if( edge<0 || edge>4) return;
		if (edge == 0) {
			for (int x=0; x<GEOPLATE_EDGELEN; x++) verts[x] = ev[x];
		} else if (edge == 1) {
			const int x = GEOPLATE_EDGELEN-1;
			for (int y=0; y<GEOPLATE_EDGELEN; y++) verts[x + y*GEOPLATE_EDGELEN] = ev[y];
		} else if (edge == 2) {
			const int y = GEOPLATE_EDGELEN-1;
			for (int x=0; x<GEOPLATE_EDGELEN; x++) verts[(GEOPLATE_EDGELEN-1)-x + y*GEOPLATE_EDGELEN] = ev[x];
		} else {
			for (int y=0; y<GEOPLATE_EDGELEN; y++) verts[0 + ((GEOPLATE_EDGELEN-1)-y)*GEOPLATE_EDGELEN] = ev[y];
		}
	}
	int GetEdgeIdxOf(GeoPlate *e) {
		//PROFILE_SCOPED()
		for (int i=0; i<4; i++) {
			if (edgeFriend[i] == e) return i;
		}
		//abort();
		return -1;
	}


	void FixEdgeNormals(const int edge, const vector3d *ev) {
		//PROFILE_SCOPED()
		int x, y;
		switch (edge) {
		case 0:
			for (x=1; x<GEOPLATE_EDGELEN-1; x++) {
				const vector3d x1 = vertices[x-1];
				const vector3d x2 = vertices[x+1];
				const vector3d y1 = ev[x];
				const vector3d y2 = vertices[x + GEOPLATE_EDGELEN];
				const vector3d norm = (x2-x1).Cross(y2-y1).Normalized();
				normals[x] = norm;
				// make color
				const vector3d p = GetSurfacePointCyl(x*GEOPLATE_FRAC, 0, m_halfLen);
				const double height = colors[x].x;
				colors[x] = geoRing->GetColor(p, height, norm);
			}
			break;
		case 1:
			x = GEOPLATE_EDGELEN-1;
			for (y=1; y<GEOPLATE_EDGELEN-1; y++) {
				const vector3d x1 = vertices[(x-1) + y*GEOPLATE_EDGELEN];
				const vector3d x2 = ev[y];
				const vector3d y1 = vertices[x + (y-1)*GEOPLATE_EDGELEN];
				const vector3d y2 = vertices[x + (y+1)*GEOPLATE_EDGELEN];
				const vector3d norm = (x2-x1).Cross(y2-y1).Normalized();
				normals[x + y*GEOPLATE_EDGELEN] = norm;
				// make color
				const vector3d p = GetSurfacePointCyl(x*GEOPLATE_FRAC, y*GEOPLATE_FRAC, m_halfLen);
				const double height = colors[x + y*GEOPLATE_EDGELEN].x;
				colors[x + y*GEOPLATE_EDGELEN] = geoRing->GetColor(p, height, norm);
			}
			break;
		case 2:
			y = GEOPLATE_EDGELEN-1;
			for (x=1; x<GEOPLATE_EDGELEN-1; x++) {
				const vector3d x1 = vertices[x-1 + y*GEOPLATE_EDGELEN];
				const vector3d x2 = vertices[x+1 + y*GEOPLATE_EDGELEN];
				const vector3d y1 = vertices[x + (y-1)*GEOPLATE_EDGELEN];
				const vector3d y2 = ev[GEOPLATE_EDGELEN-1-x];
				const vector3d norm = (x2-x1).Cross(y2-y1).Normalized();
				normals[x + y*GEOPLATE_EDGELEN] = norm;
				// make color
				const vector3d p = GetSurfacePointCyl(x*GEOPLATE_FRAC, y*GEOPLATE_FRAC, m_halfLen);
				const double height = colors[x + y*GEOPLATE_EDGELEN].x;
				colors[x + y*GEOPLATE_EDGELEN] = geoRing->GetColor(p, height, norm);
			}
			break;
		case 3:
			for (y=1; y<GEOPLATE_EDGELEN-1; y++) {
				const vector3d x1 = ev[GEOPLATE_EDGELEN-1-y];
				const vector3d x2 = vertices[1 + y*GEOPLATE_EDGELEN];
				const vector3d y1 = vertices[(y-1)*GEOPLATE_EDGELEN];
				const vector3d y2 = vertices[(y+1)*GEOPLATE_EDGELEN];
				const vector3d norm = (x2-x1).Cross(y2-y1).Normalized();
				normals[y*GEOPLATE_EDGELEN] = norm;
				// make color
				const vector3d p = GetSurfacePointCyl(0, y*GEOPLATE_FRAC, m_halfLen);
				const double height = colors[y*GEOPLATE_EDGELEN].x;
				colors[y*GEOPLATE_EDGELEN] = geoRing->GetColor(p, height, norm);
			}
			break;
		}
	}

	int GetChildIdx(GeoPlate *child) {
		//PROFILE_SCOPED()
		for (int i=0; i<4; i++) {
			if (kids[i] == child) return i;
		}
		abort();
		return -1;
	}
	
	void FixEdgeFromParentInterpolated(int edge) {
		//PROFILE_SCOPED()
		// noticeable artefacts from not doing so...
		vector3d ev[GEOPLATE_MAX_EDGELEN];
		vector3d en[GEOPLATE_MAX_EDGELEN];
		vector3d ec[GEOPLATE_MAX_EDGELEN];
		vector3d ev2[GEOPLATE_MAX_EDGELEN];
		vector3d en2[GEOPLATE_MAX_EDGELEN];
		vector3d ec2[GEOPLATE_MAX_EDGELEN];
		GetEdge(parent->vertices, edge, ev);
		GetEdge(parent->normals, edge, en);
		GetEdge(parent->colors, edge, ec);

		int kid_idx = parent->GetChildIdx(this);
		if (edge == kid_idx) {
			// use first half of edge
			for (int i=0; i<=GEOPLATE_EDGELEN/2; i++) {
				ev2[i<<1] = ev[i];
				en2[i<<1] = en[i];
				ec2[i<<1] = ec[i];
			}
		} else {
			// use 2nd half of edge
			for (int i=GEOPLATE_EDGELEN/2; i<GEOPLATE_EDGELEN; i++) {
				ev2[(i-(GEOPLATE_EDGELEN/2))<<1] = ev[i];
				en2[(i-(GEOPLATE_EDGELEN/2))<<1] = en[i];
				ec2[(i-(GEOPLATE_EDGELEN/2))<<1] = ec[i];
			}
		}
		// interpolate!!
		for (int i=1; i<GEOPLATE_EDGELEN; i+=2) {
			ev2[i] = (ev2[i-1]+ev2[i+1]) * 0.5;
			en2[i] = (en2[i-1]+en2[i+1]).Normalized();
			ec2[i] = (ec2[i-1]+ec2[i+1]) * 0.5;
		}
		
		// hacking
		/*vector3d colorIdx[4] = { 
			vector3d( 1.0, 0.0, 0.0 ),
			vector3d( 0.0, 1.0, 0.0 ),
			vector3d( 0.0, 0.0, 1.0 ),
			vector3d( 1.0, 0.0, 1.0 )
		};
		for (int i=1; i<GEOPLATE_EDGELEN; ++i) {
			if( edge>=0 && edge<4) {
				ec2[i] = colorIdx[edge];
			}
		}*/
		SetEdge(this->vertices, edge, ev2);
		SetEdge(this->normals, edge, en2);
		SetEdge(this->colors, edge, ec2);
	}

	template <int corner>
	void MakeCornerNormal(vector3d *ev, vector3d *ev2) {
		//PROFILE_SCOPED()
		int p;
		vector3d x1,x2,y1,y2;
		switch (corner) {
		case 0: {
			x1 = ev[GEOPLATE_EDGELEN-1];
			x2 = vertices[1];
			y1 = ev2[0];
			y2 = vertices[GEOPLATE_EDGELEN];
			const vector3d norm = (x2-x1).Cross(y2-y1).Normalized();
			normals[0] = norm;
			// make color
			const vector3d pt = GetSurfacePointCyl(0, 0, m_halfLen);
		//	const double height = colors[0].x;
			const double height = geoRing->GetHeight(pt);
			colors[0] = geoRing->GetColor(pt, height, norm);
			}
			break;
		case 1: {
			p = GEOPLATE_EDGELEN-1;
			x1 = vertices[p-1];
			x2 = ev2[0];
			y1 = ev[GEOPLATE_EDGELEN-1];
			y2 = vertices[p + GEOPLATE_EDGELEN];
			const vector3d norm = (x2-x1).Cross(y2-y1).Normalized();
			normals[p] = norm;
			// make color
			const vector3d pt = GetSurfacePointCyl(p*GEOPLATE_FRAC, 0, m_halfLen);
		//	const double height = colors[p].x;
			const double height = geoRing->GetHeight(pt);
			colors[p] = geoRing->GetColor(pt, height, norm);
			}
			break;
		case 2: {
			p = GEOPLATE_EDGELEN-1;
			x1 = vertices[(p-1) + p*GEOPLATE_EDGELEN];
			x2 = ev[GEOPLATE_EDGELEN-1];
			y1 = vertices[p + (p-1)*GEOPLATE_EDGELEN];
			y2 = ev2[0];
			const vector3d norm = (x2-x1).Cross(y2-y1).Normalized();
			normals[p + p*GEOPLATE_EDGELEN] = norm;
			// make color
			const vector3d pt = GetSurfacePointCyl(p*GEOPLATE_FRAC, p*GEOPLATE_FRAC, m_halfLen);
		//	const double height = colors[p + p*GEOPLATE_EDGELEN].x;
			const double height = geoRing->GetHeight(pt);
			colors[p + p*GEOPLATE_EDGELEN] = geoRing->GetColor(pt, height, norm);
			}
			break;
		case 3: {
			p = GEOPLATE_EDGELEN-1;
			x1 = ev2[0];
			x2 = vertices[1 + p*GEOPLATE_EDGELEN];
			y1 = vertices[(p-1)*GEOPLATE_EDGELEN];
			y2 = ev[GEOPLATE_EDGELEN-1];
			const vector3d norm = (x2-x1).Cross(y2-y1).Normalized();
			normals[p*GEOPLATE_EDGELEN] = norm;
			// make color
			const vector3d pt = GetSurfacePointCyl(0, p*GEOPLATE_FRAC, m_halfLen);
		//	const double height = colors[p*GEOPLATE_EDGELEN].x;
			const double height = geoRing->GetHeight(pt);
			colors[p*GEOPLATE_EDGELEN] = geoRing->GetColor(pt, height, norm);
			}
			break;
		}
	}

	void FixCornerNormalsByEdge(int edge, vector3d *ev) {
		//PROFILE_SCOPED()
		vector3d ev2[GEOPLATE_MAX_EDGELEN];
		vector3d x1, x2, y1, y2;
		/* XXX All these 'if's have an unfinished else, when a neighbour
		 * of our size doesn't exist and instead we must look at a bigger tile.
		 * But let's just leave it for the mo because it is a pain.
		 * See comment in OnEdgeFriendChanged() */
		switch (edge) {
		case 0:
			if (edgeFriend[3]) {
				int we_are = edgeFriend[3]->GetEdgeIdxOf(this);
				edgeFriend[3]->GetEdgeMinusOneVerticesFlipped(we_are, ev2);
				MakeCornerNormal<0>(ev2, ev);
			}
			if (edgeFriend[1]) {
				int we_are = edgeFriend[1]->GetEdgeIdxOf(this);
				edgeFriend[1]->GetEdgeMinusOneVerticesFlipped(we_are, ev2);
				MakeCornerNormal<1>(ev, ev2);
			}
			break;
		case 1:
			if (edgeFriend[0]) {
				int we_are = edgeFriend[0]->GetEdgeIdxOf(this);
				edgeFriend[0]->GetEdgeMinusOneVerticesFlipped(we_are, ev2);
				MakeCornerNormal<1>(ev2, ev);
			}
			if (edgeFriend[2]) {
				int we_are = edgeFriend[2]->GetEdgeIdxOf(this);
				edgeFriend[2]->GetEdgeMinusOneVerticesFlipped(we_are, ev2);
				MakeCornerNormal<2>(ev, ev2);
			}
			break;
		case 2:
			if (edgeFriend[1]) {
				int we_are = edgeFriend[1]->GetEdgeIdxOf(this);
				edgeFriend[1]->GetEdgeMinusOneVerticesFlipped(we_are, ev2);
				MakeCornerNormal<2>(ev2, ev);
			}
			if (edgeFriend[3]) {
				int we_are = edgeFriend[3]->GetEdgeIdxOf(this);
				edgeFriend[3]->GetEdgeMinusOneVerticesFlipped(we_are, ev2);
				MakeCornerNormal<3>(ev, ev2);
			}
			break;
		case 3:
			if (edgeFriend[2]) {
				int we_are = edgeFriend[2]->GetEdgeIdxOf(this);
				edgeFriend[2]->GetEdgeMinusOneVerticesFlipped(we_are, ev2);
				MakeCornerNormal<3>(ev2, ev);
			}
			if (edgeFriend[0]) {
				int we_are = edgeFriend[0]->GetEdgeIdxOf(this);
				edgeFriend[0]->GetEdgeMinusOneVerticesFlipped(we_are, ev2);
				MakeCornerNormal<0>(ev, ev2);
			}
			break;
		}
				
	}

	void GenerateEdgeNormalsAndColors() {
		//PROFILE_SCOPED()
		vector3d ev[4][GEOPLATE_MAX_EDGELEN];
		bool doneEdge[4];
		memset(doneEdge, 0, sizeof(doneEdge));
		for (int i=0; i<4; i++) {
			GeoPlate *e = edgeFriend[i];
			if (e) {
				int we_are = e->GetEdgeIdxOf(this);
				e->GetEdgeMinusOneVerticesFlipped(we_are, ev[i]);
			} else if( parent && parent->edgeFriend[i] ) {
				assert(parent->edgeFriend[i]);
				doneEdge[i] = true;
				// parent has valid edge, so take our
				// bit of that, interpolated.
				FixEdgeFromParentInterpolated(i);
				// XXX needed for corners... probably not
				// correct
				GetEdge(vertices, i, ev[i]);
			}
		}
	
		MakeCornerNormal<0>(ev[3], ev[0]);
		MakeCornerNormal<1>(ev[0], ev[1]);
		MakeCornerNormal<2>(ev[1], ev[2]);
		MakeCornerNormal<3>(ev[2], ev[3]);

		for (int i=0; i<4; i++) {
			if(!doneEdge[i]) {
				FixEdgeNormals(i, ev[i]);
			}
		}
	}

	vector3d GetSurfacePointCyl(const double x, const double y, const double halfLength) const {
		//PiAssert(y>=-0.000001);
		double theta = lerp( x, ang[1], ang[0] );//*2.0*M_PI;
		
		const vector3d topEndEdge(cos(theta), sin(theta), m_zoffset + halfLength);		// vertices at top edge of circle
		const vector3d bottomEndEdge(cos(theta), sin(theta), m_zoffset - halfLength);	// vertices at bottom edge of circle
		
		const vector3d res = lerp( y, bottomEndEdge, topEndEdge );
		//PiAssert(res.z>bottomEndEdge.z);
		return res;
	}

	// Generates full-detail vertices, and also non-edge normals and colors
	void GenerateMesh() {
		//PROFILE_SCOPED()
		vector3d *vts = vertices;
		vector3d *col = colors;
		double xfrac;
		double yfrac = 0;
		vector3d axisPt(0.0, 0.0, 0.0);

		const vector3d topEnd(0.0, 0.0, m_zoffset + m_halfLen);		// vertices at top edge of circle
		const vector3d btmEnd(0.0, 0.0, m_zoffset - m_halfLen);	// vertices at bottom edge of circle
				
		for (int y=0; y<GEOPLATE_EDGELEN; ++y) {
			xfrac = 0;
			// point along z-axis by yfrac amount
			const vector3d axisPt = lerp( yfrac, btmEnd, topEnd );
			for (int x=0; x<GEOPLATE_EDGELEN; ++x) {
				// find point on _surface_ of the cylinder
				const vector3d p = GetSurfacePointCyl(xfrac, yfrac, m_halfLen);
				// vector from axis to point-on-surface
				const vector3d cDir = (p - axisPt).Normalized();
				// height
				double height = -geoRing->GetHeight(p);
				// vertex is moved in direction of point-in-axis FROM point-on-surface by height.
				*(vts++) = p + (cDir * height);

				// remember this -- we will need it later
				(col++)->x = -height;

				xfrac += GEOPLATE_FRAC;

				normals[x + y*GEOPLATE_EDGELEN] = (cDir);
			}
			yfrac += GEOPLATE_FRAC;
		}
		assert(vts == &vertices[GEOPLATE_NUMVERTICES]);
		// Generate normals & colors for non-edge vertices since they never change
		for (int y=1; y<GEOPLATE_EDGELEN-1; y++) {
			for (int x=1; x<GEOPLATE_EDGELEN-1; x++) {
				// normal
				vector3d x1 = vertices[x-1 + y*GEOPLATE_EDGELEN];
				vector3d x2 = vertices[x+1 + y*GEOPLATE_EDGELEN];
				vector3d y1 = vertices[x + (y-1)*GEOPLATE_EDGELEN];
				vector3d y2 = vertices[x + (y+1)*GEOPLATE_EDGELEN];

				vector3d n = (x2-x1).Cross(y2-y1);
				const vector3d &norm = normals[x + y*GEOPLATE_EDGELEN] = (n.Normalized());
				// color
				vector3d p = GetSurfacePointCyl(x*GEOPLATE_FRAC, y*GEOPLATE_FRAC, m_halfLen);
				vector3d &col_r = colors[x + y*GEOPLATE_EDGELEN];
				const double height = col_r.x;
				col_r = geoRing->GetColor(p, height, norm);
			}
		}
	}

	void OnEdgeFriendChanged(int edge, GeoPlate *e) {
		//PROFILE_SCOPED()
		assert(e>(void*)0x000000FF);
		edgeFriend[edge] = e;
		vector3d ev[GEOPLATE_MAX_EDGELEN];
		int we_are = e->GetEdgeIdxOf(this);
		e->GetEdgeMinusOneVerticesFlipped(we_are, ev);
		
		vector3d axisPt(0.0, 0.0, 0.0);
		const vector3d topEnd(0.0, 0.0, m_zoffset + m_halfLen);		// vertices at top edge of circle
		const vector3d btmEnd(0.0, 0.0, m_zoffset - m_halfLen);	// vertices at bottom edge of circle

		/* now we have a valid edge, fix the edge vertices */
		if (edge == 0) {
			// point along z-axis by yfrac amount
			const vector3d axisPt = lerp( 0.0, btmEnd, topEnd );
			for (int x=0; x<GEOPLATE_EDGELEN; x++) {
				vector3d p = GetSurfacePointCyl(x * GEOPLATE_FRAC, 0, m_halfLen);
				// vector from axis to point-on-surface
				const vector3d cDir = (p - axisPt).Normalized();
				double height = -geoRing->GetHeight(p);
				vertices[x] = p + (cDir * height);
				// XXX These bounds checks in each edge case are
				// only necessary while the "All these 'if's"
				// comment in FixCOrnerNormalsByEdge stands
				if ((x>0) && (x<GEOPLATE_EDGELEN-1)) {
					colors[x].x = -height;
				}
			}
		} else if (edge == 1) {
			for (int y=0; y<GEOPLATE_EDGELEN; y++) {
				const double yfrac = y * GEOPLATE_FRAC;
				// point along z-axis by yfrac amount
				const vector3d axisPt = lerp( yfrac, btmEnd, topEnd );
				vector3d p = GetSurfacePointCyl(1.0, yfrac, m_halfLen);
				double height = -geoRing->GetHeight(p);
				// vector from axis to point-on-surface
				const vector3d cDir = (p - axisPt).Normalized();
				// vertex is moved in direction of point-in-axis FROM point-on-surface by height.
				int pos = (GEOPLATE_EDGELEN-1) + y*GEOPLATE_EDGELEN;
				vertices[pos] = p + (cDir * height);
				if ((y>0) && (y<GEOPLATE_EDGELEN-1)) {
					colors[pos].x = -height;
				}
			}
		} else if (edge == 2) {
			// point along z-axis by yfrac amount
			const vector3d axisPt = lerp( 1.0, btmEnd, topEnd );
			for (int x=0; x<GEOPLATE_EDGELEN; x++) {
				vector3d p = GetSurfacePointCyl(x * GEOPLATE_FRAC, 1.0, m_halfLen);
				// vector from axis to point-on-surface
				const vector3d cDir = (p - axisPt).Normalized();
				double height = -geoRing->GetHeight(p);
				int pos = x + (GEOPLATE_EDGELEN-1)*GEOPLATE_EDGELEN;
				vertices[pos] = p + (cDir * height);
				if ((x>0) && (x<GEOPLATE_EDGELEN-1)) {
					colors[pos].x = -height;
				}
			}
		} else {
			for (int y=0; y<GEOPLATE_EDGELEN; y++) {
				const double yfrac = y * GEOPLATE_FRAC;
				// point along z-axis by yfrac amount
				const vector3d axisPt = lerp( yfrac, btmEnd, topEnd );
				vector3d p = GetSurfacePointCyl(0, yfrac, m_halfLen);
				double height = -geoRing->GetHeight(p);
				// vector from axis to point-on-surface
				const vector3d cDir = (p - axisPt).Normalized();
				int pos = y * GEOPLATE_EDGELEN;
				vertices[pos] = p + (cDir * height);
				if ((y>0) && (y<GEOPLATE_EDGELEN-1)) {
					colors[pos].x = -height;
				}
			}
		}

		FixEdgeNormals(edge, ev);
		FixCornerNormalsByEdge(edge, ev);
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
	void NotifyEdgeFriendSplit(GeoPlate *e) {
		//PROFILE_SCOPED()
		int idx = GetEdgeIdxOf(e);
		int we_are = e->GetEdgeIdxOf(this);
		if (!kids[0]) return;
		if (-1 == we_are || -1 == idx) return;
		// match e's new kids to our own... :/
		kids[idx]->OnEdgeFriendChanged(idx, e->kids[(we_are+1)%4]);
		kids[(idx+1)%4]->OnEdgeFriendChanged(idx, e->kids[we_are]);
	}
	
	void NotifyEdgeFriendDeleted(GeoPlate *e) {
		//PROFILE_SCOPED()
		int idx = GetEdgeIdxOf(e);
		if (-1 == idx) return;
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

	GeoPlate *GetEdgeFriendForKid(int kid, int edge) {
		//PROFILE_SCOPED()
		GeoPlate *e = edgeFriend[edge];
		if (!e) return 0;
		//assert (e);// && (e->m_depth >= m_depth));
		const int we_are = e->GetEdgeIdxOf(this);
		if (-1 == we_are ) return 0;
		// neighbour patch has not split yet (is at depth of this patch), so kids of this patch do
		// not have same detail level neighbours yet
		if (edge == kid) return e->kids[(we_are+1)%4];
		else return e->kids[we_are];
	}
	
	void Render(vector3d &campos, Plane planes[6]) {
		//PROFILE_SCOPED()
		PiVerify(SDL_mutexP(m_kidsLock)==0);
		if (kids[0]) {
			for (int i=0; i<4; i++) kids[i]->Render(campos, planes);
			SDL_mutexV(m_kidsLock);
		} else {
			SDL_mutexV(m_kidsLock);
			_UpdateVBOs();
			/* frustum test! */
			for (int i=0; i<6; i++) {
				if (planes[i].DistanceToPoint(clipCentroid) <= -clipRadius) {
					return;
				}
			}
			Pi::statSceneTris += 2*(GEOPLATE_EDGELEN-1)*(GEOPLATE_EDGELEN-1);
			glEnableClientState(GL_VERTEX_ARRAY);
			glEnableClientState(GL_NORMAL_ARRAY);
			glEnableClientState(GL_COLOR_ARRAY);

			glBindBufferARB(GL_ARRAY_BUFFER, m_vbo);
			glVertexPointer(3, GL_FLOAT, sizeof(VBOVertex), 0);
			glNormalPointer(GL_FLOAT, sizeof(VBOVertex), reinterpret_cast<void *>(3*sizeof(float)));
			glColorPointer(4, GL_UNSIGNED_BYTE, sizeof(VBOVertex), reinterpret_cast<void *>(6*sizeof(float)));
			glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, indices_vbo);
			glDrawRangeElements(GL_TRIANGLES, 0, GEOPLATE_NUMVERTICES-1, VBO_COUNT_MID_IDX, GL_UNSIGNED_SHORT, reinterpret_cast<void*>(IDX_VBO_MAIN_OFFSET));
			for (int i=0; i<4; i++) {
				if (edgeFriend[i]) {
					glDrawRangeElements(GL_TRIANGLES, s_hiMinIdx[i], s_hiMaxIdx[i], VBO_COUNT_HI_EDGE, GL_UNSIGNED_SHORT, reinterpret_cast<void*>(IDX_VBO_HI_OFFSET(i)));
				} else {
					glDrawRangeElements(GL_TRIANGLES, s_loMinIdx[i], s_loMaxIdx[i], VBO_COUNT_LO_EDGE, GL_UNSIGNED_SHORT, reinterpret_cast<void*>(IDX_VBO_LO_OFFSET(i)));
				}
			}
			glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
			glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);

			glDisableClientState(GL_VERTEX_ARRAY);
			glDisableClientState(GL_NORMAL_ARRAY);
			glDisableClientState(GL_COLOR_ARRAY);

			/*static const vector3d colorIdx[4] = { 
				vector3d( 1.0, 0.0, 0.0 ),	// red
				vector3d( 0.0, 1.0, 0.0 ),	// green
				vector3d( 0.0, 0.0, 1.0 ),	// blue
				vector3d( 1.0, 0.0, 1.0 )	// purple
			};
			glLineWidth(1.0);
			for (int i=0; i<4; ++i) {
				if (edgeFriend[i]) {
					glBegin(GL_LINES);
					glColor3dv( &colorIdx[i][0] );
					glVertex3dv( &clipCentroid[0] );
					vector3d dir = (edgeFriend[i]->clipCentroid - clipCentroid).Normalized()*0.000075;
					glVertex3dv( &(clipCentroid + dir)[0] );
					glEnd();
				}
			}

			glLineWidth(5.0);
			glBegin(GL_LINES);
			//glColor3f( 1.0f, 1.0f, 1.0f );
			glColor3dv( &colorIdx[m_cIdx][0] );
			static const vector3d right(0.00001, 0.0, 0.0);
			static const vector3d up(0.0, 0.00001, 0.0);
			static const vector3d fwd(0.0, 0.0, 0.00001);
			// x axis
			glVertex3dv( &(clipCentroid-right)[0] );
			glVertex3dv( &(clipCentroid+right)[0] );
			// y axis
			glVertex3dv( &(clipCentroid-up)[0] );
			glVertex3dv( &(clipCentroid+up)[0] );
			// z axis
			glVertex3dv( &(clipCentroid-fwd)[0] );
			glVertex3dv( &(clipCentroid+fwd)[0] );
			glEnd();*/
		}
	}

	void LODUpdate(vector3d &campos) {
		//vector3d centroid = (vBE[0]+vBE[1]).Normalized();
		vector3d centroid = GetSurfacePointCyl(0.5, 0.5, m_halfLen);
		centroid = (1.0 + (-geoRing->GetHeight(centroid))) * centroid;

		bool canSplit = true;
		int nullFriends = 0;
		for (int i=0; i<4; i++) {
			if (!edgeFriend[i]) { 
				if( (++nullFriends)>1 && m_depth>0 ) {
					canSplit = false; 
					break; 
				}
			}
			if (edgeFriend[i] && (edgeFriend[i]->m_depth < m_depth)) {
				canSplit = false;
				break;
			}
		}
		if (!(canSplit && (m_depth < GEOPLATE_MAX_DEPTH) &&
		    ((campos - centroid).Length() < m_roughLength)))
			canSplit = false;
		// always split at first level
		if (!parent) canSplit = true;
		/*canSplit = false;*/

		bool canMerge = true;

		if (canSplit) {
			if (!kids[0]) {
				const double halfAng = ang[0] + 0.5 * (ang[1] - ang[0]);
				const double newLength = (m_halfLen*0.5);
				const double zoffset0 = m_zoffset + newLength;
				const double zoffset1 = m_zoffset - newLength;
				GeoPlate *_kids[4];
				_kids[0] = new GeoPlate(newLength, halfAng, ang[1], zoffset1, m_depth+1, 0);
				_kids[1] = new GeoPlate(newLength, ang[0], halfAng, zoffset1, m_depth+1, 1);
				_kids[2] = new GeoPlate(newLength, ang[0], halfAng, zoffset0, m_depth+1, 2);
				_kids[3] = new GeoPlate(newLength, halfAng, ang[1], zoffset0, m_depth+1, 3);
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
				_kids[0]->geoRing = _kids[1]->geoRing = _kids[2]->geoRing = _kids[3]->geoRing = geoRing;

				for(int i=0; i<4; ++i) {
					for(int j=0; j<4; ++j) {
						if(_kids[i]->edgeFriend[j]) {
							assert(_kids[i]->edgeFriend[j] > (void*)0x000000FF);
						}
					}
				}

#ifdef GEOPLATE_USE_THREADING
				for (int i=0; i<4; ++i) {
					PiAssert(GeoPlate::s_geoPlates[i]==0);
					GeoPlate::s_geoPlates[i] = &_kids[i];
				}
				for (int i=0; i<4; ++i) {
					PiAssert(SDL_SemValue(GeoPlate::s_geoPlateSem[i])==0);
					SDL_SemPost(GeoPlate::s_geoPlateSem[i]);
				}

				// not a good enough waiting method
				SDL_SemWait(s_geoRingSem[0]);
				SDL_SemWait(s_geoRingSem[1]);
				SDL_SemWait(s_geoRingSem[2]);
				SDL_SemWait(s_geoRingSem[3]);
#ifdef _DEBUG
				for (int i=0; i<4; ++i) {
					PiAssert(GeoPlate::s_geoPlates[i]==0);
					PiAssert(SDL_SemValue(GeoPlate::s_geoPlateSem[i])==0);
				}
#endif
#else
				for (int i=0; i<4; ++i) _kids[i]->GenerateMesh();
#endif
				PiVerify(SDL_mutexP(m_kidsLock)==0);
				for (int i=0; i<4; ++i) kids[i] = _kids[i];
				for (int i=0; i<4; ++i) {
					if(edgeFriend[i]) edgeFriend[i]->NotifyEdgeFriendSplit(this);
				}
				for (int i=0; i<4; ++i) {
					kids[i]->GenerateEdgeNormalsAndColors();
					kids[i]->UpdateVBOs();
				}
				PiVerify(SDL_mutexV(m_kidsLock)!=-1);
			}
			for (int i=0; i<4; ++i) kids[i]->LODUpdate(campos);
		} else {
			if (canMerge && kids[0]) {
				PiVerify(SDL_mutexP(m_kidsLock)==0);
				for (int i=0; i<4; ++i) { delete kids[i]; kids[i] = 0; }
				PiVerify(SDL_mutexV(m_kidsLock)!=-1);
			}
		}
	}
};
//static 
unsigned short *GeoPlate::midIndices = 0;
unsigned short *GeoPlate::loEdgeIndices[4];
unsigned short *GeoPlate::hiEdgeIndices[4];
GLuint GeoPlate::indices_vbo;
VBOVertex *GeoPlate::vbotemp;

#ifdef GEOPLATE_USE_THREADING
GeoPlate**	GeoPlate::s_geoPlates[4]	= {0};
SDL_mutex*	GeoPlate::s_geoPlateLock[4] = {0};
SDL_sem*	GeoPlate::s_geoPlateSem[4]	= {0};
SDL_sem*	GeoPlate::s_geoRingSem[4]	= {0};
SDL_Thread*	GeoPlate::s_geoRingThread[4]= {0};
#endif // GEOPLATE_USE_THREADING

#define PLANET_AMBIENT	0.0
static std::list<GeoRing*> s_allGeoRings;
SDL_mutex *s_allGeoRingsLock;

/* Thread that updates geoRing level of detail thingies */
int GeoRing::UpdateLODThread(void *data)
{
	PROFILE_THREAD_SCOPED()
	for(;;) {
		SDL_mutexP(s_allGeoRingsLock);
		for(std::list<GeoRing*>::iterator i = s_allGeoRings.begin();
				i != s_allGeoRings.end(); ++i) {
			if ((*i)->m_runUpdateThread) (*i)->_UpdateLODs();
		}
		SDL_mutexV(s_allGeoRingsLock);

		SDL_Delay(10);
	}
	return 0;
}

void GeoRing::_UpdateLODs()
{
	//PROFILE_SCOPED()
	for (size_t i=0; i<m_plates.size(); i++) {
		m_plates[i]->LODUpdate(m_tempCampos);
	}
	m_runUpdateThread = 0;
}

/* This is to stop threads keeping on iterating over the s_allGeoRings list,
 * which may have been destroyed by exit() (does on lunix anyway...)
 */
static void _LockoutThreadsBeforeExit()
{
	SDL_mutexP(s_allGeoRingsLock);
	Profiler::dumphtml();
}

void GeoRing::Init()
{
	//PROFILE_SCOPED()
	s_geoRingSurfaceShader[0] = new GeoRingShader("geoRing", "#define NUM_LIGHTS 1\n");
	s_geoRingSurfaceShader[1] = new GeoRingShader("geoRing", "#define NUM_LIGHTS 2\n");
	s_geoRingSurfaceShader[2] = new GeoRingShader("geoRing", "#define NUM_LIGHTS 3\n");
	s_geoRingSurfaceShader[3] = new GeoRingShader("geoRing", "#define NUM_LIGHTS 4\n");
	/*s_geoRingSkyShader[0] = new GeoRingShader("geosphere_sky", "#define NUM_LIGHTS 1\n");
	s_geoRingSkyShader[1] = new GeoRingShader("geosphere_sky", "#define NUM_LIGHTS 2\n");
	s_geoRingSkyShader[2] = new GeoRingShader("geosphere_sky", "#define NUM_LIGHTS 3\n");
	s_geoRingSkyShader[3] = new GeoRingShader("geosphere_sky", "#define NUM_LIGHTS 4\n");*/
	s_allGeoRingsLock = SDL_CreateMutex();
	OnChangeDetailLevel();
#ifdef GEORING_USE_THREADING
	SDL_CreateThread(&GeoRing::UpdateLODThread, 0);
#endif /* GEORING_USE_THREADING */
	atexit(&_LockoutThreadsBeforeExit);
}

void GeoRing::OnChangeDetailLevel()
{
	//PROFILE_SCOPED()
	SDL_mutexP(s_allGeoRingsLock);
	for(std::list<GeoRing*>::iterator i = s_allGeoRings.begin(); i != s_allGeoRings.end(); ++i) {
		// remove the plates
		for (size_t p=0; p<(*i)->m_plates.size(); p++) {
			if ((*i)->m_plates[p]) {
				delete (*i)->m_plates[p];
			}
		}
		(*i)->m_plates.clear();
	}

	switch (Pi::detail.planets) {
		case 0: GEOPLATE_EDGELEN = 15; break;
		case 1: GEOPLATE_EDGELEN = 25; break;
		case 2: GEOPLATE_EDGELEN = 45; break;
		case 3: GEOPLATE_EDGELEN = 65; break;
		default:
		case 4: GEOPLATE_EDGELEN = GEOPLATE_MAX_EDGELEN; break;
	}
	assert(GEOPLATE_EDGELEN <= GEOPLATE_MAX_EDGELEN);
	GeoPlate::Init();
	for(std::list<GeoRing*>::iterator i = s_allGeoRings.begin(); i != s_allGeoRings.end(); ++i) {
		(*i)->BuildFirstPatches();
	}
	SDL_mutexV(s_allGeoRingsLock);
}

#define GEOSPHERE_TYPE	(m_sbody->type)

GeoRing::GeoRing(const SBody *body): m_style(body)
{
	//PROFILE_SCOPED()
	m_vbosToDestroyLock = SDL_CreateMutex();
	m_runUpdateThread = 0;
	m_sbody = body;

	SDL_mutexP(s_allGeoRingsLock);
	s_allGeoRings.push_back(this);
	SDL_mutexV(s_allGeoRingsLock);

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
	//PROFILE_SCOPED()
	SDL_mutexP(s_allGeoRingsLock);
	s_allGeoRings.remove(this);
	SDL_mutexV(s_allGeoRingsLock);

	for (size_t i=0; i<m_plates.size(); i++) {
		if (m_plates[i]) {
			delete m_plates[i];
		}
	}
	DestroyVBOs();
	SDL_DestroyMutex(m_vbosToDestroyLock);
}

void GeoRing::AddVBOToDestroy(GLuint vbo)
{
	//PROFILE_SCOPED()
	SDL_mutexP(m_vbosToDestroyLock);
	m_vbosToDestroy.push_back(vbo);
	SDL_mutexV(m_vbosToDestroyLock);
}

void GeoRing::DestroyVBOs()
{
	//PROFILE_SCOPED()
	SDL_mutexP(m_vbosToDestroyLock);
	for (std::list<GLuint>::iterator i = m_vbosToDestroy.begin();
			i != m_vbosToDestroy.end(); ++i) {
		glDeleteBuffersARB(1, &(*i));
	}
	m_vbosToDestroy.clear();
	SDL_mutexV(m_vbosToDestroyLock);
}

GeoPlate* GeoRing::FindGeoPlateByIndex(const int idx) const
{
	const int nPlates = m_plates.size();
	if( idx<0 ) {
		const int negidx = (nPlates + idx);
		return m_plates[negidx];
	} else if( idx>=nPlates ) {
		const int wrapidx = idx % nPlates;
		return m_plates[wrapidx];
	}
	// default
	return m_plates[idx];
}

void GeoRing::BuildFirstPatches(const int numSegments)
{
	//PROFILE_SCOPED()
	std::vector<vector3d>	points;
	std::vector<double>		angles;
    for (int i = 0; i <= numSegments; ++i) {
		double angle = ((M_PI * 360.0) / 180.0) * i / numSegments;
		double sinval = sin(angle);
		double cosval = cos(angle);
		vector3d vp(sinval, cosval, 0.0 );
		points.push_back( vp.Normalized() );
		angles.push_back( angle );
    }

	// calculate the width of the orbital
	mRingWidth = (points[1] - points[0]).Length();

	// build the terrain plates
	m_plates.clear();
	for( int i=0 ; i<points.size()-1 ; ++i ) {
		m_plates.push_back( new GeoPlate(mRingWidth, angles[i], angles[i+1], 0.0, 0, 0) );
	}

	// set edge friends
	for (size_t i=0; i<m_plates.size(); i++) {
		m_plates[i]->geoRing = this;

		// trailing edge?
		m_plates[i]->edgeFriend[0] = 0;
		m_plates[i]->edgeFriend[1] = FindGeoPlateByIndex(i-1);	// backward?
		m_plates[i]->edgeFriend[2] = 0;
		m_plates[i]->edgeFriend[3] = FindGeoPlateByIndex(i+1);	// forward?
	}

	// create mesh(es)
	for (size_t i=0; i<m_plates.size(); i++) m_plates[i]->GenerateMesh();
	for (size_t i=0; i<m_plates.size(); i++) m_plates[i]->GenerateEdgeNormalsAndColors();
	for (size_t i=0; i<m_plates.size(); i++) m_plates[i]->UpdateVBOs();

	// hacking
	/*vector3d colorIdx[4] = { 
		vector3d( 1.0, 0.0, 0.0 ),	// red
		vector3d( 0.0, 1.0, 0.0 ),	// green
		vector3d( 0.0, 0.0, 1.0 ),	// blue
		vector3d( 1.0, 0.0, 1.0 )	// purple
	};
	for (size_t i=0; i<m_plates.size(); i++) {
		for (int y=0; y<GEOPLATE_EDGELEN; y++) {
			for (int x=0; x<GEOPLATE_EDGELEN; x++) {
				m_plates[i]->colors[x + y*GEOPLATE_EDGELEN] = colorIdx[ m_plates[i]->m_cIdx ];
			}
		}
	}*/
}

static const float g_ambient[4] = { 0, 0, 0, 1.0 };

static void DrawAtmosphereSurface(const vector3d &campos, float rad)
{
	//PROFILE_SCOPED()
	const int LAT_SEGS = 20;
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

	/* Tri-fan above viewer */
	glBegin(GL_TRIANGLE_FAN);
	glVertex3f(0.0f, 1.0f, 0.0f);
	for (int i=0; i<=LONG_SEGS; i++) {
		glVertex3f(sin(latDiff)*sinCosTable[i][0], cos(latDiff), -sin(latDiff)*sinCosTable[i][1]);
	}
	glEnd();

	/* and wound latitudinal strips */
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

	glPopMatrix();
}

void GeoRing::Render(vector3d campos, const float radius, const float scale) {
	//PROFILE_SCOPED()
	Plane planes[6];
	GetFrustum(planes);
	const float atmosRadius = ATMOSPHERE_RADIUS;
	
	// no frustum test of entire geoRing, since Space::Render does this
	// for each body using its GetBoundingRadius() value
	if (Render::AreShadersEnabled()) {
		Color atmosCol;
		double atmosDensity;
		matrix4x4d modelMatrix;
		glGetDoublev (GL_MODELVIEW_MATRIX, &modelMatrix[0]);
		vector3d center = modelMatrix * vector3d(0.0, 0.0, 0.0);
		
		GetAtmosphereFlavor(&atmosCol, &atmosDensity);
		atmosDensity *= 0.00005;

		/*if (atmosDensity != 0.0) {
			GeoRingShader *shader = s_geoRingSkyShader[Render::State::GetNumLights()-1];
			Render::State::UseProgram(shader);
			shader->set_geosphereScale(scale);
			shader->set_geosphereAtmosTopRad(atmosRadius*radius/scale);
			shader->set_geosphereAtmosFogDensity(atmosDensity);
			shader->set_atmosColor(atmosCol.r, atmosCol.g, atmosCol.b, atmosCol.a);
			shader->set_geosphereCenter(center.x, center.y, center.z);
			
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE);
			// make atmosphere sphere slightly bigger than required so
			// that the edges of the pixel shader atmosphere jizz doesn't
			// show ugly polygonal angles
			DrawAtmosphereSurface(campos, atmosRadius*1.01);
			glDisable(GL_BLEND);
		}*/

		GeoRingShader *shader = s_geoRingSurfaceShader[Render::State::GetNumLights()-1];
		Render::State::UseProgram(shader);
		shader->set_geosphereScale(scale);
		shader->set_geosphereAtmosTopRad(atmosRadius*radius/scale);
		shader->set_geosphereAtmosFogDensity(atmosDensity);
		shader->set_atmosColor(atmosCol.r, atmosCol.g, atmosCol.b, atmosCol.a);
		shader->set_geosphereCenter(center.x, center.y, center.z);
	}

	if (0==m_plates.size()) {
		BuildFirstPatches();
	}

	const float black[4] = { 0,0,0,0 };
	float ambient[4];// = { 0.1, 0.1, 0.1, 1.0 };

	// save old global ambient
	float oldAmbient[4];
	glGetFloatv(GL_LIGHT_MODEL_AMBIENT, oldAmbient);

	// give planet some ambient lighting if the viewer is close to it
	/*{
		double camdist = campos.Length();
		camdist = 0.1 / (camdist*camdist);
		// why the fuck is this returning 0.1 when we are sat on the planet??
		// XXX oh well, it is the value we want anyway...
		ambient[0] = ambient[1] = ambient[2] = float(camdist);
		ambient[3] = 1.0f;
	}*/
	
	glLightModelfv (GL_LIGHT_MODEL_AMBIENT, ambient);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glMaterialfv (GL_FRONT, GL_SPECULAR, black);
	glMaterialfv (GL_FRONT, GL_EMISSION, black);
	glEnable(GL_COLOR_MATERIAL);

	glPolygonMode(GL_FRONT, GL_FILL);
	for (size_t i=0; i<m_plates.size(); i++) {
		m_plates[i]->Render(campos, planes);
	}
	Render::State::UseProgram(0);

	glDisable(GL_COLOR_MATERIAL);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, oldAmbient);

	// if the update thread has deleted any geopatches, destroy the vbos
	// associated with them
	DestroyVBOs();
	
	if (!m_runUpdateThread) {
		this->m_tempCampos = campos;
		m_runUpdateThread = 1;
	}
#ifndef GEORING_USE_THREADING
	m_tempCampos = campos;
	_UpdateLODs();
#endif /* !GEORING_USE_THREADING */
}

double GeoRing::GetDistFromSurface(const vector3d p) {
	//return DBL_MAX;
	// one way to find angle (radians)
	double radang = atan2(p.y,p.x);
	//double radang = atan2(p.x, p.y);
	if( radang<0.0 ) {
		radang = (2*M_PI) + radang;
	}

	// iterate through plates until we find the one our angles/z are within
	PlateIter iter = m_plates.begin();
	for( ; iter!=m_plates.end() ; ++iter ) {
		// scale the p.z into the 0.0 to 1.0 range
		const double realZ = (p.z + (*iter)->m_halfLen) / (2.0*(*iter)->m_halfLen);
		if( (*iter)->IsWithinAng(radang,p.z) ) {
			// then get distance from surface
			return (*iter)->DistFromSurface(radang, realZ);
		}
	}

	/*const vector3d CP1(0.0, 0.0, mRingWidth);		// vertex at middle of top end
	const vector3d CP2(0.0, 0.0, -mRingWidth);	// vertex at middle of bottom end
	const vector3d CN1 = (CP2 - CP1).Normalized();
	const vector3d CN2 = -CN1;
	const double fDistanceToPlane1 = (p-CP1).Dot( CN1 );
	const double fDistanceToPlane2 = (p-CP2).Dot( CN2 );
	if (fDistanceToPlane1 < 0.0)	
		return DBL_MAX;	
	if (fDistanceToPlane2 < 0.0) 
		return DBL_MAX; 
	const vector3d TempP =  p - (CN1 * fDistanceToPlane1);
	const double fDistanceFromCenter = (TempP-CP1).Length();
	if (fDistanceFromCenter > 1.0) 
		return DBL_MAX;	*/

	return DBL_MAX; // All tests passed, point is in cylinder
}
