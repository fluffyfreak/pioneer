#include "libs.h"
#include "GeoRing.h"
#include "perlin.h"
#include "Pi.h"
#include "StarSystem.h"
#include "render/Render.h"

#include "profiler/Profiler.h"

#include "vcacheopt.h"

#define GEOPLATE_USE_THREADING

// tri edge lengths
#define GEOPLATE_SUBDIVIDE_AT_CAMDIST	1.5 //5.0
#define GEOPLATE_MAX_DEPTH	15
// must be an odd number
#define GEOPLATE_EDGELEN_DEFAULT	7//15
#define GEOPLATE_NUMVERTICES	(GEOPLATE_EDGELEN*GEOPLATE_EDGELEN)
#define GEORING_USE_THREADING

int GEOPLATE_EDGELEN = GEOPLATE_EDGELEN_DEFAULT;	// 15;
static const int GEOPLATE_MAX_EDGELEN = 55;
static double GEOPLATE_FRAC;
static double GEOPLATEHULL_FRAC;
static double GEOPLATEWALL_FRAC;

#define lerp(t, a, b) ( a + t * (b - a) )

#define PRINT_VECTOR(_v) printf("%f,%f,%f\n", (_v).x, (_v).y, (_v).z);

#define SAFE_DELETE(d) delete d; d=0;
#define SAFE_DELETE_ARRAY(d) delete[] d; d=0;

SHADER_CLASS_BEGIN(GeoRingShader)
	SHADER_UNIFORM_VEC4(atmosColor)
	SHADER_UNIFORM_FLOAT(georingScale)
	SHADER_UNIFORM_FLOAT(georingAtmosTopRad)
	SHADER_UNIFORM_VEC3(georingCenter)
	SHADER_UNIFORM_FLOAT(georingAtmosFogDensity)
SHADER_CLASS_END()

static GeoRingShader *s_geoRingSurfaceShader[4], *s_geoRingSkyShader[4];

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

// hold the 16 possible terrain edge connections
const int NUM_INDEX_LISTS = 16;
typedef struct {
	std::vector<unsigned short> v;
} VecShortHolder;

#define NUM_VBE 2

class BaseGeo {
public:
	BaseGeo(const double halfLength, const double startAng, const double endAng, 
		const double yoffset, const int depth) : vertices(NULL), normals(NULL), colors(NULL), geoRing(NULL), m_vbo(0)
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
		
		const vector3d topEndEdge(sin(theta), m_yoffset + halfLength, cos(theta));		// vertices at top edge of circle
		const vector3d bottomEndEdge(sin(theta), m_yoffset - halfLength, cos(theta));	// vertices at bottom edge of circle
		
		const vector3d res = lerp( y, bottomEndEdge, topEndEdge );
		return res;
	}

	void SetGeoRingPtr( GeoRing *ptr ) { PiAssert(ptr!=NULL); geoRing = ptr; }

protected:
	double ang[NUM_VBE];
	double m_halfLen;
	double m_yoffset;		// offset from the centre i.e. newyoffset = (m_yoffset + m_halfLen);
	vector3d *vertices;
	vector3d *normals;
	vector3d *colors;
	GeoRing *geoRing;
	GLuint m_vbo;
	double m_roughLength;
	vector3d clipCentroid;
	double clipRadius;
};

class GeoPlateHull : public BaseGeo {
public:
	static unsigned short *indices;
	static GLuint indices_vbo;
	static VBOVertex *vbotemp;

	// params
	// v0, v1 - define points on the line describing the loop of the ring/orbital.
	// depth - 0 is the topmost plate with each depth+1 describing it's depth within the tree.
	GeoPlateHull(const double halfLength, const double startAng, const double endAng, 
		const double yoffset, const int depth) : BaseGeo(halfLength, startAng, endAng, yoffset, depth) 
	{
		normals = new vector3d[GEOPLATE_NUMVERTICES];
		vertices = new vector3d[GEOPLATE_NUMVERTICES];
		colors = new vector3d[GEOPLATE_NUMVERTICES];
	}

	~GeoPlateHull() {
		delete vertices;
		delete normals;
		delete colors;
		geoRing->AddVBOToDestroy(m_vbo);
	}

	static void Init() {
		GEOPLATEHULL_FRAC = 1.0 / double(GEOPLATE_EDGELEN-1);

		if (indices) {
			delete [] indices;
			if (indices_vbo) {
				glDeleteBuffersARB(1, &indices_vbo);
			}
			delete [] vbotemp;
		}
		{
			vbotemp = new VBOVertex[GEOPLATE_NUMVERTICES];
			unsigned short *idx;
			indices = new unsigned short[VBO_COUNT_ALL_IDX];
			idx = indices;
			for (int x=0; x<GEOPLATE_EDGELEN-1; x++) {
				for (int y=0; y<GEOPLATE_EDGELEN-1; y++) {
					idx[0] = x + GEOPLATE_EDGELEN*y;		PiAssert(idx[0] < GEOPLATE_NUMVERTICES);
					idx[1] = x+1 + GEOPLATE_EDGELEN*y;		PiAssert(idx[1] < GEOPLATE_NUMVERTICES);
					idx[2] = x + GEOPLATE_EDGELEN*(y+1);	PiAssert(idx[2] < GEOPLATE_NUMVERTICES);
					idx+=3;

					idx[0] = x+1 + GEOPLATE_EDGELEN*y;		PiAssert(idx[0] < GEOPLATE_NUMVERTICES);
					idx[1] = x+1 + GEOPLATE_EDGELEN*(y+1);	PiAssert(idx[1] < GEOPLATE_NUMVERTICES);
					idx[2] = x + GEOPLATE_EDGELEN*(y+1);	PiAssert(idx[2] < GEOPLATE_NUMVERTICES);
					PiAssert(idx < indices+VBO_COUNT_ALL_IDX);
					idx+=3;
				}
			}

			glGenBuffersARB(1, &indices_vbo);
			glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, indices_vbo);
			glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned short)*VBO_COUNT_ALL_IDX, 0, GL_STATIC_DRAW);
			glBufferSubDataARB(GL_ELEMENT_ARRAY_BUFFER, 0, sizeof(unsigned short)*VBO_COUNT_ALL_IDX, indices);
			glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);
		}
	}

	void UpdateVBOs() {
		if (!m_vbo) glGenBuffersARB(1, &m_vbo);
		glBindBufferARB(GL_ARRAY_BUFFER, m_vbo);
		glBufferDataARB(GL_ARRAY_BUFFER, sizeof(VBOVertex)*GEOPLATE_NUMVERTICES, 0, GL_DYNAMIC_DRAW);
		for (int i=0; i<GEOPLATE_NUMVERTICES; i++)
		{
			clipRadius = std::max(clipRadius, (vertices[i]-clipCentroid).Length());
			VBOVertex *pData = vbotemp + i;
			pData->x = float(vertices[i].x - clipCentroid.x);
			pData->y = float(vertices[i].y - clipCentroid.y);
			pData->z = float(vertices[i].z - clipCentroid.z);
			pData->nx = float(normals[i].x);
			pData->ny = float(normals[i].y);
			pData->nz = float(normals[i].z);
			pData->col[0] = static_cast<unsigned char>(Clamp(colors[i].x*255.0, 0.0, 255.0));
			pData->col[1] = static_cast<unsigned char>(Clamp(colors[i].y*255.0, 0.0, 255.0));
			pData->col[2] = static_cast<unsigned char>(Clamp(colors[i].z*255.0, 0.0, 255.0));
			pData->col[3] = 255;
		}
		glBufferDataARB(GL_ARRAY_BUFFER, sizeof(VBOVertex)*GEOPLATE_NUMVERTICES, vbotemp, GL_DYNAMIC_DRAW);
		glBindBufferARB(GL_ARRAY_BUFFER, 0);
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
				PiAssert(x!=GEOPLATE_EDGELEN);
				PiAssert(y!=GEOPLATE_EDGELEN);
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
				double col = 0.5;
				colors[x + y*GEOPLATE_EDGELEN] = vector3d(col, col, 1.0);
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

		// Generate bottom row normals
		for (int y=0; y<GEOPLATE_EDGELEN-1; ++y) {
			const int x=GEOPLATE_EDGELEN-1;
			// normal
			vector3d xy = vertices[x + y*GEOPLATE_EDGELEN];
			vector3d x1 = vertices[x-1 + y*GEOPLATE_EDGELEN];
			vector3d y1 = vertices[x + (y+1)*GEOPLATE_EDGELEN];

			vector3d n = (xy-x1).Cross(y1-xy);
			normals[x + y*GEOPLATE_EDGELEN] = (n.Normalized());
		}

		// Generate last col normals
		{
			const int y=GEOPLATE_EDGELEN-1;
			for (int x=1; x<GEOPLATE_EDGELEN; ++x) {
				// normal
				vector3d xy = vertices[x + y*GEOPLATE_EDGELEN];
				vector3d x1 = vertices[x-1 + y*GEOPLATE_EDGELEN];
				vector3d y1 = vertices[x + (y-1)*GEOPLATE_EDGELEN];

				vector3d n = (xy-x1).Cross(xy-y1);
				normals[x + y*GEOPLATE_EDGELEN] = (n.Normalized());
			}
		}

		// corner normals
		{
			const int y=GEOPLATE_EDGELEN-1;
			const int x=0;
			{
				// normal
				vector3d xy = vertices[x + y*GEOPLATE_EDGELEN];
				vector3d x1 = vertices[x+1 + y*GEOPLATE_EDGELEN];
				vector3d y1 = vertices[x + (y-1)*GEOPLATE_EDGELEN];

				vector3d n = (x1-xy).Cross(xy-y1);
				normals[x + y*GEOPLATE_EDGELEN] = (n.Normalized());
			}
		}
		{
			const int y=0;
			const int x=GEOPLATE_EDGELEN-1;
			{
				// normal
				vector3d xy = vertices[x + y*GEOPLATE_EDGELEN];
				vector3d x1 = vertices[x-1 + y*GEOPLATE_EDGELEN];
				vector3d y1 = vertices[x + (y+1)*GEOPLATE_EDGELEN];

				vector3d n = (xy-x1).Cross(y1-xy);
				normals[x + y*GEOPLATE_EDGELEN] = (n.Normalized());
			}
		}
	}

	void Render(vector3d &campos, Plane planes[6]) {
		// frustum test!
		for (int i=0; i<6; i++) {
			if (planes[i].DistanceToPoint(clipCentroid) <= -clipRadius) {
				return;
			}
		}

		vector3d relpos = clipCentroid - campos;
		glPushMatrix();
		glTranslated(relpos.x, relpos.y, relpos.z);

		Pi::statSceneTris += 2*(GEOPLATE_EDGELEN-1)*(GEOPLATE_EDGELEN-1);
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);
		glEnableClientState(GL_COLOR_ARRAY);

		glBindBufferARB(GL_ARRAY_BUFFER, m_vbo);
		glVertexPointer(3, GL_FLOAT, sizeof(VBOVertex), 0);
		glNormalPointer(GL_FLOAT, sizeof(VBOVertex), reinterpret_cast<void *>(3*sizeof(float)));
		glColorPointer(4, GL_UNSIGNED_BYTE, sizeof(VBOVertex), reinterpret_cast<void *>(6*sizeof(float)));
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, indices_vbo);
		glDrawElements(GL_TRIANGLES, VBO_COUNT_ALL_IDX, GL_UNSIGNED_SHORT, 0);
		glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_COLOR_ARRAY);
		glPopMatrix();
	}

	void RenderNormals(vector3d &campos, Plane planes[6]) {
		for (int i=0; i<6; i++) {
			if (planes[i].DistanceToPoint(clipCentroid) <= -clipRadius) {
				return;
			}
		}
		glLineWidth(2.0);
		glBegin(GL_LINES);
		glColor3f( 1.0f, 1.0f, 1.0f );
		for (int y=0; y<GEOPLATE_EDGELEN; ++y) {	// across the width
			for (int x=0; x<GEOPLATE_EDGELEN; ++x) {	// along the length (circumference)
				vector3d xy = vertices[x + y*GEOPLATE_EDGELEN];
				vector3d nxy = (normals[x + y*GEOPLATE_EDGELEN])*0.01;

				// normal
				glVertex3dv( &(xy)[0] );
				glVertex3dv( &(xy+nxy)[0] );
			}
		}
		glEnd();
	}
};
//static 
unsigned short *	GeoPlateHull::indices = 0;
GLuint				GeoPlateHull::indices_vbo = 0;
VBOVertex *			GeoPlateHull::vbotemp= 0 ;

// must be an odd number
#define GEOPLATE_WALL_LEN 7
#define GEOPLATE_NUM_WALL_VERTICES	(GEOPLATE_WALL_LEN*GEOPLATE_WALL_LEN)

class GeoPlateWall : public BaseGeo {
public:
	static unsigned short *indices;
	static GLuint indices_vbo;
	static VBOVertex *vbotemp;
	bool bInnerWall;

	// params
	// v0, v1 - define points on the line describing the loop of the ring/orbital.
	// depth - 0 is the topmost plate with each depth+1 describing it's depth within the tree.
	GeoPlateWall(const bool isInnerWall, const double halfLength, const double startAng, const double endAng, 
		const double yoffset, const int depth) : BaseGeo(halfLength, startAng, endAng, yoffset, depth), bInnerWall(isInnerWall)
	{
		normals = new vector3d[GEOPLATE_NUM_WALL_VERTICES];
		vertices = new vector3d[GEOPLATE_NUM_WALL_VERTICES];
		colors = new vector3d[GEOPLATE_NUM_WALL_VERTICES];
	}

	~GeoPlateWall() {
		delete vertices;
		delete normals;
		delete colors;
		geoRing->AddVBOToDestroy(m_vbo);
	}

	static void Init() {
		GEOPLATEWALL_FRAC = 1.0 / double(GEOPLATE_WALL_LEN-1);

		if (indices) {
			delete [] indices;
			if (indices_vbo) {
				glDeleteBuffersARB(1, &indices_vbo);
			}
			delete [] vbotemp;
		}
		{
			vbotemp = new VBOVertex[GEOPLATE_NUM_WALL_VERTICES];
			unsigned short *idx;
			indices = new unsigned short[VBO_COUNT_ALL_IDX];
			idx = indices;
			for (int x=0; x<GEOPLATE_WALL_LEN-1; x++) {
				for (int y=0; y<GEOPLATE_WALL_LEN-1; y++) {
					idx[0] = x + GEOPLATE_WALL_LEN*y;		PiAssert(idx[0] < GEOPLATE_NUM_WALL_VERTICES);
					idx[1] = x+1 + GEOPLATE_WALL_LEN*y;		PiAssert(idx[1] < GEOPLATE_NUM_WALL_VERTICES);
					idx[2] = x + GEOPLATE_WALL_LEN*(y+1);	PiAssert(idx[2] < GEOPLATE_NUM_WALL_VERTICES);
					idx+=3;

					idx[0] = x+1 + GEOPLATE_WALL_LEN*y;		PiAssert(idx[0] < GEOPLATE_NUM_WALL_VERTICES);
					idx[1] = x+1 + GEOPLATE_WALL_LEN*(y+1);	PiAssert(idx[1] < GEOPLATE_NUM_WALL_VERTICES);
					idx[2] = x + GEOPLATE_WALL_LEN*(y+1);	PiAssert(idx[2] < GEOPLATE_NUM_WALL_VERTICES);
					PiAssert(idx < indices+VBO_COUNT_ALL_IDX);
					idx+=3;
				}
			}

			glGenBuffersARB(1, &indices_vbo);
			glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, indices_vbo);
			glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned short)*VBO_COUNT_ALL_IDX, 0, GL_STATIC_DRAW);
			glBufferSubDataARB(GL_ELEMENT_ARRAY_BUFFER, 0, sizeof(unsigned short)*VBO_COUNT_ALL_IDX, indices);
			glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);
		}
	}

	void UpdateVBOs() {
		if (!m_vbo) glGenBuffersARB(1, &m_vbo);
		glBindBufferARB(GL_ARRAY_BUFFER, m_vbo);
		glBufferDataARB(GL_ARRAY_BUFFER, sizeof(VBOVertex)*GEOPLATE_NUM_WALL_VERTICES, 0, GL_DYNAMIC_DRAW);
		for (int i=0; i<GEOPLATE_NUM_WALL_VERTICES; i++)
		{
			clipRadius = std::max(clipRadius, (vertices[i]-clipCentroid).Length());
			VBOVertex *pData = vbotemp + i;
			pData->x = float(vertices[i].x - clipCentroid.x);
			pData->y = float(vertices[i].y - clipCentroid.y);
			pData->z = float(vertices[i].z - clipCentroid.z);
			pData->nx = float(normals[i].x);
			pData->ny = float(normals[i].y);
			pData->nz = float(normals[i].z);
			pData->col[0] = static_cast<unsigned char>(Clamp(colors[i].x*255.0, 0.0, 255.0));
			pData->col[1] = static_cast<unsigned char>(Clamp(colors[i].y*255.0, 0.0, 255.0));
			pData->col[2] = static_cast<unsigned char>(Clamp(colors[i].z*255.0, 0.0, 255.0));
			pData->col[3] = 255;
		}
		glBufferDataARB(GL_ARRAY_BUFFER, sizeof(VBOVertex)*GEOPLATE_NUM_WALL_VERTICES, vbotemp, GL_DYNAMIC_DRAW);
		glBindBufferARB(GL_ARRAY_BUFFER, 0);
	}

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

				double col = 0.5;
				colors[x + y*GEOPLATE_WALL_LEN] = vector3d(col, col, 1.0);
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

	void Render(vector3d &campos, Plane planes[6]) {
		// frustum test!
		for (int i=0; i<6; i++) {
			if (planes[i].DistanceToPoint(clipCentroid) <= -clipRadius) {
				return;
			}
		}

		vector3d relpos = clipCentroid - campos;
		glPushMatrix();
		glTranslated(relpos.x, relpos.y, relpos.z);

		Pi::statSceneTris += 2*(GEOPLATE_WALL_LEN-1)*(GEOPLATE_WALL_LEN-1);
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);
		glEnableClientState(GL_COLOR_ARRAY);

		glBindBufferARB(GL_ARRAY_BUFFER, m_vbo);
		glVertexPointer(3, GL_FLOAT, sizeof(VBOVertex), 0);
		glNormalPointer(GL_FLOAT, sizeof(VBOVertex), reinterpret_cast<void *>(3*sizeof(float)));
		glColorPointer(4, GL_UNSIGNED_BYTE, sizeof(VBOVertex), reinterpret_cast<void *>(6*sizeof(float)));
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, indices_vbo);
		glDrawElements(GL_TRIANGLES, VBO_COUNT_ALL_IDX, GL_UNSIGNED_SHORT, 0);
		glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_COLOR_ARRAY);
		glPopMatrix();
	}

	void RenderNormals(vector3d &campos, Plane planes[6]) {
		for (int i=0; i<6; i++) {
			if (planes[i].DistanceToPoint(clipCentroid) <= -clipRadius) {
				return;
			}
		}
		glLineWidth(2.0);
		glBegin(GL_LINES);
		glColor3f( 1.0f, 1.0f, 1.0f );
		for (int y=0; y<GEOPLATE_WALL_LEN; ++y) {	// across the width
			for (int x=0; x<GEOPLATE_WALL_LEN; ++x) {	// along the length (circumference)
				vector3d xy = vertices[x + y*GEOPLATE_WALL_LEN];
				vector3d nxy = (normals[x + y*GEOPLATE_WALL_LEN])*0.01;

				// normal
				glVertex3dv( &(xy)[0] );
				glVertex3dv( &(xy+nxy)[0] );
			}
		}
		glEnd();
	}
};
//static 
unsigned short *	GeoPlateWall::indices = 0;
GLuint				GeoPlateWall::indices_vbo = 0;
VBOVertex *			GeoPlateWall::vbotemp= 0 ;

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

#define NUM_PLATE_INDICES 16
class GeoPlate : public BaseGeo {
public:
	vector3d vbe[NUM_VBE];
	GeoPlate *kids[4];
	GeoPlate *parent;
	GeoPlate *edgeFriend[4]; // [0]=v01, [1]=v12, [2]=v23, [3]=30 - original from GeoPatch...
	int m_depth;
	int m_cIdx;
	SDL_mutex *m_kidsLock;
	bool m_needUpdateVBOs;

	static unsigned short *midIndices;
	static unsigned short *loEdgeIndices[4];
	static unsigned short *hiEdgeIndices[4];
	static GLuint indices_vbo;
	static GLuint indices_list[NUM_PLATE_INDICES];
	static GLuint indices_tri_count;
	static GLuint indices_tri_counts[NUM_PLATE_INDICES];
	static VBOVertex *vbotemp;

	static GeoPlate** s_geoPlates[4];
	static SDL_mutex *s_geoPlateLock[4];
	static SDL_sem* s_geoPlateSem[4];
	static SDL_sem* s_geoRingSem[4];
	static SDL_Thread* s_geoRingThread[4];

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
	GeoPlate(const double halfLength, const vector3d &startVBE, const vector3d &endVBE, 
		const double startAng, const double endAng, 
		const double yoffset, const int depth, const int cIdx) : BaseGeo(halfLength, startAng, endAng, yoffset, depth) 
		, parent(0), m_kidsLock(0), m_needUpdateVBOs(false)
	{
		for( int i=0; i<4; ++i )
		{
			kids[i] = 0;
			edgeFriend[i] = 0;
		}

		m_kidsLock = SDL_CreateMutex();
		vbe[0] = startVBE;
		vbe[1] = endVBE;
		m_depth = depth;
		m_cIdx = cIdx;
		m_needUpdateVBOs = false;
		normals = new vector3d[GEOPLATE_NUMVERTICES];
		vertices = new vector3d[GEOPLATE_NUMVERTICES];
		colors = new vector3d[GEOPLATE_NUMVERTICES];
	}

	~GeoPlate() {
		SDL_DestroyMutex(m_kidsLock);
		for (int i=0; i<4; i++) {
			if (edgeFriend[i]) edgeFriend[i]->NotifyEdgeFriendDeleted(this);
		}
		for (int i=0; i<4; i++) if (kids[i]) delete kids[i];
		SAFE_DELETE_ARRAY(vertices);
		SAFE_DELETE_ARRAY(normals);
		SAFE_DELETE_ARRAY(colors);
		geoRing->AddVBOToDestroy(m_vbo);
	}

	void updateIndexBufferId() {
		GLuint acc = 0;
		for (int i=0; i<4; i++) {
			GLuint x = (edgeFriend[i]) ? 1 : 0;
			acc = acc | (x<<i);
		}
		PiAssert(acc<NUM_INDEX_LISTS);
		indices_vbo = indices_list[acc];
		indices_tri_count = indices_tri_counts[acc];
	}

	template<typename T>
	static int getIndices(std::vector<T> &pl, const bool *isHigh)
	{
		// calculate how many tri's there are
		int tri_count = (VBO_COUNT_MID_IDX / 3); 
		for( int i=0; i<4; ++i ) {
			if( isHigh[i] ) {
				tri_count += (VBO_COUNT_HI_EDGE / 3);
			} else {
				tri_count += (VBO_COUNT_LO_EDGE / 3);
			}
		}

		// pre-allocate enough space
		pl.reserve(tri_count);

		// add all of the middle indices
		for(int i=0; i<VBO_COUNT_MID_IDX; ++i) {
			pl.push_back(midIndices[i]);
		}
		// selectively add the HI or LO detail indices
		for (int i=0; i<4; i++) {
			if( isHigh[i] ) {
				for(int j=0; j<VBO_COUNT_HI_EDGE; ++j) {
					pl.push_back(hiEdgeIndices[i][j]);
				}
			} else {
				for(int j=0; j<VBO_COUNT_LO_EDGE; ++j) {
					pl.push_back(loEdgeIndices[i][j]);
				}
			}
		}

		return tri_count;
	}

	static void Init() {
		GEOPLATE_FRAC = 1.0 / double(GEOPLATE_EDGELEN-1);

		if (midIndices) {
			SAFE_DELETE_ARRAY(midIndices);
			for (int i=0; i<4; i++) {
				SAFE_DELETE_ARRAY(loEdgeIndices[i]);
				SAFE_DELETE_ARRAY(hiEdgeIndices[i]);
			}
			SAFE_DELETE_ARRAY(vbotemp);
		}

		if (indices_vbo) {
			glDeleteBuffersARB(1, &indices_vbo);
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

			// these will hold the optimised indices
			VecShortHolder pl_short[NUM_INDEX_LISTS];
			{
				// an initial bunch of values that will get used to sort out the above and following real lists
				typedef struct {
					bool e[4];
				} edgebools;
				edgebools e[NUM_INDEX_LISTS] = {
					{0000,0000,0000,0000},
					{true,0000,0000,0000},
					{0000,true,0000,0000},
					{true,true,0000,0000},
					{0000,0000,true,0000},
					{true,0000,true,0000},
					{0000,true,true,0000},
					{true,true,true,0000},
					{0000,0000,0000,true},
					{true,0000,0000,true},
					{0000,true,0000,true},
					{true,true,0000,true},
					{0000,0000,true,true},
					{true,0000,true,true},
					{0000,true,true,true},
					{true,true,true,true}
				};

				// Calculate the index based on the four binary ops
				// then generate the entry at that index. 
				// Means I don't have to work it out by hand ;)
				unsigned char acc = 0;
				edgebools test_e[NUM_INDEX_LISTS] = {0};
				for( int i=0; i<NUM_INDEX_LISTS; ++i ) {
					acc=0;
					for( size_t j=0; j<4; ++j ) {
						size_t x = e[i].e[j] ? 1 : 0;
						acc = acc | (x<<j);
					}
					assert(acc<16);
					test_e[acc].e[0] = (acc & 0x01) ? true : false;
					test_e[acc].e[1] = (acc & 0x02) ? true : false;
					test_e[acc].e[2] = (acc & 0x04) ? true : false;
					test_e[acc].e[3] = (acc & 0x08) ? true : false;
				}

				// populate the N indices lists from the arrays built during InitTerrainIndices()
				for( int i=0; i<NUM_INDEX_LISTS; ++i ) {
					indices_tri_counts[i] = getIndices<unsigned short>(pl_short[i].v, test_e[i].e);
				}

				// iterate over each index list and optimize it
				for( int i=0; i<NUM_INDEX_LISTS; ++i ) {
					int tri_count = indices_tri_counts[i];
					VertexCacheOptimizerUShort vco;
					VertexCacheOptimizerUShort::Result res = vco.Optimize(&pl_short[i].v[0], tri_count);
					PiAssert(0 == res);
				}
			}

			// everything should be hunky-dory for setting up as OpenGL index buffers now.
			for( int i=0; i<NUM_INDEX_LISTS; ++i ) {
				glGenBuffersARB(1, &indices_list[i]);
				glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, indices_list[i]);
				glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned short)*indices_tri_counts[i]*3, &(pl_short[i].v[0]), GL_STATIC_DRAW);
				glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);
			}

			// default it to the last entry which uses the hi-res borders
			indices_vbo			= indices_list[NUM_INDEX_LISTS-1];
			indices_tri_count	= indices_tri_counts[NUM_INDEX_LISTS-1];

			// if enabled we create the threads which can build the terrain
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

		if (midIndices) {
			SAFE_DELETE_ARRAY(midIndices);
			for (int i=0; i<4; i++) {
				SAFE_DELETE_ARRAY(loEdgeIndices[i]);
				SAFE_DELETE_ARRAY(hiEdgeIndices[i]);
			}
		}
	}

	void UpdateVBOs() {
		m_needUpdateVBOs = true;
	}

	void _UpdateVBOs() {
		if (m_needUpdateVBOs) {
			if (!m_vbo) glGenBuffersARB(1, &m_vbo);
			m_needUpdateVBOs = false;
			glBindBufferARB(GL_ARRAY_BUFFER, m_vbo);
			glBufferDataARB(GL_ARRAY_BUFFER, sizeof(VBOVertex)*GEOPLATE_NUMVERTICES, 0, GL_DYNAMIC_DRAW);
			for (int i=0; i<GEOPLATE_NUMVERTICES; i++)
			{
				clipRadius = std::max(clipRadius, (vertices[i]-clipCentroid).Length());
				VBOVertex *pData = vbotemp + i;
				pData->x = float(vertices[i].x - clipCentroid.x);
				pData->y = float(vertices[i].y - clipCentroid.y);
				pData->z = float(vertices[i].z - clipCentroid.z);
				pData->nx = float(normals[i].x);
				pData->ny = float(normals[i].y);
				pData->nz = float(normals[i].z);
				pData->col[0] = static_cast<unsigned char>(Clamp(colors[i].x*255.0, 0.0, 255.0));
				pData->col[1] = static_cast<unsigned char>(Clamp(colors[i].y*255.0, 0.0, 255.0));
				pData->col[2] = static_cast<unsigned char>(Clamp(colors[i].z*255.0, 0.0, 255.0));
				pData->col[3] = 255;
			}
			glBufferDataARB(GL_ARRAY_BUFFER, sizeof(VBOVertex)*GEOPLATE_NUMVERTICES, vbotemp, GL_DYNAMIC_DRAW);
			glBindBufferARB(GL_ARRAY_BUFFER, 0);
		}
	}	
	/* not quite edge, since we share edge vertices so that would be
	 * fucking pointless. one position inwards. used to make edge normals
	 * for adjacent tiles */
	void GetEdgeMinusOneVerticesFlipped(int edge, vector3d *ev) {
		PiAssert(edge!=-1);
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
		PiAssert(edge!=-1);
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
		PiAssert(edge!=-1);
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
		for (int i=0; i<4; i++) {
			if (edgeFriend[i] == e) return i;
		}
		//abort();
		return -1;
	}


	void FixEdgeNormals(const int edge, const vector3d *ev) {
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
				const vector3d p = GetSurfacePointCyl(x*GEOPLATE_FRAC, 0.0, m_halfLen);
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
				const vector3d p = GetSurfacePointCyl(0.0, y*GEOPLATE_FRAC, m_halfLen);
				const double height = colors[y*GEOPLATE_EDGELEN].x;
				colors[y*GEOPLATE_EDGELEN] = geoRing->GetColor(p, height, norm);
			}
			break;
		}
	}

	int GetChildIdx(GeoPlate *child) {
		for (int i=0; i<4; i++) {
			if (kids[i] == child) return i;
		}
		abort();
		return -1;
	}
	
	void FixEdgeFromParentInterpolated(int edge) {
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

		SetEdge(this->vertices, edge, ev2);
		SetEdge(this->normals, edge, en2);
		SetEdge(this->colors, edge, ec2);
	}

	template <int corner>
	void MakeCornerNormal(vector3d *ev, vector3d *ev2) {
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
			} else {
				// what this does is gets the points from our own edge - 
				// - and just copies them, it is most certainly not correct,
				// - however it looks effective enough for the normals.
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

	// Generates full-detail vertices, and also non-edge normals and colors
	void GenerateMesh() {
		vector3d *vts = vertices;
		vector3d *col = colors;
		double xfrac;
		double zfrac = 0;
		vector3d axisPt(0.0, 0.0, 0.0);

		const vector3d topEnd(0.0, m_yoffset + m_halfLen, 0.0);		// vertices at top edge of circle
		const vector3d btmEnd(0.0, m_yoffset - m_halfLen, 0.0);		// vertices at bottom edge of circle
				
		for (int y=0; y<GEOPLATE_EDGELEN; ++y) {
			xfrac = 0;
			// point along z-axis by zfrac amount
			const vector3d axisPt = lerp( zfrac, btmEnd, topEnd );
			for (int x=0; x<GEOPLATE_EDGELEN; ++x) {
				// find point on _surface_ of the cylinder
				const vector3d p = GetSurfacePointCyl(xfrac, zfrac, m_halfLen);
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
			zfrac += GEOPLATE_FRAC;
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
		assert(e>(void*)0x000000FF);
		edgeFriend[edge] = e;
		vector3d ev[GEOPLATE_MAX_EDGELEN];
		int we_are = e->GetEdgeIdxOf(this);
		e->GetEdgeMinusOneVerticesFlipped(we_are, ev);
		
		vector3d axisPt(0.0, 0.0, 0.0);
		const vector3d topEnd(0.0, m_yoffset + m_halfLen, 0.0);		// vertices at top edge of circle
		const vector3d btmEnd(0.0, m_yoffset - m_halfLen, 0.0);	// vertices at bottom edge of circle

		/* now we have a valid edge, fix the edge vertices */
		if (edge == 0) {
			// point along z-axis by zfrac amount
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
				const double zfrac = y * GEOPLATE_FRAC;
				// point along z-axis by zfrac amount
				const vector3d axisPt = lerp( zfrac, btmEnd, topEnd );
				vector3d p = GetSurfacePointCyl(1.0, zfrac, m_halfLen);
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
			// point along z-axis by zfrac amount
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
				const double zfrac = y * GEOPLATE_FRAC;
				// point along z-axis by zfrac amount
				const vector3d axisPt = lerp( zfrac, btmEnd, topEnd );
				vector3d p = GetSurfacePointCyl(0, zfrac, m_halfLen);
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
		int idx = GetEdgeIdxOf(e);
		int we_are = e->GetEdgeIdxOf(this);
		if (!kids[0]) return;
		if (-1 == we_are || -1 == idx) return;
		// match e's new kids to our own... :/
		kids[idx]->OnEdgeFriendChanged(idx, e->kids[(we_are+1)%4]);
		kids[(idx+1)%4]->OnEdgeFriendChanged(idx, e->kids[we_are]);
	}
	
	void NotifyEdgeFriendDeleted(GeoPlate *e) {
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
		GeoPlate *e = edgeFriend[edge];
		if (!e) return 0;
		const int we_are = e->GetEdgeIdxOf(this);
		if (-1 == we_are ) return 0;
		// neighbour patch has not split yet (is at depth of this patch), so kids of this patch do
		// not have same detail level neighbours yet
		if (edge == kid) return e->kids[(we_are+1)%4];
		else return e->kids[we_are];
	}
	
	void Render(vector3d &campos, Plane planes[6]) {
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

			vector3d relpos = clipCentroid - campos;
			glPushMatrix();
			glTranslated(relpos.x, relpos.y, relpos.z);

			Pi::statSceneTris += 2*(GEOPLATE_EDGELEN-1)*(GEOPLATE_EDGELEN-1);
			glEnableClientState(GL_VERTEX_ARRAY);
			glEnableClientState(GL_NORMAL_ARRAY);
			glEnableClientState(GL_COLOR_ARRAY);

			// update the indices used for rendering
			updateIndexBufferId();

			glBindBufferARB(GL_ARRAY_BUFFER, m_vbo);
			glVertexPointer(3, GL_FLOAT, sizeof(VBOVertex), 0);
			glNormalPointer(GL_FLOAT, sizeof(VBOVertex), reinterpret_cast<void *>(3*sizeof(float)));
			glColorPointer(4, GL_UNSIGNED_BYTE, sizeof(VBOVertex), reinterpret_cast<void *>(6*sizeof(float)));
			glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, indices_vbo);
			glDrawElements(GL_TRIANGLES, indices_tri_count*3, GL_UNSIGNED_SHORT, 0);
			glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
			glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);

			glDisableClientState(GL_VERTEX_ARRAY);
			glDisableClientState(GL_NORMAL_ARRAY);
			glDisableClientState(GL_COLOR_ARRAY);
			glPopMatrix();
		}
	}

	void LODUpdate(vector3d &campos) {
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
		    ((campos - centroid).Length() < m_roughLength))) {
			canSplit = false;
		}
		// always split at first level
		//if (!parent) canSplit = true;

		bool canMerge = true;

		if (canSplit) {
			if (!kids[0]) {
				const vector3d halfVBE = vbe[0] + 0.5 * (vbe[1] - vbe[0]);

				const double halfAng = ang[0] + 0.5 * (ang[1] - ang[0]);
				const double newLength = (m_halfLen*0.5);
				const double yoffset0 = m_yoffset + newLength;
				const double yoffset1 = m_yoffset - newLength;
				GeoPlate *_kids[4];
				_kids[0] = new GeoPlate(newLength, halfVBE, vbe[1], halfAng, ang[1], yoffset1, m_depth+1, 0);
				_kids[1] = new GeoPlate(newLength, vbe[0], halfVBE, ang[0], halfAng, yoffset1, m_depth+1, 1);
				_kids[2] = new GeoPlate(newLength, vbe[0], halfVBE, ang[0], halfAng, yoffset0, m_depth+1, 2);
				_kids[3] = new GeoPlate(newLength, halfVBE, vbe[1], halfAng, ang[1], yoffset0, m_depth+1, 3);
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
unsigned short *GeoPlate::loEdgeIndices[4] = {0};
unsigned short *GeoPlate::hiEdgeIndices[4] = {0};
GLuint GeoPlate::indices_vbo = 0;
GLuint GeoPlate::indices_list[NUM_PLATE_INDICES] = {0};
GLuint GeoPlate::indices_tri_count = 0;
GLuint GeoPlate::indices_tri_counts[NUM_PLATE_INDICES] = {0};
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
}

void GeoRing::Init()
{
	s_geoRingSurfaceShader[0] = new GeoRingShader("geoRing", "#define NUM_LIGHTS 1\n");
	s_geoRingSurfaceShader[1] = new GeoRingShader("geoRing", "#define NUM_LIGHTS 2\n");
	s_geoRingSurfaceShader[2] = new GeoRingShader("geoRing", "#define NUM_LIGHTS 3\n");
	s_geoRingSurfaceShader[3] = new GeoRingShader("geoRing", "#define NUM_LIGHTS 4\n");
	s_geoRingSkyShader[0] = new GeoRingShader("geoRing_sky", "#define NUM_LIGHTS 1\n");
	s_geoRingSkyShader[1] = new GeoRingShader("geoRing_sky", "#define NUM_LIGHTS 2\n");
	s_geoRingSkyShader[2] = new GeoRingShader("geoRing_sky", "#define NUM_LIGHTS 3\n");
	s_geoRingSkyShader[3] = new GeoRingShader("geoRing_sky", "#define NUM_LIGHTS 4\n");
	s_allGeoRingsLock = SDL_CreateMutex();
	OnChangeDetailLevel();
#ifdef GEORING_USE_THREADING
	SDL_CreateThread(&GeoRing::UpdateLODThread, 0);
#endif /* GEORING_USE_THREADING */
	atexit(&_LockoutThreadsBeforeExit);
}

void GeoRing::OnChangeDetailLevel()
{
	SDL_mutexP(s_allGeoRingsLock);
	for(std::list<GeoRing*>::iterator i = s_allGeoRings.begin(); i != s_allGeoRings.end(); ++i) {
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
	for(std::list<GeoRing*>::iterator i = s_allGeoRings.begin(); i != s_allGeoRings.end(); ++i) {
		(*i)->BuildFirstPatches();
	}
	SDL_mutexV(s_allGeoRingsLock);
}

#define GEOSPHERE_TYPE	(m_sbody->type)

GeoRing::GeoRing(const SBody *body): m_style(body)
{
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
	SDL_mutexP(s_allGeoRingsLock);
	s_allGeoRings.remove(this);
	SDL_mutexV(s_allGeoRingsLock);

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
	DestroyVBOs();
	SDL_DestroyMutex(m_vbosToDestroyLock);
}

void GeoRing::AddVBOToDestroy(GLuint vbo)
{
	SDL_mutexP(m_vbosToDestroyLock);
	m_vbosToDestroy.push_back(vbo);
	SDL_mutexV(m_vbosToDestroyLock);
}

void GeoRing::DestroyVBOs()
{
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
	const int numSegments = (NUM_SEGMENTS * m_sbody->radius.ToDouble());

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
		m_plates.push_back( new GeoPlate(mRingWidth, points[i], points[i+1], angles[i], angles[i+1], 0.0, 0, 0) );
	}

	// set edge friends
	for (size_t i=0; i<m_plates.size(); i++) {
		m_plates[i]->SetGeoRingPtr( this );

		// The winding for these edges is not at all clear.
		// This is probably due to me adapting it for squares on the -
		// - outside of a spherised-cube to the inside of a cylinderised-lidless-cube!!!
		m_plates[i]->edgeFriend[0] = 0;
		m_plates[i]->edgeFriend[1] = FindGeoPlateByIndex(i-1);	// backward?
		m_plates[i]->edgeFriend[2] = 0;
		m_plates[i]->edgeFriend[3] = FindGeoPlateByIndex(i+1);	// forward?
	}

	// create mesh(es)
	for (size_t i=0; i<m_plates.size(); i++) m_plates[i]->GenerateMesh();
	for (size_t i=0; i<m_plates.size(); i++) m_plates[i]->GenerateEdgeNormalsAndColors();
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

static const float g_ambient[4] = { 0, 0, 0, 1.0 };

static void DrawAtmosphereSurface(const vector3d &campos, float rad)
{
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
	Plane planes[6];
	glPushMatrix();
 	glTranslated(-campos.x, -campos.y, -campos.z);
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
			shader->set_georingScale(scale);
			shader->set_georingAtmosTopRad(atmosRadius*radius/scale);
			shader->set_georingAtmosFogDensity(atmosDensity);
			shader->set_atmosColor(atmosCol.r, atmosCol.g, atmosCol.b, atmosCol.a);
			shader->set_georingCenter(center.x, center.y, center.z);
			
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
		shader->set_georingScale(scale);
		shader->set_georingAtmosTopRad(atmosRadius*radius/scale);
		shader->set_georingAtmosFogDensity(atmosDensity);
		shader->set_atmosColor(atmosCol.r, atmosCol.g, atmosCol.b, atmosCol.a);
		shader->set_georingCenter(center.x, center.y, center.z);
	}

	glPopMatrix();

	if (0==m_plates.size()) {
		BuildFirstPatches();
	}

	const float black[4] = { 0,0,0,0 };
	float ambient[4] = { 0.1, 0.1, 0.1, 1.0 };

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

	{
		for (size_t i=0; i<m_wallInner.size(); ++i) {
			m_wallInner[i]->Render(campos, planes);
		}
	}

	{
		for (size_t i=0; i<m_wallOuter.size(); ++i) {
			m_wallOuter[i]->Render(campos, planes);
		}
	}
	
	{
		for (size_t i=0; i<m_hull.size(); ++i) {
			m_hull[i]->Render(campos, planes);
		}
	}
	
	{
		for (size_t i=0; i<m_plates.size(); ++i) {
			m_plates[i]->Render(campos, planes);
		}
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
