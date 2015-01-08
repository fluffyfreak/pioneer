//
// Marching Cubes Example Program 
// by Cory Bloyd (corysama@yahoo.com)
//
// A simple, portable and complete implementation of the Marching Cubes
// and Marching Tetrahedrons algorithms in a single source file.
// There are many ways that this code could be made faster, but the 
// intent is for the code to be easy to understand.
//
// For a description of the algorithm go to
// http://astronomy.swin.edu.au/pbourke/modelling/polygonise/
//
// This code is public domain.
//

#include "MarchingGeom.h"

#include "graphics/VertexArray.h"
#include "graphics/VertexBuffer.h"
#include "graphics/RenderState.h"
#include "graphics/Material.h"

// include and use the Marching Cube & Tetrahedron data
#include "MarchingGeomData.h"
using namespace MarchingGeomData;

int     iDataSetSize =  16;
float   fStepSize = 1.0 / iDataSetSize;
float   fTargetValue = 48.0;
float   fTime = 0.0;
#define NUM_SRC_POINTS 3
vector3f  sSourcePoint[NUM_SRC_POINTS];
bool bSpin = false;
bool bMove = true;
bool bLight = true;

float mul = 10.0f;
float scaleRes = 60.0f;

void vPrintHelp();
void vSetTime(float fTime);

float fSample1(const float fX, const float fY, const float fZ);
float fSample2(const float fX, const float fY, const float fZ);
float fSample3(const float fX, const float fY, const float fZ);
float fSample4(const float fX, const float fY, const float fZ);
float(*fSample)(const float x, const float y, const float z) = fSample1;

void vMarchingCubes(Graphics::VertexArray&);
void vMarchCube1(Graphics::VertexArray&, float fX, float fY, float fZ, float fScale);
void vMarchCube2(Graphics::VertexArray&, float fX, float fY, float fZ, float fScale);
void(*vMarchCube)(Graphics::VertexArray&, float x, float y, float z, float fScale) = vMarchCube1;

using namespace Graphics;

MarchingGeometry::MarchingGeometry(Graphics::Renderer* renderer) : m_renderer(renderer), m_fTime(-0.025f)
{
}

void MarchingGeometry::Init()
{
	Graphics::RenderStateDesc rsd;
	m_renderState = m_renderer->CreateRenderState(rsd);

	MaterialDescriptor desc;
	desc.vertexColors = true;
	desc.lighting = true;
	m_material.Reset(m_renderer->CreateMaterial(desc));
	m_material->diffuse = Color::WHITE;
	m_material->emissive = Color::WHITE;
	m_material->specular = Color::WHITE;

	Update();
}

void MarchingGeometry::Draw()
{
	if (m_vb.Valid()) {
		m_renderer->DrawBuffer(m_vb.Get(), m_renderState, m_material.Get());
	}
}

static bool bUpdateAnyway = true;
void MarchingGeometry::Update()
{
	m_fTime += 0.0125f;
	vSetTime(m_fTime);

	if (m_vb.Valid() && fSample == fSample4 && !bUpdateAnyway) {
		return;
	}

	// Create skybox geometry
	Graphics::VertexArray mc(ATTRIB_POSITION | ATTRIB_NORMAL | ATTRIB_DIFFUSE);
	vMarchingCubes(mc);

	if (!m_vb.Get() || (mc.GetNumVerts() != m_vb->GetVertexCount()))
	{
		//create buffer and upload data
		Graphics::VertexBufferDesc vbd;
		vbd.attrib[0].semantic = Graphics::ATTRIB_POSITION;
		vbd.attrib[0].format = Graphics::ATTRIB_FORMAT_FLOAT3;
		vbd.attrib[1].semantic = Graphics::ATTRIB_NORMAL;
		vbd.attrib[1].format = Graphics::ATTRIB_FORMAT_FLOAT3;
		vbd.attrib[2].semantic = Graphics::ATTRIB_DIFFUSE;
		vbd.attrib[2].format = Graphics::ATTRIB_FORMAT_UBYTE4;
		vbd.numVertices = mc.GetNumVerts();
		vbd.usage = Graphics::BUFFER_USAGE_DYNAMIC;

		m_vb.Reset(m_renderer->CreateVertexBuffer(vbd));
	}

	m_vb->Populate(mc);
}

//fGetOffset finds the approximate point of intersection of the surface
// between two points with the values fValue1 and fValue2
float fGetOffset(const float fValue1, const float fValue2, const float fValueDesired)
{
	const float fDelta = fValue2 - fValue1;

	if(fDelta == 0.0f)
	{
		return 0.5f;
	}
	return (fValueDesired - fValue1)/fDelta;
}


//vGetColor generates a color from a given position and normal of a point
void vGetColor(vector3f &rfColor, const vector3f &rfPosition, const vector3f &rfNormal)
{
	const float x = rfNormal.x;
	const float y = rfNormal.y;
	const float z = rfNormal.z;
	rfColor.x = (x > 0.0 ? x : 0.0) + (y < 0.0 ? -0.5*y : 0.0) + (z < 0.0 ? -0.5*z : 0.0);
	rfColor.y = (y > 0.0 ? y : 0.0) + (z < 0.0 ? -0.5*z : 0.0) + (x < 0.0 ? -0.5*x : 0.0);
	rfColor.z = (z > 0.0 ? z : 0.0) + (x < 0.0 ? -0.5*x : 0.0) + (y < 0.0 ? -0.5*y : 0.0);
}

void vNormalizeVector(vector3f &rfVectorResult, const vector3f &rfVectorSource)
{
	const float fOldLength = sqrtf((rfVectorSource.x * rfVectorSource.x) +
						(rfVectorSource.y * rfVectorSource.y) +
						(rfVectorSource.z * rfVectorSource.z) );

	if(fOldLength == 0.0)
	{
		rfVectorResult.x = rfVectorSource.x;
		rfVectorResult.y = rfVectorSource.y;
		rfVectorResult.z = rfVectorSource.z;
	}
	else
	{
		const float fScale = 1.0 / fOldLength;
		rfVectorResult.x = rfVectorSource.x*fScale;
		rfVectorResult.y = rfVectorSource.y*fScale;
		rfVectorResult.z = rfVectorSource.z*fScale;
	}
}


//Generate a sample data set.  fSample1(), fSample2() and fSample3() define three scalar fields whose
// values vary by the X,Y and Z coordinates and by the fTime value set by vSetTime()
void vSetTime(const float fNewTime)
{
	float fOffset;
	int iSourceNum;

	for (iSourceNum = 0; iSourceNum < NUM_SRC_POINTS; iSourceNum++)
	{
		sSourcePoint[iSourceNum].x = 0.5;
		sSourcePoint[iSourceNum].y = 0.5;
		sSourcePoint[iSourceNum].z = 0.5;
	}

	fTime = fNewTime;
	fOffset = 1.0 + sinf(fTime);
	sSourcePoint[0].x *= fOffset;
	sSourcePoint[1].y *= fOffset;
	sSourcePoint[2].z *= fOffset;
}

//fSample1 finds the distance of (x, y, z) from three moving points
float fSample1(const float x, const float y, const float z)
{
	double fResult = 0.0;
	double fDx, fDy, fDz;
	fDx = x - sSourcePoint[0].x;
	fDy = y - sSourcePoint[0].y;
	fDz = z - sSourcePoint[0].z;
	fResult += 0.5/(fDx*fDx + fDy*fDy + fDz*fDz);

	fDx = x - sSourcePoint[1].x;
	fDy = y - sSourcePoint[1].y;
	fDz = z - sSourcePoint[1].z;
	fResult += 1.0/(fDx*fDx + fDy*fDy + fDz*fDz);

	fDx = x - sSourcePoint[2].x;
	fDy = y - sSourcePoint[2].y;
	fDz = z - sSourcePoint[2].z;
	fResult += 1.5/(fDx*fDx + fDy*fDy + fDz*fDz);

	return fResult;
}

//fSample2 finds the distance of (x, y, z) from three moving lines
float fSample2(const float x, const float y, const float z)
{
	double fResult = 0.0;
	double fDx, fDy, fDz;
	fDx = x - sSourcePoint[0].x;
	fDy = y - sSourcePoint[0].y;
	fResult += 0.5/(fDx*fDx + fDy*fDy);

	fDx = x - sSourcePoint[1].x;
	fDz = z - sSourcePoint[1].z;
	fResult += 0.75/(fDx*fDx + fDz*fDz);

	fDy = y - sSourcePoint[2].y;
	fDz = z - sSourcePoint[2].z;
	fResult += 1.0/(fDy*fDy + fDz*fDz);

	return fResult;
}


//fSample3 defines a height field by plugging the distance from the center into the sin and cos functions
float fSample3(const float x, const float y, const float z)
{
	float fHeight = 20.0*(fTime + sqrt((0.5-x)*(0.5-x) + (0.5-y)*(0.5-y)));
	fHeight = 1.5 + 0.1*(sinf(fHeight) + cosf(fHeight));
	return (fHeight - z)*50.0;
}

namespace Worley {
	inline vector3f floor(const vector3f& a) { 
		return vector3f(
			::floor(a.x), 
			::floor(a.y), 
			::floor(a.z)); 
	}
	inline float fract(const float x) { return x - ::floor(x); }
	inline float rnd(float x) { return fract(1000.*sin(234.56*x)); }
	vector3f rnd3(float x) { return vector3f(rnd(x), rnd(x + .1), rnd(x + .2)); }
	inline float hash(float x, float y, float z) { return (x + 432.432*y - 1178.65*z); }
	inline float hash(vector3f v) { return v.Dot(vector3f(1., 32.432, -1178.65)); }

	class vector4f {
	public:
		float x, y, z, w;
		vector4f(float _x, float _y, float _z, float _w) : x(_x), y(_y), z(_z), w(_w) {}
	};

	vector4f WorleyVec4(vector3f uvw)
	{
		const vector3f uvwi = floor(uvw);							// cell coords
		float dmin = 1e9, d2min = 1e9, nmin = -1.;

		for (int i = -1; i <= 1; i++) {						// visit neighborhood
			for (int j = -1; j <= 1; j++) {						// to find the closest point
				for (int k = -1; k <= 1; k++)
				{
					const vector3f c = uvwi + vector3f(float(i), float(j), float(k));	// neighbor cells
					const float n = hash(c);	 										// cell ID
					const vector3f p = c + rnd3(n + .1);								// random point in cell
					const float d = (p - uvw).Length();									// dist to point
					if (d < dmin) {
						d2min = dmin; dmin = d; nmin = n;
					}		// 2 closest dists
					else if (d < d2min) {
						d2min = d;
					}
				}
			}
		}
		return vector4f(dmin, d2min, d2min - dmin, nmin);			// 2 closest dists + closest ID
	}

	float Worleyf(vector3f uvw) 
	{
		const vector3f uvwi = floor(uvw);							// cell coords
		float dmin = 1e9, d2min = 1e9, nmin = -1.;

		for (int i = -1; i <= 1; i++) {						// visit neighborhood
			for (int j = -1; j <= 1; j++) {						// to find the closest point
				for (int k = -1; k <= 1; k++)
				{
					const vector3f c = uvwi + vector3f(float(i), float(j), float(k));	// neighbor cells
					const float n = hash(c);	 										// cell ID
					const vector3f p = c + rnd3(n + .1);								// random point in cell
					const float d = (p - uvw).Length();									// dist to point
					dmin = std::min(d, dmin);
				}
			}
		}
		return dmin;
	}
}

//fSample4 defines a height field by plugging the distance from the center into the sin and cos functions
float fSample4(const float x, const float y, const float z)
{
	//Worley::vector4f v4 = WorleyVec4::Worley(vector3f(x, y, z) * mul);
	const float redx = Worley::Worleyf(vector3f(x * mul, y * mul, z * mul));
	return redx * scaleRes;
}


//vGetNormal() finds the gradient of the scalar field at a point
//This gradient can be used as a very accurate vertx normal for lighting calculations
void vGetNormal(vector3f &rfNormal, const float x, const float y, const float z)
{
	rfNormal.x = fSample(x-0.01, y, z) - fSample(x+0.01, y, z);
	rfNormal.y = fSample(x, y-0.01, z) - fSample(x, y+0.01, z);
	rfNormal.z = fSample(x, y, z-0.01) - fSample(x, y, z+0.01);
	vNormalizeVector(rfNormal, rfNormal);
}


//vMarchCube1 performs the Marching Cubes algorithm on a single cube
void vMarchCube1(Graphics::VertexArray& va, const float x, const float y, const float z, const float fScale)
{
	int iCorner, iVertex, iVertexTest, iEdge, iTriangle, iFlagIndex, iEdgeFlags;
	float fOffset;
	vector3f sColor;
	float afCubeValue[8];
	vector3f asEdgeVertex[12];
	vector3f asEdgeNorm[12];

	//Make a local copy of the values at the cube's corners
	for(iVertex = 0; iVertex < 8; iVertex++)
	{
		afCubeValue[iVertex] = fSample(x + a2fVertexOffset[iVertex][0]*fScale,
		y + a2fVertexOffset[iVertex][1]*fScale,
		z + a2fVertexOffset[iVertex][2]*fScale);
	}

	//Find which vertices are inside of the surface and which are outside
	iFlagIndex = 0;
	for(iVertexTest = 0; iVertexTest < 8; iVertexTest++)
	{
		if (afCubeValue[iVertexTest] <= fTargetValue) {
			iFlagIndex |= 1 << iVertexTest;
		}
	}

	//Find which edges are intersected by the surface
	iEdgeFlags = aiCubeEdgeFlags[iFlagIndex];

	//If the cube is entirely inside or outside of the surface, then there will be no intersections
	if(iEdgeFlags == 0) 
	{
		return;
	}

	//Find the point of intersection of the surface with each edge
	//Then find the normal to the surface at those points
	for(iEdge = 0; iEdge < 12; iEdge++)
	{
		//if there is an intersection on this edge
		if(iEdgeFlags & (1<<iEdge))
		{
			fOffset = fGetOffset(afCubeValue[ a2iEdgeConnection[iEdge][0] ], 
			afCubeValue[ a2iEdgeConnection[iEdge][1] ], fTargetValue);

			asEdgeVertex[iEdge].x = x + (a2fVertexOffset[ a2iEdgeConnection[iEdge][0] ][0]  +  fOffset * a2fEdgeDirection[iEdge][0]) * fScale;
			asEdgeVertex[iEdge].y = y + (a2fVertexOffset[ a2iEdgeConnection[iEdge][0] ][1]  +  fOffset * a2fEdgeDirection[iEdge][1]) * fScale;
			asEdgeVertex[iEdge].z = z + (a2fVertexOffset[ a2iEdgeConnection[iEdge][0] ][2]  +  fOffset * a2fEdgeDirection[iEdge][2]) * fScale;

			vGetNormal(asEdgeNorm[iEdge], asEdgeVertex[iEdge].x, asEdgeVertex[iEdge].y, asEdgeVertex[iEdge].z);
		}
	}

	//Draw the triangles that were found.  There can be up to five per cube
	for(iTriangle = 0; iTriangle < 5; iTriangle++)
	{
		if(a2iTriangleConnectionTable[iFlagIndex][3*iTriangle] < 0)
			break;

		for(iCorner = 0; iCorner < 3; iCorner++)
		{
			iVertex = a2iTriangleConnectionTable[iFlagIndex][3*iTriangle+iCorner];

			vGetColor(sColor, asEdgeVertex[iVertex], asEdgeNorm[iVertex]);
			va.Add(asEdgeVertex[iVertex], Color4ub(Color4f(sColor.x, sColor.y, sColor.z)), asEdgeNorm[iVertex]);
		}
	}
}

//vMarchTetrahedron performs the Marching Tetrahedrons algorithm on a single tetrahedron
void vMarchTetrahedron(Graphics::VertexArray& va, const vector3f *pasTetrahedronPosition, const float *pafTetrahedronValue)
{
	int iEdge, iVert0, iVert1, iEdgeFlags, iTriangle, iCorner, iVertex, iFlagIndex = 0;
	float fOffset, fInvOffset, fValue = 0.0;
	vector3f asEdgeVertex[6];
	vector3f asEdgeNorm[6];
	vector3f sColor;

	//Find which vertices are inside of the surface and which are outside
	for(iVertex = 0; iVertex < 4; iVertex++)
	{
		if (pafTetrahedronValue[iVertex] <= fTargetValue) {
			iFlagIndex |= 1 << iVertex;
		}
	}

	//Find which edges are intersected by the surface
	iEdgeFlags = aiTetrahedronEdgeFlags[iFlagIndex];

	//If the tetrahedron is entirely inside or outside of the surface, then there will be no intersections
	if(iEdgeFlags == 0)
	{
		return;
	}

	//Find the point of intersection of the surface with each edge
	// Then find the normal to the surface at those points
	for(iEdge = 0; iEdge < 6; iEdge++)
	{
		//if there is an intersection on this edge
		if(iEdgeFlags & (1<<iEdge))
		{
			iVert0 = a2iTetrahedronEdgeConnection[iEdge][0];
			iVert1 = a2iTetrahedronEdgeConnection[iEdge][1];
			fOffset = fGetOffset(pafTetrahedronValue[iVert0], pafTetrahedronValue[iVert1], fTargetValue);
			fInvOffset = 1.0 - fOffset;

			asEdgeVertex[iEdge].x = fInvOffset*pasTetrahedronPosition[iVert0].x  +  fOffset*pasTetrahedronPosition[iVert1].x;
			asEdgeVertex[iEdge].y = fInvOffset*pasTetrahedronPosition[iVert0].y  +  fOffset*pasTetrahedronPosition[iVert1].y;
			asEdgeVertex[iEdge].z = fInvOffset*pasTetrahedronPosition[iVert0].z  +  fOffset*pasTetrahedronPosition[iVert1].z;
			
			vGetNormal(asEdgeNorm[iEdge], asEdgeVertex[iEdge].x, asEdgeVertex[iEdge].y, asEdgeVertex[iEdge].z);
		}
	}

	//Draw the triangles that were found.  There can be up to 2 per tetrahedron
	for(iTriangle = 0; iTriangle < 2; iTriangle++)
	{
		if(a2iTetrahedronTriangles[iFlagIndex][3*iTriangle] < 0)
			break;

		for(iCorner = 0; iCorner < 3; iCorner++)
		{
			iVertex = a2iTetrahedronTriangles[iFlagIndex][3*iTriangle+iCorner];

			vGetColor(sColor, asEdgeVertex[iVertex], asEdgeNorm[iVertex]);
			
			va.Add(asEdgeVertex[iVertex], Color4ub(Color4f(sColor.x, sColor.y, sColor.z)), asEdgeNorm[iVertex]);
		}
	}
}



//vMarchCube2 performs the Marching Tetrahedrons algorithm on a single cube by making six calls to vMarchTetrahedron
void vMarchCube2(Graphics::VertexArray& va, const float x, const float y, const float z, const float fScale)
{
	int iVertex, iTetrahedron, iVertexInACube;
	vector3f asCubePosition[8];
	float  afCubeValue[8];
	vector3f asTetrahedronPosition[4];
	float  afTetrahedronValue[4];

	//Make a local copy of the cube's corner positions
	for(iVertex = 0; iVertex < 8; iVertex++)
	{
		asCubePosition[iVertex].x = x + a2fVertexOffset[iVertex][0]*fScale;
		asCubePosition[iVertex].y = y + a2fVertexOffset[iVertex][1]*fScale;
		asCubePosition[iVertex].z = z + a2fVertexOffset[iVertex][2]*fScale;
	}

	//Make a local copy of the cube's corner values
	for(iVertex = 0; iVertex < 8; iVertex++)
	{
		afCubeValue[iVertex] = fSample( asCubePosition[iVertex].x,
										asCubePosition[iVertex].y,
										asCubePosition[iVertex].z);
	}

	for(iTetrahedron = 0; iTetrahedron < 6; iTetrahedron++)
	{
		for(iVertex = 0; iVertex < 4; iVertex++)
		{
			iVertexInACube = a2iTetrahedronsInACube[iTetrahedron][iVertex];
			asTetrahedronPosition[iVertex].x = asCubePosition[iVertexInACube].x;
			asTetrahedronPosition[iVertex].y = asCubePosition[iVertexInACube].y;
			asTetrahedronPosition[iVertex].z = asCubePosition[iVertexInACube].z;
			afTetrahedronValue[iVertex] = afCubeValue[iVertexInACube];
		}
		vMarchTetrahedron(va, asTetrahedronPosition, afTetrahedronValue);
	}
}


//vMarchingCubes iterates over the entire dataset, calling vMarchCube on each cube
void vMarchingCubes(Graphics::VertexArray& va)
{
	int iX, iY, iZ;
	for (iX = 0; iX < iDataSetSize; iX++) {
		for (iY = 0; iY < iDataSetSize; iY++) {
			for (iZ = 0; iZ < iDataSetSize; iZ++)
			{
				vMarchCube(va, iX*fStepSize, iY*fStepSize, iZ*fStepSize, fStepSize);
			}
		}
	}
}
