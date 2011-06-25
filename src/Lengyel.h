/* Copyright (C) Eric Lengyel, 2001. 
 * All rights reserved worldwide.
 *
 * This software is provided "as is" without express or implied
 * warranties. You may freely copy and compile this source into
 * applications you distribute provided that the copyright text
 * below is included in the resulting source code, for example:
 * "Portions Copyright (C) Eric Lengyel, 2001"
 */
 
 #ifndef __LENGYEL_H__
 #define __LENGYEL_H__

#include "vector2.h"
#include "vector3.h"
#include "vector4.h"


const long maxDecalVertices = 256;
const float decalEpsilon = 0.25F;


struct ColorRGBA
{
	float		red;
	float		green;
	float		blue;
	float		alpha;
	
	ColorRGBA() {}
	
	ColorRGBA(float r, float g, float b, float a)
	{
		red = r;
		green = g;
		blue = b;
		alpha = a;
	}
};

struct Triangle
{
	unsigned short	index[3];
};


inline float DotProduct(const vector4f& p, const vector3f& q)
{
	return (p.x * q.x + p.y * q.y + p.z * q.z + p.w);
}


class Decal
{
	private:
		
		vector3f		decalCenter;
		vector3f		decalNormal;
		
		vector4f		leftPlane;
		vector4f		rightPlane;
		vector4f		bottomPlane;
		vector4f		topPlane;
		vector4f		frontPlane;
		vector4f		backPlane;
		
		long		decalVertexCount;
		long		decalTriangleCount;
		
		vector3f	vertexArray[maxDecalVertices];
		vector2f	texcoordArray[maxDecalVertices];
		ColorRGBA	colorArray[maxDecalVertices];
		Triangle	triangleArray[maxDecalVertices];
		
		bool AddPolygon(long vertexCount, const vector3f *vertex, const vector3f *normal);
		void ClipMesh(long triangleCount, const Triangle *triangle, const vector3f *vertex, const vector3f *normal);
		long ClipPolygon(long vertexCount, const vector3f *vertex, const vector3f *normal, vector3f *newVertex, vector3f *newNormal) const;
		static long ClipPolygonAgainstPlane(const vector4f& plane, long vertexCount, const vector3f *vertex, const vector3f *normal, vector3f *newVertex, vector3f *newNormal);
	
	public:
		
		Decal(const vector3f& center, const vector3f& normal, const vector3f& tangent, float width, float height, float depth);
};

#endif // __LENGYEL_H__