/* Copyright (C) Eric Lengyel, 2001. 
 * All rights reserved worldwide.
 *
 * This software is provided "as is" without express or implied
 * warranties. You may freely copy and compile this source into
 * applications you distribute provided that the copyright text
 * below is included in the resulting source code, for example:
 * "Portions Copyright (C) Eric Lengyel, 2001"
 */

#include "Lengyel.h"

float Dot(const vector4f &a, const vector3f &b) { 
	return a.x*b.x + a.y*b.y + a.z*b.z; 
}

Decal::Decal(const vector3f& center, const vector3f& normal, const vector3f& tangent, float width, float height, float depth)
{
	decalCenter = center;
	decalNormal = normal;
	
	vector3f binormal = normal.Cross(tangent);
	
	// Calculate boundary planes
	float d = center.Dot(tangent);
	leftPlane = vector4f(tangent.x, tangent.y, tangent.z, width * 0.5F - d);
	rightPlane = vector4f(-tangent.x, -tangent.y, -tangent.z, width * 0.5F + d);
	
	d = center.Dot(binormal);
	bottomPlane = vector4f(binormal.x, binormal.y, binormal.z, height * 0.5F - d);
	topPlane = vector4f(-binormal.x, -binormal.y, -binormal.z, height * 0.5F + d);
	
	d = center.Dot(normal);
	frontPlane = vector4f(-normal.x, -normal.y, -normal.z, depth + d);
	backPlane = vector4f(normal.x, normal.y, normal.z, depth - d);
	
	// Begin with empty mesh
	decalVertexCount = 0;
	decalTriangleCount = 0;
	
	// Add this point, determine which surfaces may be affected by this decal
	// and call ClipMesh() for each one.
	//ClipMesh(long triangleCount, const Triangle *triangle, const vector3f *vertex, const vector3f *normal)
	
	// Assign texture mapping coordinates
	float one_over_w = 1.0F / width;
	float one_over_h = 1.0F / height;
	for (long a = 0; a < decalVertexCount; a++)
	{
		vector3f v = vertexArray[a] - center;
		float s = v.Dot(tangent) * one_over_w + 0.5F;
		float t = v.Dot(binormal) * one_over_h + 0.5F;
		texcoordArray[a] = vector2f(s, t);
	}
}

bool Decal::AddPolygon(long vertexCount, const vector3f *vertex, const vector3f *normal)
{
	long count = decalVertexCount;
	if (count + vertexCount >= maxDecalVertices) return (false);
	
	// Add polygon as a triangle fan
	Triangle *triangle = triangleArray + decalTriangleCount;
	decalTriangleCount += vertexCount - 2;
	for (long a = 2; a < vertexCount; a++)
	{
		triangle->index[0] = (unsigned short) count;
		triangle->index[1] = (unsigned short) (count + a - 1);
		triangle->index[2] = (unsigned short) (count + a);
		triangle++;
	}
	
	// Assign vertex colors
	float f = 1.0F / (1.0F - decalEpsilon);
	for (long b = 0; b < vertexCount; b++)
	{
		vertexArray[count] = vertex[b];
		const vector3f& n = normal[b];
		float alpha = (decalNormal.Dot(n) / n.Length() - decalEpsilon) * f;
		colorArray[count] = ColorRGBA(1.0F, 1.0F, 1.0F, (alpha > 0.0F) ? alpha : 0.0F);
		count++;
	}
	
	decalVertexCount = count;
	return (true);
}

void Decal::ClipMesh(long triangleCount, const Triangle *triangle, const vector3f *vertex, const vector3f *normal)
{
	vector3f		newVertex[9];
	vector3f		newNormal[9];
	
	// Clip one triangle at a time
	for (long a = 0; a < triangleCount; a++)
	{
		long i1 = triangle->index[0];
		long i2 = triangle->index[1];
		long i3 = triangle->index[2];
		
		const vector3f& v1 = vertex[i1];
		const vector3f& v2 = vertex[i2];
		const vector3f& v3 = vertex[i3];
		
		vector3f cross = (v2 - v1).Cross(v3 - v1);
		if (decalNormal.Dot(cross) > decalEpsilon * cross.Length())
		{
			newVertex[0] = v1;
			newVertex[1] = v2;
			newVertex[2] = v3;
			
			newNormal[0] = normal[i1];
			newNormal[1] = normal[i2];
			newNormal[2] = normal[i3];
			
			long count = ClipPolygon(3, newVertex, newNormal, newVertex, newNormal);
			if ((count != 0) && (!AddPolygon(count, newVertex, newNormal))) break;
		}
		
		triangle++;
	}
}

long Decal::ClipPolygon(long vertexCount, const vector3f *vertex, const vector3f *normal, vector3f *newVertex, vector3f *newNormal) const
{
	vector3f		tempVertex[9];
	vector3f		tempNormal[9];
	
	// Clip against all six planes
	long count = ClipPolygonAgainstPlane(leftPlane, vertexCount, vertex, normal, tempVertex, tempNormal);
	if (count != 0)
	{
		count = ClipPolygonAgainstPlane(rightPlane, count, tempVertex, tempNormal, newVertex, newNormal);
		if (count != 0)
		{
			count = ClipPolygonAgainstPlane(bottomPlane, count, newVertex, newNormal, tempVertex, tempNormal);
			if (count != 0)
			{
				count = ClipPolygonAgainstPlane(topPlane, count, tempVertex, tempNormal, newVertex, newNormal);
				if (count != 0)
				{
					count = ClipPolygonAgainstPlane(backPlane, count, newVertex, newNormal, tempVertex, tempNormal);
					if (count != 0)
					{
						count = ClipPolygonAgainstPlane(frontPlane, count, tempVertex, tempNormal, newVertex, newNormal);
					}
				}
			}
		}
	}
	
	return (count);
}

long Decal::ClipPolygonAgainstPlane(const vector4f& plane, long vertexCount, const vector3f *vertex, const vector3f *normal, vector3f *newVertex, vector3f *newNormal)
{
	bool	negative[10];
	
	// Classify vertices
	long negativeCount = 0;
	for (long a = 0; a < vertexCount; a++)
	{
		bool neg = (Dot(plane, vertex[a]) < 0.0F);
		negative[a] = neg;
		negativeCount += neg;
	}
	
	// Discard this polygon if it's completely culled
	if (negativeCount == vertexCount) return (0);
	
	long count = 0;
	for (long b = 0; b < vertexCount; b++)
	{
		// c is the index of the previous vertex
		long c = (b != 0) ? b - 1 : vertexCount - 1;
		
		if (negative[b])
		{
			if (!negative[c])
			{
				// Current vertex is on negative side of plane,
				// but previous vertex is on positive side.
				const vector3f& v1 = vertex[c];
				const vector3f& v2 = vertex[b];
				float t = Dot(plane, v1) / (plane.x * (v1.x - v2.x) + plane.y * (v1.y - v2.y) + plane.z * (v1.z - v2.z));
				newVertex[count] = v1 * (1.0F - t) + v2 * t;
				
				const vector3f& n1 = normal[c];
				const vector3f& n2 = normal[b];
				newNormal[count] = n1 * (1.0F - t) + n2 * t;
				count++;
			}
		}
		else
		{
			if (negative[c])
			{
				// Current vertex is on positive side of plane,
				// but previous vertex is on negative side.
				const vector3f& v1 = vertex[b];
				const vector3f& v2 = vertex[c];
				float t = Dot(plane, v1) / (plane.x * (v1.x - v2.x) + plane.y * (v1.y - v2.y) + plane.z * (v1.z - v2.z));
				newVertex[count] = v1 * (1.0F - t) + v2 * t;
				
				const vector3f& n1 = normal[b];
				const vector3f& n2 = normal[c];
				newNormal[count] = n1 * (1.0F - t) + n2 * t;
				count++;
			}
			
			// Include current vertex
			newVertex[count] = vertex[b];
			newNormal[count] = normal[b];
			count++;
		}
	}
	
	// Return number of vertices in clipped polygon
	return (count);
}
