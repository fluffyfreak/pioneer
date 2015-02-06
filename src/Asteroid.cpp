// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Asteroid.h"
#include "graphics/Texture.h"
#include "collider/Weld.h"
#include "gameconsts.h"
#include "utils.h"
#include <set>

using namespace Graphics;
using namespace Graphics::Drawables;

//------------------------------------------------------------

namespace
{
	struct PosNormTangentUVVert {
		vector3f pos;
		vector3f norm;
		vector2f uv0;
		vector3f tangent;
		bool operator==(const PosNormTangentUVVert &a) const {
			return (pos == a.pos) && (norm == a.norm) && (uv0 == a.uv0) && (tangent == a.tangent);
		}
	};

	// Lengyel, Eric. “Computing Tangent Space Basis Vectors for an Arbitrary Mesh”. 
	// Terathon Software 3D Graphics Library, 2001. 
	// http://www.terathon.com/code/tangent.html
	static void CalculateTangentArray(
		std::vector<PosNormTangentUVVert> &vertices,
		const std::vector<Uint16> &indices)
	{
		PROFILE_SCOPED()
		const size_t vertexCount = vertices.size();
		vector3f *tan1 = new vector3f[vertexCount * 2];
		vector3f *tan2 = tan1 + vertexCount;
		//ZeroMemory(tan1, vertexCount * sizeof(vector3f) * 2);
		memset(tan1, 0, vertexCount * sizeof(vector3f) * 2);

		const size_t triangleCount = indices.size() / 3;
		for (size_t a = 0; a < triangleCount; a++)
		{
			const size_t i1 = indices[(a * 3) + 0];
			const size_t i2 = indices[(a * 3) + 1];
			const size_t i3 = indices[(a * 3) + 2];

			const vector3f& v1 = vertices[i1].pos;
			const vector3f& v2 = vertices[i2].pos;
			const vector3f& v3 = vertices[i3].pos;

			const vector2f& w1 = vertices[i1].uv0;
			const vector2f& w2 = vertices[i2].uv0;
			const vector2f& w3 = vertices[i3].uv0;

			const float x1 = v2.x - v1.x;
			const float x2 = v3.x - v1.x;
			const float y1 = v2.y - v1.y;
			const float y2 = v3.y - v1.y;
			const float z1 = v2.z - v1.z;
			const float z2 = v3.z - v1.z;

			const float s1 = w2.x - w1.x;
			const float s2 = w3.x - w1.x;
			const float t1 = w2.y - w1.y;
			const float t2 = w3.y - w1.y;

			// handle the divide by zero case!
			const float div(s1 * t2 - s2 * t1);
			const float r = (is_equal_exact(div,0.0f)) ? 1.0 : (1.0F / div);

			const vector3f sdir(
				(t2 * x1 - t1 * x2) * r,
				(t2 * y1 - t1 * y2) * r,
				(t2 * z1 - t1 * z2) * r);
			const vector3f tdir(
				(s1 * x2 - s2 * x1) * r,
				(s1 * y2 - s2 * y1) * r,
				(s1 * z2 - s2 * z1) * r);

			tan1[i1] += sdir;
			tan1[i2] += sdir;
			tan1[i3] += sdir;

			tan2[i1] += tdir;
			tan2[i2] += tdir;
			tan2[i3] += tdir;
		}

		for (size_t a = 0; a < vertexCount; a++)
		{
			const vector3f& n = vertices[a].norm;
			const vector3f& t = tan1[a];

			// Gram-Schmidt orthogonalize
			vertices[a].tangent = (t - n * n.Dot(t)).Normalized();
		}

		delete[] tan1;
	}

	void CalculatePerVertexNormals(std::vector<PosNormTangentUVVert> &vertices, std::vector<Uint16> &indices, std::vector < std::vector<size_t> > &faceIndices) 
	{
		PROFILE_SCOPED()
		for (size_t verti = 0; verti < vertices.size(); verti++)
		{
			vector3f sumNorm(0.0f);
			// Calculate and add the normal of every face that contains this vertex
			for (auto fi : faceIndices[verti]) {
				const size_t faceidx = (fi * 3);
				const vector3f v01 = (vertices[indices[faceidx + 0]].pos - vertices[indices[faceidx + 1]].pos).Normalized();
				const vector3f v02 = (vertices[indices[faceidx + 0]].pos - vertices[indices[faceidx + 2]].pos).Normalized();
				sumNorm += v01.Cross(v02);
			}
			vertices[verti].norm = sumNorm.Normalized();
		}
	}
}

static const Sint32 MAX_SUBDIVS = 5;
static const float ICOSX = 0.525731112119133f;
static const float ICOSZ = 0.850650808352039f;

static const vector3f icosahedron_vertices[12] = {
	vector3f(-ICOSX, 0.0, ICOSZ), vector3f(ICOSX, 0.0, ICOSZ), vector3f(-ICOSX, 0.0, -ICOSZ), vector3f(ICOSX, 0.0, -ICOSZ),
	vector3f(0.0, ICOSZ, ICOSX), vector3f(0.0, ICOSZ, -ICOSX), vector3f(0.0, -ICOSZ, ICOSX), vector3f(0.0, -ICOSZ, -ICOSX),
	vector3f(ICOSZ, ICOSX, 0.0), vector3f(-ICOSZ, ICOSX, 0.0), vector3f(ICOSZ, -ICOSX, 0.0), vector3f(-ICOSZ, -ICOSX, 0.0)
};

static const Sint32 icosahedron_faces[20][3] = {
	{0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},
	{8,10,1}, {8,3,10},{5,3,8}, {5,2,3}, {2,7,3},
	{7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6},
	{6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11}
};
#pragma optimize("",off)
Asteroid::Asteroid(Renderer *renderer, RefCountedPtr<Material> mat, Graphics::RenderState *state,
	const TDeformations &deformations, const Sint32 subdivs, const float scale)
{
	PROFILE_SCOPED();
	
	m_material = mat;
	m_renderState = state;

	const Sint32 subdivsLocal = Clamp(subdivs, 0, MAX_SUBDIVS);
	const float scaleLocal = fabs(scale);
	matrix4x4f trans = matrix4x4f::Identity();
	trans.Scale(scaleLocal, scaleLocal, scaleLocal);

	//reserve some data
	VertexArray vts(ATTRIB_POSITION | ATTRIB_NORMAL | ATTRIB_UV0 | ATTRIB_TANGENT, 65536);
	std::vector<Uint16> indices;

	//initial vertices
	Sint32 vi[12];
	for (Sint32 i = 0; i < 12; i++) {
		const vector3f &v = icosahedron_vertices[i];
		vi[i] = AddVertex(vts, trans * v, v);
	}

	//subdivide
	for (Sint32 i = 0; i < 20; i++) {
		Subdivide(vts, indices, trans,
			icosahedron_vertices[icosahedron_faces[i][0]],
			icosahedron_vertices[icosahedron_faces[i][1]],
			icosahedron_vertices[icosahedron_faces[i][2]],
			vi[icosahedron_faces[i][0]],
			vi[icosahedron_faces[i][1]],
			vi[icosahedron_faces[i][2]],
			subdivsLocal);
	}

	std::vector<PosNormTangentUVVert> vertices;
	vertices.resize(vts.GetNumVerts());
	for (Uint32 i = 0; i < vts.GetNumVerts(); i++)
	{
		vertices[i].pos = vts.position[i];
		vertices[i].norm = vts.normal[i];
		vertices[i].tangent = vts.tangent[i];
		vertices[i].uv0 = vts.uv0[i];
	}

	// eliminate duplicate vertices
	{
		std::vector<Uint32> xrefs;
		nv::Weld<PosNormTangentUVVert> weld;
		weld(vertices, xrefs);

		// Remap faces.
		const size_t faceCount = indices.size() / 3;
		for (size_t f = 0; f < faceCount; f++)
		{
			const size_t idx = (f * 3);
			indices[idx + 0] = xrefs.at(indices[idx + 0]);
			indices[idx + 1] = xrefs.at(indices[idx + 1]);
			indices[idx + 2] = xrefs.at(indices[idx + 2]);
		}
	}

	// setup a random number generator to choose the indices
	const size_t num_deformations = deformations.size();
	const size_t s = (num_deformations * 2) + 2;
	std::vector<Uint32> _init;
	_init.reserve(s);
	_init.push_back(num_deformations);
	for (auto i : deformations) {
		_init.push_back(Uint32(i.radius * 100.0f));
		_init.push_back(Uint32(i.offset * 100.0f));
	}
	_init.push_back(UNIVERSE_SEED);
	Random rand(&_init[0], s);


	std::vector < std::vector<size_t> > faceIndices;
	faceIndices.resize(vertices.size());
	// Create vertex to face listing
	const size_t faceCount = indices.size() / 3;
	for (size_t verti = 0; verti < vertices.size(); verti++)
	{
		// Find all of the faces that use this vertex
		for (size_t f = 0; f < faceCount; f++)
		{
			const size_t idx = (f * 3);
			if ((indices[idx + 0] == verti) || (indices[idx + 1] == verti) || (indices[idx + 2] == verti))
			{
				faceIndices[verti].push_back(f);
			}
		}
	}

	// radius, depth, and the number of impacts
	std::set<Uint16> usedIndices;
	const size_t num_indices = indices.size();
	for (auto i : deformations) {
		// get the vertex that forms the basis of our deformation
		Uint16 prospective = indices[(rand.Int32() % num_indices)];
		while (usedIndices.find(prospective) != usedIndices.end()) {
			prospective = indices[(rand.Int32() % num_indices)];
		}
		const Uint16 idx = prospective;
		usedIndices.insert(idx);
		const PosNormTangentUVVert &vert = vertices[idx];

		// find all vertices within the target radius
		const float squareRad = (i.radius * scaleLocal) * (i.radius * scaleLocal);
		std::vector<std::pair<Uint32, float>> idxDistSqr;
		for (Uint32 vidx = 0; vidx < vertices.size(); vidx++) {
			// skip the current vertex
			if (vidx == idx)
				continue;

			// test, store
			const PosNormTangentUVVert &cur = vertices[vidx];
			const float distSqr = (vert.pos - cur.pos).LengthSqr();
			if (distSqr < squareRad) {
				idxDistSqr.push_back(std::make_pair(vidx, distSqr));
			}
		}
		// add the current verteex too of course
		idxDistSqr.push_back(std::make_pair(idx, 0.0f));

		// move the vertices by the offset amount in the direction of the "vert" normal, scaled by distance from centre of the radius
		for (auto vipair : idxDistSqr) {
			const float sclOffset = (i.offset * scaleLocal) * (1.0f - (vipair.second / squareRad));
			const vector3f dirOffset = (vert.norm * sclOffset);
			vertices[vipair.first].pos += dirOffset;
		}

		// regenerate vertex normals after each deformation so the new normals affect future iterations
		CalculatePerVertexNormals(vertices, indices, faceIndices);
	}

	// regenerate vertex normals a final time
	CalculatePerVertexNormals(vertices, indices, faceIndices);

	CalculateTangentArray(vertices, indices);

	// calculate the slope and "height"
	float min = 2.0f;
	float max = 0.0f;
	for (size_t verti = 0; verti < vertices.size(); verti++)
	{
		const vector3f pos = vertices[verti].pos.Normalized();
		const float dot = 1.0f - pos.Dot(vertices[verti].norm.Normalized());
		min = std::min(min, dot);
		max = std::max(max, dot);
		vertices[verti].uv0.x = Clamp(dot, 0.0f, 1.0f);
		vertices[verti].uv0.y = Clamp(pos.Length() - scale, 0.0f, 1.0f);
	}
	assert(min <= max);
	Output("min (%5.2f), max (%5.2f)\n", min, max);

	//Create vtx & index buffers and copy data
	VertexBufferDesc vbd;
	vbd.attrib[0].semantic = ATTRIB_POSITION;
	vbd.attrib[0].format   = ATTRIB_FORMAT_FLOAT3;
	vbd.attrib[1].semantic = ATTRIB_NORMAL;
	vbd.attrib[1].format   = ATTRIB_FORMAT_FLOAT3;
	vbd.attrib[2].semantic = ATTRIB_UV0;
	vbd.attrib[2].format   = ATTRIB_FORMAT_FLOAT2;
	vbd.attrib[3].semantic = ATTRIB_TANGENT;
	vbd.attrib[3].format   = ATTRIB_FORMAT_FLOAT3;
	vbd.numVertices = vts.GetNumVerts();
	vbd.usage = BUFFER_USAGE_STATIC;
	m_vertexBuffer.reset(renderer->CreateVertexBuffer(vbd));
	
	// vertices
	PosNormTangentUVVert* vtxPtr = m_vertexBuffer->Map<PosNormTangentUVVert>(Graphics::BUFFER_MAP_WRITE);
	assert(m_vertexBuffer->GetDesc().stride == sizeof(PosNormTangentUVVert));
	for (size_t i = 0; i<vertices.size(); i++)
	{
		vtxPtr[i].pos = vertices[i].pos;
		vtxPtr[i].norm = vertices[i].norm;
		vtxPtr[i].tangent = vertices[i].tangent;
		vtxPtr[i].uv0 = vertices[i].uv0;
	}
	m_vertexBuffer->Unmap();

	// indices
	m_indexBuffer.reset(renderer->CreateIndexBuffer(indices.size(), BUFFER_USAGE_STATIC));
	Uint16 *idxPtr = m_indexBuffer->Map(Graphics::BUFFER_MAP_WRITE);
	for (auto it : indices) {
		*idxPtr = it;
		idxPtr++;
	}
	m_indexBuffer->Unmap();

#ifdef DEBUG_RENDER_NORMALS
	static float s_lineScale(0.125f);
	m_normLines.reset(new Graphics::Drawables::Lines());
	std::vector<vector3f> lines;
	lines.reserve(vertices.size() * 2);
	for (size_t i = 0; i<vertices.size(); i++)
	{
		lines.push_back(vertices[i].pos);
		lines.push_back(vertices[i].pos + (vertices[i].norm * s_lineScale));
	}
	m_normLines->SetData(vertices.size() * 2, &lines[0], Color::GREEN);

	m_tangentLines.reset(new Graphics::Drawables::Lines());
	lines.clear();
	for (size_t i = 0; i<vertices.size(); i++)
	{
		lines.push_back(vertices[i].pos);
		lines.push_back(vertices[i].pos + (vertices[i].tangent * s_lineScale));
	}
	m_tangentLines->SetData(vertices.size() * 2, &lines[0], Color::RED);

	m_bitangentLines.reset(new Graphics::Drawables::Lines());
	lines.clear();
	for (size_t i = 0; i<vertices.size(); i++)
	{
		lines.push_back(vertices[i].pos);
		lines.push_back(vertices[i].pos + (vertices[i].norm.Cross(vertices[i].tangent) * s_lineScale));
	}
	m_bitangentLines->SetData(vertices.size() * 2, &lines[0], Color::BLUE);
#endif // DEBUG_RENDER_NORMALS
}

void Asteroid::Draw(Renderer *r)
{
	PROFILE_SCOPED()
	r->DrawBufferIndexed(m_vertexBuffer.get(), m_indexBuffer.get(), m_renderState, m_material.Get());
#ifdef DEBUG_RENDER_NORMALS
	m_normLines->Draw(r, m_renderState, Graphics::LINE_SINGLE);
	m_tangentLines->Draw(r, m_renderState, Graphics::LINE_SINGLE);
	m_bitangentLines->Draw(r, m_renderState, Graphics::LINE_SINGLE);
#endif // DEBUG_RENDER_NORMALS
}

Sint32 Asteroid::AddVertex(VertexArray &vts, const vector3f &v, const vector3f &n)
{
	PROFILE_SCOPED()
	vts.position.push_back(v);
	vts.normal.push_back(n);
	//http://www.mvps.org/directx/articles/spheremap.htm
	vts.uv0.push_back(vector2f(asinf(n.x) / M_PI + 0.5f, asinf(n.y) / M_PI + 0.5f));
	vts.tangent.push_back(n);
	return vts.GetNumVerts() - 1;
}

void Asteroid::AddTriangle(std::vector<Uint16> &indices, const Sint32 i1, const Sint32 i2, const Sint32 i3)
{
	PROFILE_SCOPED()
	indices.push_back(i1);
	indices.push_back(i2);
	indices.push_back(i3);
}

void Asteroid::Subdivide(VertexArray &vts, std::vector<Uint16> &indices,
	const matrix4x4f &trans, const vector3f &v1, const vector3f &v2, const vector3f &v3,
	const Sint32 i1, const Sint32 i2, const Sint32 i3, const Sint32 depth)
{
	PROFILE_SCOPED()
	if (depth == 0) {
		AddTriangle(indices, i1, i3, i2);
		return;
	}

	const vector3f v12 = (v1+v2).Normalized();
	const vector3f v23 = (v2+v3).Normalized();
	const vector3f v31 = (v3+v1).Normalized();
	const Sint32 i12 = AddVertex(vts, trans * v12, v12);
	const Sint32 i23 = AddVertex(vts, trans * v23, v23);
	const Sint32 i31 = AddVertex(vts, trans * v31, v31);
	Subdivide(vts, indices, trans, v1, v12, v31, i1, i12, i31, depth-1);
	Subdivide(vts, indices, trans, v2, v23, v12, i2, i23, i12, depth-1);
	Subdivide(vts, indices, trans, v3, v31, v23, i3, i31, i23, depth-1);
	Subdivide(vts, indices, trans, v12, v23, v31, i12, i23, i31, depth-1);
}
//------------------------------------------------------------
