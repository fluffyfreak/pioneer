// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Asteroid.h"
#include "Graphics/Texture.h"
#include "collider/Weld.h"
#include "gameconsts.h"
#include <set>

using namespace Graphics;
using namespace Graphics::Drawables;

//------------------------------------------------------------

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
	VertexArray vts(ATTRIB_POSITION | ATTRIB_NORMAL | ATTRIB_UV0, 65536);
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
	
	struct PosNormUVVert {
		vector3f pos;
		vector3f norm;
		vector2f uv;
		bool operator==(const PosNormUVVert &a) const {
			return (pos == a.pos) && (norm == a.norm) && (uv == a.uv);
		}
	};
	std::vector<PosNormUVVert> vertices;
	vertices.resize(vts.GetNumVerts());
	for (Uint32 i = 0; i < vts.GetNumVerts(); i++)
	{
		vertices[i].pos = vts.position[i];
		vertices[i].norm = vts.normal[i];
		vertices[i].uv = vts.uv0[i];
	}

	// eliminate duplicate vertices
	{
		std::vector<Uint32> xrefs;
		nv::Weld<PosNormUVVert> weld;
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
		const PosNormUVVert &vert = vertices[idx];

		// find all vertices within the target radius
		const float squareRad = (i.radius * scaleLocal) * (i.radius * scaleLocal);
		std::vector<std::pair<Uint32, float>> idxDistSqr;
		for (Uint32 vi = 0; vi < vertices.size(); vi++) {
			// skip the current vertex
			if (vi == idx)
				continue;

			// test, store
			const PosNormUVVert &cur = vertices[vi];
			const float distSqr = (vert.pos - cur.pos).LengthSqr();
			if (distSqr < squareRad) {
				idxDistSqr.push_back(std::make_pair(vi, distSqr));
			}
		}
		// add the current verteex too of course
		idxDistSqr.push_back(std::make_pair(idx, 0.0f));

		// move the vertices by the offset amount in the direction of the "vert" normal, scaled by distance from centre of the radius
		for (auto vi : idxDistSqr) {
			const float sclOffset = (i.offset * scaleLocal) * (1.0f - (vi.second / squareRad));
			const vector3f dirOffset = (vert.norm * sclOffset);
			vertices[vi.first].pos += dirOffset;
		}

		// regenerate vertex normals after each deformation so the new normals affect future iterations
		for (size_t verti = 0; verti < vertices.size(); verti++)
		{
			vector3f sumNorm(0.0f);
			// Calculate and add the normal of every face that contains this vertex
			for (auto fi : faceIndices[verti]) {
				const size_t idx = (fi * 3);
				const vector3f v01 = (vertices[indices[idx + 0]].pos - vertices[indices[idx + 1]].pos).Normalized();
				const vector3f v02 = (vertices[indices[idx + 0]].pos - vertices[indices[idx + 2]].pos).Normalized();
				sumNorm += v01.Cross(v02);
			}
			vertices[verti].norm = sumNorm.Normalized();
		}
	}

	// regenerate vertex normals a final time
	for (size_t verti = 0; verti < vertices.size(); verti++)
	{
		vector3f sumNorm(0.0f);
		// Calculate and add the normal of every face that contains this vertex
		for (auto fi : faceIndices[verti]) {
			const size_t idx = (fi * 3);
			const vector3f v01 = (vertices[indices[idx + 0]].pos - vertices[indices[idx + 1]].pos).Normalized();
			const vector3f v02 = (vertices[indices[idx + 0]].pos - vertices[indices[idx + 2]].pos).Normalized();
			sumNorm += v01.Cross(v02);
		}
		vertices[verti].norm = sumNorm.Normalized();
	}

	//Create vtx & index buffers and copy data
	VertexBufferDesc vbd;
	vbd.attrib[0].semantic = ATTRIB_POSITION;
	vbd.attrib[0].format   = ATTRIB_FORMAT_FLOAT3;
	vbd.attrib[1].semantic = ATTRIB_NORMAL;
	vbd.attrib[1].format   = ATTRIB_FORMAT_FLOAT3;
	vbd.attrib[2].semantic = ATTRIB_UV0;
	vbd.attrib[2].format   = ATTRIB_FORMAT_FLOAT2;
	vbd.numVertices = vts.GetNumVerts();
	vbd.usage = BUFFER_USAGE_STATIC;
	m_vertexBuffer.reset(renderer->CreateVertexBuffer(vbd));
	
	// vertices
	PosNormUVVert* vtxPtr = m_vertexBuffer->Map<PosNormUVVert>(Graphics::BUFFER_MAP_WRITE);
	assert(m_vertexBuffer->GetDesc().stride == sizeof(PosNormUVVert));
	for (size_t i = 0; i<vertices.size(); i++)
	{
		vtxPtr[i].pos = vertices[i].pos;
		vtxPtr[i].norm = vertices[i].norm;
		vtxPtr[i].uv = vertices[i].uv;
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
}

void Asteroid::Draw(Renderer *r)
{
	PROFILE_SCOPED()
	r->DrawBufferIndexed(m_vertexBuffer.get(), m_indexBuffer.get(), m_renderState, m_material.Get());
}

Sint32 Asteroid::AddVertex(VertexArray &vts, const vector3f &v, const vector3f &n)
{
	PROFILE_SCOPED()
	vts.position.push_back(v);
	vts.normal.push_back(n);
	//http://www.mvps.org/directx/articles/spheremap.htm
	vts.uv0.push_back(vector2f(asinf(n.x)/M_PI+0.5f, asinf(n.y)/M_PI+0.5f));
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
