// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Asteroid.h"
#include "graphics/Texture.h"
#include "collider/Weld.h"
#include "gameconsts.h"
#include "StringF.h"
#include "utils.h"
#include "Easing.h"
#include "perlin.h"
#include "scenegraph/LOD.h"
#include "scenegraph/Model.h"
#include "scenegraph/Node.h"
#include "scenegraph/StaticGeometry.h"
#include <set>

using namespace Graphics;
using namespace Graphics::Drawables;

//------------------------------------------------------------

typedef std::pair<Uint32, float> TIdxDist;
typedef std::vector<TIdxDist> TIdxDistSqr;

//------------------------------------------------------------

namespace Spatial
{
	static const size_t kNUM_CHILDREN = 8;
	class Octree {
	public:
		struct OctNode;
		typedef bool(*callback)(OctNode &o, const float sqrRad, const vector3f &pos, TIdxDistSqr &idxDistSqr, const bool bRemove);

		Octree(const float dimensions, const size_t maxDepth) : m_maxDepth(maxDepth) {
			PROFILE_SCOPED();
			m_root.reset(new OctNode(vector3f(0.0f), dimensions, 0, maxDepth));
		}

		void Build(std::vector<vector3f> &verts) {
			PROFILE_SCOPED();
			Clear();
			const size_t numVerts = verts.size();
			for (size_t i = 0; i < numVerts; i++) {
				m_root->Insert(verts[i], i);
			}
		}

		void Insert(const vector3f &pos, const size_t idx) {
			PROFILE_SCOPED();
			m_root->Insert(pos, idx);
		}

		bool Traverse(callback proc, const float sqrRad, const vector3f &pos, TIdxDistSqr &idxDistSqr, const bool bRemove) {
			PROFILE_SCOPED();
			return m_root->Traverse(proc, sqrRad, pos, idxDistSqr, bRemove);
		}

		void Clear() {
			PROFILE_SCOPED();
			m_root->Clear();
		}

		struct OctNode {
			OctNode(const vector3f &centre, const float radius, const size_t depth, const size_t maxDepth)
				: m_centre(centre), m_radius(radius), m_depth(depth), m_maxDepth(maxDepth)
			{
				PROFILE_SCOPED();
				if (m_depth < m_maxDepth) {
					const float half_rad = (radius * 0.5f);
					const size_t nextDepth = (depth + 1);
					for (size_t i = 0; i < kNUM_CHILDREN; i++) {
						m_children[i] = new OctNode(centre + (s_boundsOffsetTable[i] * radius), half_rad, nextDepth, m_maxDepth);
					}
				} else {
					for (size_t i = 0; i < kNUM_CHILDREN; i++) {
						m_children[i] = nullptr;
					}
				}
			}

			~OctNode() {
				PROFILE_SCOPED();
				for (size_t i = 0; i < kNUM_CHILDREN; i++) {
					delete m_children[i];
					m_children[i] = nullptr;
				}
			}

			void Insert(const vector3f &pos, const size_t idx)
			{
				PROFILE_SCOPED();

				if (m_children[0]) {
					// pass onto children
					// calculate the points determinant
					size_t code = 0;
					if (pos.x > m_centre.x) code |= 1;
					if (pos.y > m_centre.y) code |= 2;
					if (pos.z > m_centre.z) code |= 4;
					// pass
					m_children[code]->Insert(pos, idx);
				} else {
					// leaf node
					m_points.push_back(TPoint(pos, idx));
				}
			}

			bool Traverse(callback proc, const float sqrRad, const vector3f &pos, TIdxDistSqr &idxDistSqr, const bool bRemove)
			{
				PROFILE_SCOPED();

				// recursively traverse my children
				if (m_children[0]) {
					for (size_t i = 0; i < kNUM_CHILDREN; i++) {
						m_children[i]->Traverse(proc, sqrRad, pos, idxDistSqr, bRemove);
					}
				} else {
					// Call the callback on this node
					return proc(*this, sqrRad, pos, idxDistSqr, bRemove);
				}

				return true;
			}

			void Clear()
			{
				PROFILE_SCOPED();

				m_points.clear();
				if (m_depth < m_maxDepth) {
					for (size_t i = 0; i < kNUM_CHILDREN; i++) {
						m_children[i]->Clear();
					}
				}
			}

			// once set nothing changes these
			const vector3f m_centre;
			const float m_radius;
			const size_t m_depth;
			const size_t m_maxDepth;

			OctNode *m_children[kNUM_CHILDREN];

			struct TPoint {
				TPoint(const vector3f &pos, const size_t idx) : m_pos(pos), m_idx(idx) {}
				vector3f m_pos;
				size_t m_idx;
			};
			std::vector<TPoint> m_points; // the actual data
		};

	private:
		std::unique_ptr<OctNode> m_root;
		const size_t m_maxDepth;

		static const vector3f s_boundsOffsetTable[kNUM_CHILDREN];
	};

	// define the offset table
	const vector3f Octree::s_boundsOffsetTable[kNUM_CHILDREN] =
	{
		{ -0.5, -0.5, -0.5 },
		{ +0.5, -0.5, -0.5 },
		{ -0.5, +0.5, -0.5 },
		{ +0.5, +0.5, -0.5 },
		{ -0.5, -0.5, +0.5 },
		{ +0.5, -0.5, +0.5 },
		{ -0.5, +0.5, +0.5 },
		{ +0.5, +0.5, +0.5 }
	};
}

//------------------------------------------------------------

namespace AsteroidFunctions
{
	static const vector3d v2(2.0f, 2.0f, 2.0f);
	float fbm(const vector3f &pIN)
	{
		vector3d p(pIN);
		float n = 0.0;
		float scale = 1.0;
		for (int i = 0; i<8; i++) {
			n += noise(p) * scale;
			p *= v2;
			scale *= 0.5;
		}
		return n;
	}

	float Warp(const vector3f &p)
	{
		return fbm(p + vector3f(fbm(p + vector3f(fbm(p)))));
	}

	struct PosNormTangentUVVert {
		vector3f pos;
		vector3f norm;
		vector2f uv0;
		vector3f tangent;
		bool operator==(const PosNormTangentUVVert &a) const {
			return (pos == a.pos) && (norm == a.norm) && (uv0 == a.uv0) && (tangent == a.tangent);
		}
	};

	struct VertFaces {
		VertFaces() : count(0) {}
		void Add(const size_t idx) {
			for (size_t i = 0; i < count; i++) {
				if (faces[i] == idx)
					return;
			}
			faces[count++] = idx;
		}
		size_t faces[6];
		size_t count;
	};
	typedef std::vector<VertFaces> TFaceIndices;

	// Lengyel, Eric. “Computing Tangent Space Basis Vectors for an Arbitrary Mesh”. 
	// Terathon Software 3D Graphics Library, 2001. 
	// http://www.terathon.com/code/tangent.html
	static void CalculateTangentArray(VertexArray &vts, const std::vector<Uint16> &indices)
	{
		PROFILE_SCOPED();
		const size_t vertexCount = vts.GetNumVerts();
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

			const vector3f& v1 = vts.position[i1];
			const vector3f& v2 = vts.position[i2];
			const vector3f& v3 = vts.position[i3];

			const vector2f& w1 = vts.uv0[i1];
			const vector2f& w2 = vts.uv0[i2];
			const vector2f& w3 = vts.uv0[i3];

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
			const vector3f& n = vts.normal[a];
			const vector3f& t = tan1[a];

			// Gram-Schmidt orthogonalize
			vts.tangent[a] = (t - n * n.Dot(t)).Normalized();
		}

		delete[] tan1;
	}

	void CalculatePerVertexNormals(VertexArray &vts, std::vector<Uint16> &indices, TFaceIndices &faceIndices)
	{
		PROFILE_SCOPED();
		for (size_t verti = 0; verti < vts.GetNumVerts(); verti++)
		{
			assert(faceIndices[verti].count > 0);
			vector3f sumNorm(0.0f);
			// Calculate and add the normal of every face that contains this vertex
			for (size_t fi = 0; fi < faceIndices[verti].count; fi++) {
				const size_t faceidx = (faceIndices[verti].faces[fi] * 3);
				const vector3f v01 = (vts.position[indices[faceidx + 0]] - vts.position[indices[faceidx + 1]]).Normalized();
				const vector3f v02 = (vts.position[indices[faceidx + 0]] - vts.position[indices[faceidx + 2]]).Normalized();
				sumNorm += v01.Cross(v02);
			}
			vts.normal[verti] = sumNorm.Normalized();
		}
	}

	void CalculatePerVertexNormalsLimited(VertexArray &vts, std::vector<Uint16> &indices, TFaceIndices &faceIndices, std::vector<std::pair<Uint32, float>> &idxDistSqr)
	{
		PROFILE_SCOPED();
		for (size_t idx = 0; idx < idxDistSqr.size(); idx++)
		{
			const size_t verti = idxDistSqr[idx].first;
			vector3f sumNorm(0.0f);
			// Calculate and add the normal of every face that contains this vertex
			for (size_t fi = 0; fi < faceIndices[verti].count; fi++) {
				const size_t faceidx = (faceIndices[verti].faces[fi] * 3);
				const vector3f v01 = (vts.position[indices[faceidx + 0]] - vts.position[indices[faceidx + 1]]).Normalized();
				const vector3f v02 = (vts.position[indices[faceidx + 0]] - vts.position[indices[faceidx + 2]]).Normalized();
				sumNorm += v01.Cross(v02);
			}
			vts.normal[verti] = sumNorm.Normalized();
		}
	}

	bool ProcessNode(Spatial::Octree::OctNode &o, const float sqrRad, const vector3f &pos, TIdxDistSqr &idxDistSqr, const bool bRemove)
	{
		PROFILE_SCOPED();
		if (((pos - o.m_centre).Length() - (o.m_radius * 2.0f)) < sqrRad) {
			// push back all verts within this node - process them later
			std::vector<Spatial::Octree::OctNode::TPoint>::iterator i = o.m_points.begin();
			while(i != o.m_points.end()) 
			{
				const float distSqr = (pos - i->m_pos).LengthSqr();
				if (distSqr < sqrRad) {
					idxDistSqr.push_back(std::make_pair(i->m_idx, distSqr));
					if (bRemove) {
						o.m_points.erase(i);
					} else {
						++i;
					}
				} else {
					++i;
				}
			}
			return true;
		} 
		return false;
	}

	void findIndexAndDistFrom(Spatial::Octree *pOctree, const float sqrRad, const vector3f &pos, TIdxDistSqr &idxDistSqr, const bool bRemove = true)
	{
		PROFILE_SCOPED();
		assert(pOctree);
		pOctree->Traverse(&ProcessNode, sqrRad, pos, idxDistSqr, bRemove);
	}

	void findIndexAndDistFromBrute(const float sqrRad, const vector3f &pos, const std::vector<vector3f> &verts, TIdxDistSqr &idxDistSqr)
	{
		PROFILE_SCOPED();

		// find all vertices within the target radius
		const size_t numVerts = verts.size();
		for (Uint32 vidx = 0; vidx < numVerts; vidx++) {
			// test, store
			const vector3f &cur = verts[vidx];
			const float distSqr = (pos - cur).LengthSqr();
			if (distSqr < sqrRad) {
				idxDistSqr.push_back(std::make_pair(vidx, distSqr));
			}
		}
	}

	void CreateVertexToFaceList(const size_t numverts, const std::vector<Uint16> &indices, TFaceIndices &faceIndices)
	{
		PROFILE_SCOPED();

		// Create vertex to face listing
		const size_t faceCount = indices.size() / 3;
	
		// Find all of the faces that use this vertex
		faceIndices.resize(numverts);
		for (size_t f = 0; f < faceCount; f++)
		{
			const size_t idx = (f * 3);
			faceIndices[indices[idx + 0]].Add(f);
			faceIndices[indices[idx + 1]].Add(f);
			faceIndices[indices[idx + 2]].Add(f);
		}
#if 0
		TFaceIndices localIndices;
		localIndices.resize(numverts);
		// Create vertex to face listing
		for (size_t verti = 0; verti < numverts; verti++)
		{
			// Find all of the faces that use this vertex
			for (size_t f = 0; f < faceCount; f++)
			{
				const size_t idx = (f * 3);
				if ((indices[idx + 0] == verti) || (indices[idx + 1] == verti) || (indices[idx + 2] == verti))
				{
					localIndices[verti].Add(f);
				}
			}
		}

		// compare
		if( localIndices.size() == faceIndices.size() ) {
			for(size_t ic = 0; ic<localIndices.size(); ic++) {
				if( localIndices[ic].count == faceIndices[ic].count ) {
					for(size_t il=0; il<localIndices[ic].count; il++ ) {
						if(localIndices[ic].faces[il] == faceIndices[ic].faces[il]) {
							Output("ok\n");
						} else {
							assert(false);
							Output("bad\n");
						}
					}
				} else {
					assert(false);
					Output("bad\n");
				}
			}
		} else {
			assert(false);
			Output("bad\n");
		}
#endif
	}

	// convert vts (SOA) to vertex format (AOS)
	void ConvertSOAToAOS(const VertexArray &vtsIn, std::vector<PosNormTangentUVVert> &vOut)
	{
		PROFILE_SCOPED();
		// convert vts (SOA) to vertex format (AOS)
		const size_t numVerts = vtsIn.GetNumVerts();
		vOut.resize(numVerts);
		for (size_t i = 0; i < numVerts; i++)
		{
			vOut[i].pos = vtsIn.position[i];
			vOut[i].norm = vtsIn.normal[i];
			vOut[i].tangent = vtsIn.tangent[i];
			vOut[i].uv0 = vtsIn.uv0[i];
		}
	}

	// convert vts (SOA) to vertex format (AOS)
	void ConvertAOSToSOA(const std::vector<PosNormTangentUVVert> &vIn, VertexArray &vtsOut)
	{
		PROFILE_SCOPED();
		// convert vertex format (AOS) to vts (SOA)
		vtsOut.Clear();
		const size_t numVerts = vIn.size();
		for (size_t i = 0; i < numVerts; i++)
		{
			const PosNormTangentUVVert &v = vIn[i];
			vtsOut.Add(v.pos, v.norm, v.uv0, v.tangent);
		}
	}

	Uint32 EliminateDuplicateVertices(std::vector<PosNormTangentUVVert> &vIn, std::vector<Uint32> &indices)
	{
		PROFILE_SCOPED();
		std::vector<Uint32> xrefs;
		nv::Weld<PosNormTangentUVVert> weld;
		const Uint32 outCount = weld(vIn, xrefs);

		// Remap faces.
		const size_t faceCount = indices.size() / 3;
		for (size_t f = 0; f < faceCount; f++)
		{
			const size_t idx = (f * 3);
			indices[idx + 0] = xrefs.at(indices[idx + 0]);
			indices[idx + 1] = xrefs.at(indices[idx + 1]);
			indices[idx + 2] = xrefs.at(indices[idx + 2]);
		}

		return outCount;
	}

	// Where should this live?
	inline vector3f lerp(const vector3f &start, const vector3f &end, const float n) {
		return (start * (1.0f - n)) + (end * n);
	}
}

using namespace AsteroidFunctions;

static const Sint32 MAX_SUBDIVS = 6;
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


Asteroid::TLOD::TLOD() : m_pOctree(nullptr) 
{
}

Asteroid::TLOD::~TLOD() 
{ 
	if (m_pOctree) { 
		delete m_pOctree; 
	} 
}


void Asteroid::GenerateInitialMesh(VertexArray &vts, std::vector<Uint32> &indices, const matrix4x4f &trans, const Sint32 subdivsLocal)
{
	PROFILE_SCOPED();
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
}

void Asteroid::BuildLods(const size_t subdivsLocal, const matrix4x4f &trans)
{
	PROFILE_SCOPED_DESC("LOD Generation");
	m_lods.resize(subdivsLocal);
	for (Sint32 s = subdivsLocal; s > 0; s--)
	{
		PROFILE_SCOPED_DESC("LOD");
		TLOD *pLOD = new TLOD;
		pLOD->m_pva.reset(new VertexArray(ATTRIB_POSITION | ATTRIB_NORMAL | ATTRIB_UV0 | ATTRIB_TANGENT, 131072));
		m_lods[s - 1] = pLOD;

		// 32-bit indices during generation, before the number of vertices is reduced via welding	
		std::vector<Uint32> indices32;
		indices32.reserve(250000);
		GenerateInitialMesh(*(pLOD->m_pva), indices32, trans, s);

		{
			PROFILE_SCOPED_DESC("LOD Weld");
			// convert vts (SOA) to vertex format (AOS)
			std::vector<PosNormTangentUVVert> vertices;
			ConvertSOAToAOS(*(pLOD->m_pva), vertices);
			// eliminate duplicate vertices
			const Uint32 outCount = EliminateDuplicateVertices(vertices, indices32);
			const size_t numIndices = indices32.size();
			pLOD->m_indices.resize(numIndices);
			for (size_t i32 = 0; i32 < numIndices; i32++) {
				pLOD->m_indices[i32] = indices32[i32];
			}
			// convert back to vts (SOA) from vertex format (AOS)
			ConvertAOSToSOA(vertices, *(pLOD->m_pva));
		}

		if (nullptr == pLOD->m_pOctree) {
			const size_t count = pLOD->m_pva->GetNumVerts();
			const size_t maxDepth = (count < 20000) ? 1 : (count < 60000) ? 2 : 3;
			pLOD->m_pOctree = new Spatial::Octree(m_scale * 4.0f, maxDepth);
		}
		assert(pLOD->m_pOctree);
		pLOD->m_pOctree->Build(pLOD->m_pva->position);

		if (s != subdivsLocal && m_lods[s])
		{
			PROFILE_SCOPED_DESC("LOD Calculate weights");

			// Sorty-McSorty-pants
			struct DistSortStruct{
				static bool DistSortPredicate(const TIdxDist &a, const TIdxDist &b)
				{
					return a.second < b.second;
				}
			};

			const Graphics::VertexArray *va = pLOD->m_pva.get();
			const size_t numVerts = va->GetNumVerts();
			pLOD->m_mappedIdx.reserve(numVerts);
			for (size_t n = 0; n < numVerts; n++) {
				// get position of current vertex
				const vector3f &pos = va->position[n];

				// find all of the closest indices - in the next LOD up (more detail) than us
				TIdxDistSqr idxDistSqr; idxDistSqr.reserve(256);
#if 1
				const bool bRemove = false;
				findIndexAndDistFrom(m_lods[s]->m_pOctree, 0.05f, pos, idxDistSqr, bRemove);
#else
				findIndexAndDistFromBrute(0.05f, pos, m_lods[s]->m_pva->position, idxDistSqr);
#endif
				std::sort(idxDistSqr.begin(), idxDistSqr.end(), DistSortStruct::DistSortPredicate);

				// calculate closest weights
				const Uint16 mappedIdx = idxDistSqr[0].first;
				assert(idxDistSqr[0].second < 0.0001f);
				pLOD->m_mappedIdx.push_back(mappedIdx);
			}
		}
	}
}

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
	m_scale = scaleLocal;

	// build each LOD mesh, weld it to reduce vertices
	BuildLods(subdivsLocal, trans);

	// reserve some data
	VertexArray &vts = *(m_lods[subdivsLocal-1]->m_pva);
	std::vector<Uint16> &indices(m_lods[subdivsLocal - 1]->m_indices); // we only allow 16bit indices ingame so this will hold the final indices

	// setup a random number generator to choose the indices
	const size_t num_deformations = deformations.size();
	const size_t s = (num_deformations * 2) + 2;
	std::vector<Uint32> _init;
	{
		PROFILE_SCOPED_DESC("_init Random");
		_init.reserve(s);
		_init.push_back(num_deformations);
		for (auto i : deformations) {
			_init.push_back(Uint32(i.radius * 100.0f));
			_init.push_back(Uint32(i.offset * 100.0f));
		}
		_init.push_back(UNIVERSE_SEED);
	}
	Random rand(&_init[0], s);

	// This saves a lot of time later on by building the mapping now
	TFaceIndices faceIndices;
	CreateVertexToFaceList(vts.GetNumVerts(), indices, faceIndices);

	Spatial::Octree *pOct = m_lods[subdivsLocal - 1]->m_pOctree;
	assert(pOct);
	pOct->Build(vts.position);

	// radius, depth, and the number of impacts
	const size_t num_indices = indices.size();
	std::set<Uint16> unusedIndices;
	for (auto i : indices) {
		unusedIndices.insert(i);
	}
	float minScl = scaleLocal * 2.0f;
	float maxScl = 0.0f;
	for (auto i : deformations) {
		PROFILE_SCOPED_DESC("Deformation");

		// get the vertex that forms the basis of our deformation
		std::set<Uint16>::iterator it(unusedIndices.begin());
		std::advance(it, (rand.Int32() % unusedIndices.size()));
		const Uint16 idx = *it;
		unusedIndices.erase(it);

		const vector3f &vert = vts.position[idx];
		const vector3f &norm = vts.normal[idx];

		// find all vertices within the target radius
		const float squareRad = (i.radius * scaleLocal) * (i.radius * scaleLocal);
		TIdxDistSqr idxDistSqr; idxDistSqr.reserve(256);
#if 1
		findIndexAndDistFrom(pOct, squareRad, vert, idxDistSqr);
#else
		findIndexAndDistFromBrute(squareRad, vert, vts.position, idxDistSqr);
#endif

		// move the vertices by the offset amount in the direction of the "vert" normal, scaled by distance from centre of the radius
		for (auto vipair : idxDistSqr) {
			PROFILE_SCOPED_DESC("Modify");
			const float sclOffset = (i.offset * scaleLocal) * (1.0f - (vipair.second / squareRad));
			// blend between the centre -> vertex normal and the surface normal to reduce chance of interpenetrating geometry
			const vector3f dirOffset = (norm * sclOffset);
			const vector3f pos = (vts.position[vipair.first].Normalized() * sclOffset);;
			vts.position[vipair.first] += lerp(pos, dirOffset, 0.75f);
			
			const float scl = vts.position[vipair.first].Length() - scaleLocal;
			minScl = std::min(minScl, scl);
			maxScl = std::max(maxScl, scl);

			// update the octree
			pOct->Insert(vts.position[vipair.first], vipair.first);
		}

		// regenerate vertex normals after each deformation so the new normals affect future iterations
		CalculatePerVertexNormalsLimited(vts, indices, faceIndices, idxDistSqr);
	}

	// Deform the surface using some warping noise
	float offsetScaling = (scaleLocal * 0.05f); // TODO: this shoud be dependent on the size of the asteroid being generated?
	const size_t numVerts = vts.GetNumVerts();
	for (size_t verti = 0; verti < numVerts; verti++)
	{
		PROFILE_SCOPED_DESC("Warping");
		const vector3f pos = vts.position[verti].Normalized();
		const vector3f nrm = vts.normal[verti].Normalized();
		const float dot = Clamp(1.0f - pos.Dot(nrm), 0.0f, 1.0f);
		// Scale the deformation according the slope flat surafces should be smooth, vertical get all the noise.
		const float offset = (dot * Clamp(Warp(pos), 0.0f, 1.0f)) * offsetScaling;
		vts.position[verti] = vts.position[verti] + (nrm * offset);
	}

	// regenerate all of the vertex normals
	CalculatePerVertexNormals(vts, indices, faceIndices);

	CalculateTangentArray(vts, indices);

	// calculate the final slope and "height"
	float minDot = 2.0f;
	float maxDot = 0.0f;
	for (size_t verti = 0; verti < numVerts; verti++)
	{
		PROFILE_SCOPED_DESC("UVSlopeHeight");
		const vector3f pos = vts.position[verti].Normalized();
		const vector3f nrm = vts.normal[verti].Normalized();
		const float dot = pos.Dot(nrm);
		//assert(dot >= -1.0f && dot <= 1.0f); - commented out because the fast normalise sometimes produces values greater than 1.0f
		const float scl = (((pos.Length() - scale) - minScl) / (minScl - maxScl));
		minDot = std::min(minDot, dot);
		maxDot = std::max(maxDot, dot);
		vts.uv0[verti].x = 1.0f - Clamp(dot, 0.0f, 1.0f);
		vts.uv0[verti].y = Clamp(scl, 0.0f, 1.0f);
	}
	Output("min dot (%5.2f), max dot (%5.2f)\n", minDot, maxDot);
	Output("min scl (%5.2f), max scl (%5.2f)\n", minScl, maxScl);

	// TODO: decorate

	// TODO: generate/copy data for LODs
	CopyLODsData(subdivsLocal);

	// TODO: create collision mesh/GoemTree etc
	BuildModels(renderer);
}

void Asteroid::CopyLODsData(const size_t subdivs)
{
	PROFILE_SCOPED();
	for (Sint32 s = subdivs - 2; s >= 0; s--)
	{
		PROFILE_SCOPED_DESC("Mesh Update");
		if (s != subdivs && m_lods[s])
		{
			TLOD *pLODCurr = m_lods[s];
			TLOD *pLODNext = m_lods[s + 1];
			const size_t numMap = pLODCurr->m_mappedIdx.size();
			for (size_t mi = 0; mi < numMap; mi++)
			{
				const size_t mIdx = pLODCurr->m_mappedIdx[mi];
				pLODCurr->m_pva->position[mi] = pLODNext->m_pva->position[mIdx];
				pLODCurr->m_pva->normal[mi] = pLODNext->m_pva->normal[mIdx];
				pLODCurr->m_pva->uv0[mi] = pLODNext->m_pva->uv0[mIdx];
				pLODCurr->m_pva->tangent[mi] = pLODNext->m_pva->tangent[mIdx];
			}
		}
	}

}

void Asteroid::BuildModels(Graphics::Renderer *renderer)
{
	PROFILE_SCOPED();
	using namespace SceneGraph;
	m_model.reset( new Model(renderer, "Asteroid") );

	LOD *lodNode = new LOD(renderer);
	m_model->GetRoot()->AddChild(lodNode);

	// material
	m_model->AddMaterial("AsteroidMat", m_material);

	size_t pixelLod = 50;

	for (std::vector<TLOD*>::const_iterator lod = m_lods.begin(); lod != m_lods.end(); ++lod)
	{
		CreateAndPopulateRenderBuffers(*lod, renderer);

		// multiple lods might use the same mesh
		std::vector<RefCountedPtr<StaticGeometry> > geoms;
		RefCountedPtr<StaticGeometry> geom(new StaticGeometry(renderer));
		geom->SetName(stringf("sgMesh%0{u}", 0));
		geom->AddMesh((*lod)->m_vb, (*lod)->m_ib, m_material);

		Graphics::RenderStateDesc rsd;
		geom->SetRenderState(renderer->CreateRenderState(rsd));
		geoms.push_back(geom);

		if (lodNode) {
			lodNode->AddLevel(pixelLod, geom.Get());
			pixelLod += 50;
		}
	}

	// Load collision meshes
	// They are added at the top level of the model root as CollisionGeometry nodes
	/*for (std::vector<std::string>::const_iterator it = def.collisionDefs.begin();
	it != def.collisionDefs.end(); ++it)
	{
	try {
	LoadCollision(*it);
	}
	catch (LoadingError &err) {
	throw (LoadingError(stringf("%0:\n%1", *it, err.what())));
	}
	}

	// Run CollisionVisitor to create the initial CM and its GeomTree.
	// If no collision mesh is defined, a simple bounding box will be generated
	Output("CreateCollisionMesh for : (%s)\n", m_model->m_name.c_str());
	m_model->CreateCollisionMesh();*/

	// Do an initial animation update to get all the animation transforms correct
	m_model->UpdateAnimations();
}

void Asteroid::CreateAndPopulateRenderBuffers(TLOD *pLOD, Graphics::Renderer *r)
{
	PROFILE_SCOPED();

	RefCountedPtr<Graphics::VertexBuffer> &vb = pLOD->m_vb;
	RefCountedPtr<Graphics::IndexBuffer> &ib = pLOD->m_ib;
	const Graphics::VertexArray &vts = *(pLOD->m_pva);
	const size_t numVerts = vts.GetNumVerts();

	//Create vtx & index buffers and copy data
	VertexBufferDesc vbd;
	vbd.attrib[0].semantic = ATTRIB_POSITION;
	vbd.attrib[0].format = ATTRIB_FORMAT_FLOAT3;
	vbd.attrib[1].semantic = ATTRIB_NORMAL;
	vbd.attrib[1].format = ATTRIB_FORMAT_FLOAT3;
	vbd.attrib[2].semantic = ATTRIB_UV0;
	vbd.attrib[2].format = ATTRIB_FORMAT_FLOAT2;
	vbd.attrib[3].semantic = ATTRIB_TANGENT;
	vbd.attrib[3].format = ATTRIB_FORMAT_FLOAT3;
	vbd.numVertices = numVerts;
	vbd.usage = BUFFER_USAGE_STATIC;
	vb.Reset(r->CreateVertexBuffer(vbd));

	// vertices
	PosNormTangentUVVert* vtxPtr = vb->Map<PosNormTangentUVVert>(Graphics::BUFFER_MAP_WRITE);
	assert(vb->GetDesc().stride == sizeof(PosNormTangentUVVert));
	for (size_t i = 0; i<numVerts; i++)
	{
		vtxPtr[i].pos = vts.position[i];
		vtxPtr[i].norm = vts.normal[i];
		vtxPtr[i].tangent = vts.tangent[i];
		vtxPtr[i].uv0 = vts.uv0[i];
	}
	vb->Unmap();

	// indices
	ib.Reset(r->CreateIndexBuffer(pLOD->m_indices.size(), BUFFER_USAGE_STATIC));
	Uint16 *idxPtr = ib->Map(Graphics::BUFFER_MAP_WRITE);
	for (auto it : pLOD->m_indices) {
		*idxPtr = it;
		idxPtr++;
	}
	ib->Unmap();

#ifdef DEBUG_RENDER_NORMALS
	static float s_lineScale(0.125f);
	pLOD->m_normLines.reset(new Graphics::Drawables::Lines());
	std::vector<vector3f> lines;
	lines.reserve(numVerts * 2);
	for (size_t i = 0; i<vts.GetNumVerts(); i++)
	{
		lines.push_back(vts.position[i]);
		lines.push_back(vts.position[i] + (vts.normal[i] * s_lineScale));
	}
	pLOD->m_normLines->SetData(numVerts * 2, &lines[0], Color::GREEN);

	pLOD->m_tangentLines.reset(new Graphics::Drawables::Lines());
	lines.clear();
	for (size_t i = 0; i<numVerts; i++)
	{
		lines.push_back(vts.position[i]);
		lines.push_back(vts.position[i] + (vts.tangent[i] * s_lineScale));
	}
	pLOD->m_tangentLines->SetData(numVerts * 2, &lines[0], Color::RED);

	pLOD->m_bitangentLines.reset(new Graphics::Drawables::Lines());
	lines.clear();
	for (size_t i = 0; i<numVerts; i++)
	{
		lines.push_back(vts.position[i]);
		lines.push_back(vts.position[i] + (vts.normal[i].Cross(vts.tangent[i]) * s_lineScale));
	}
	pLOD->m_bitangentLines->SetData(numVerts * 2, &lines[0], Color::BLUE);
#endif // DEBUG_RENDER_NORMALS
}

void Asteroid::Draw(Renderer *r, const matrix4x4f &trans)
{
	PROFILE_SCOPED();
	static const size_t idx = m_lods.size() - 1;
	const TLOD *pLOD = m_lods[idx];
	if (pLOD) {
		m_model->Render(trans);
		//r->DrawBufferIndexed(pLOD->m_vb.Get(), pLOD->m_ib.Get(), m_renderState, m_material.Get());
#ifdef DEBUG_RENDER_NORMALS
		pLOD->m_normLines->Draw(r, m_renderState, Graphics::LINE_SINGLE);
		pLOD->m_tangentLines->Draw(r, m_renderState, Graphics::LINE_SINGLE);
		pLOD->m_bitangentLines->Draw(r, m_renderState, Graphics::LINE_SINGLE);
#endif // DEBUG_RENDER_NORMALS
	}
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

void Asteroid::AddTriangle(std::vector<Uint32> &indices, const Sint32 i1, const Sint32 i2, const Sint32 i3)
{
	PROFILE_SCOPED()
	indices.push_back(i1);
	indices.push_back(i2);
	indices.push_back(i3);
}

void Asteroid::Subdivide(VertexArray &vts, std::vector<Uint32> &indices,
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
