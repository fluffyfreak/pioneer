// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _ASTEROID_H
#define _ASTEROID_H

#include "graphics/Drawables.h"
#include "libs.h"

//------------------------------------------------------------
//#define DEBUG_RENDER_NORMALS

namespace Spatial {
	class Octree;
};
namespace SceneGraph {
	class Model;
};

// Three dimensional sphere (subdivided icosahedron) with normals
// and spherical texture coordinates.
class Asteroid {
public:
	struct TDeform {
		float radius;
		float offset; // depth/height
	};
	typedef std::vector<TDeform> TDeformations;

	struct TLOD {
		TLOD();
		~TLOD();
		std::unique_ptr<Graphics::VertexArray> m_pva;
		std::vector<Uint16> m_indices;
		std::vector<Uint16> m_mappedIdx;
		Spatial::Octree *m_pOctree;
		// renderable bits
		RefCountedPtr<Graphics::VertexBuffer> m_vb;
		RefCountedPtr<Graphics::IndexBuffer> m_ib;
#ifdef DEBUG_RENDER_NORMALS
		std::unique_ptr<Graphics::Drawables::Lines> m_normLines;
		std::unique_ptr<Graphics::Drawables::Lines> m_tangentLines;
		std::unique_ptr<Graphics::Drawables::Lines> m_bitangentLines;
#endif // DEBUG_RENDER_NORMALS
	};

	// subdivisions must be 0-5
	Asteroid(Graphics::Renderer *, RefCountedPtr<Graphics::Material> material, Graphics::RenderState *, const TDeformations &deformations, const Sint32 subdivisions = 0, const float scale = 1.f);
	virtual void Draw(Graphics::Renderer *r, const matrix4x4f &trans);

	RefCountedPtr<Graphics::Material> GetMaterial() const { return m_material; }

private:
	std::unique_ptr<SceneGraph::Model> m_model;
	RefCountedPtr<Graphics::Material> m_material;
	Graphics::RenderState *m_renderState;
	std::vector<TLOD *> m_lods;
	float m_scale;

	// high level lod builder
	void BuildLods(const size_t subdivsLocal, const matrix4x4f &trans);
	// Starts the icoshedron generation
	void GenerateInitialMesh(Graphics::VertexArray &vts, std::vector<Uint32> &indices, const matrix4x4f &trans, const Sint32 subdivsLocal);
	// add a new vertex, return the index
	Sint32 AddVertex(Graphics::VertexArray &, const vector3f &v, const vector3f &n);
	// add three vertex indices to form a triangle
	void AddTriangle(std::vector<Uint32> &, const Sint32 i1, const Sint32 i2, const Sint32 i3);
	void Subdivide(Graphics::VertexArray &, std::vector<Uint32> &,
		const matrix4x4f &trans, const vector3f &v1, const vector3f &v2, const vector3f &v3,
		const Sint32 i1, const Sint32 i2, const Sint32 i3, const Sint32 depth);

	// create the render buffers
	void CreateAndPopulateRenderBuffers(TLOD *pLOD, Graphics::Renderer *r);

	// copy data for LODs
	void CopyLODsData(const size_t subdivs);

	// Create the SceneGraph models
	void BuildModels(Graphics::Renderer *r);
};
//------------------------------------------------------------

#endif
