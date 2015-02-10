// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _ASTEROID_H
#define _ASTEROID_H

#include "libs.h"
#include "CollMesh.h"
#include "collider/Geom.h"
#include "graphics/Drawables.h"

//------------------------------------------------------------
//#define DEBUG_RENDER_NORMALS

// Three dimensional sphere (subdivided icosahedron) with normals
// and spherical texture coordinates.
class Asteroid {
public:
	struct TDeform {
		float radius;
		float offset; // depth/height
	};
	typedef std::vector<TDeform> TDeformations;

	// subdivisions must be 0-5
	Asteroid(Graphics::Renderer*, RefCountedPtr<Graphics::Material> material, Graphics::RenderState*, const TDeformations &deformations, const Sint32 subdivisions = 0, const float scale = 1.f);
	virtual void Draw(Graphics::Renderer *r);

	RefCountedPtr<Graphics::Material> GetMaterial() const { return m_material; }

private:
#ifdef DEBUG_RENDER_NORMALS
	std::unique_ptr<Graphics::Drawables::Lines> m_normLines;
	std::unique_ptr<Graphics::Drawables::Lines> m_tangentLines;
	std::unique_ptr<Graphics::Drawables::Lines> m_bitangentLines;
#endif // DEBUG_RENDER_NORMALS

	std::unique_ptr<Graphics::VertexBuffer> m_vertexBuffer;
	std::unique_ptr<Graphics::IndexBuffer> m_indexBuffer;
	RefCountedPtr<Graphics::Material> m_material;
	Graphics::RenderState *m_renderState;

	RefCountedPtr<CollMesh> m_collMesh;
	std::unique_ptr<Geom> m_geom; //static geom

	// Starts the icoshedron generation
	void GenerateInitialMesh(Graphics::VertexArray &vts, std::vector<Uint16> &indices, const matrix4x4f &trans, const Sint32 subdivsLocal);
	// add a new vertex, return the index
	Sint32 AddVertex(Graphics::VertexArray&, const vector3f &v, const vector3f &n);
	// add three vertex indices to form a triangle
	void AddTriangle(std::vector<Uint16>&, const Sint32 i1, const Sint32 i2, const Sint32 i3);
	void Subdivide(Graphics::VertexArray&, std::vector<Uint16>&,
		const matrix4x4f &trans, const vector3f &v1, const vector3f &v2, const vector3f &v3,
		const Sint32 i1, const Sint32 i2, const Sint32 i3, const Sint32 depth);
	// create the render buffers
	void CreateAndPopulateRenderBuffers(const Graphics::VertexArray &vts, std::vector<Uint16> &indices, Graphics::Renderer *r);
};
//------------------------------------------------------------

#endif
