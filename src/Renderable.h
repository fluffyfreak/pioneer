// Copyright © 2008-2013 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _RENDERABLE_H
#define _RENDERABLE_H

#include "libs.h"
#include "graphics/Types.h"

namespace Graphics {
	class Material;
	class Renderer;
	class RenderState;
	class VertexArray;
}

namespace SceneGraph { class Model; }

// A graphic is a combination of a transform and something renderable
// (vertex array + material, StaticMesh)
// A graphic is not shared, but the renderables may be

class Renderable {
public:
	virtual ~Renderable() { }

	Graphics::Renderer *GetRenderer() const { return m_renderer; }

	void SetTransform(const matrix4x4f &transform) { m_transform = transform; }
	const matrix4x4f &GetTransform() const { return m_transform; }

	virtual void Draw(Graphics::RenderState *rs) = 0;

protected:
	Renderable(Graphics::Renderer *r) : m_renderer(r), m_transform(matrix4x4f::Identity()) {}

	Graphics::Renderer *m_renderer;
	matrix4x4f m_transform;
};

//takes ownership of vertex array
class TriangleGraphic : public Renderable {
public:
	TriangleGraphic(Graphics::Renderer*, Graphics::VertexArray*, RefCountedPtr<Graphics::Material>, Graphics::PrimitiveType = Graphics::TRIANGLES);

	RefCountedPtr<Graphics::Material> &GetMaterial() { return m_material; }
	void SetMaterial(RefCountedPtr<Graphics::Material> m) { m_material = m; }

	Graphics::VertexArray *GetVertexArray() const;
	void SetVertexArray(Graphics::VertexArray*);

	virtual void Draw(Graphics::RenderState *rs) override;

private:
	std::unique_ptr<Graphics::VertexArray> m_vertexArray;
	RefCountedPtr<Graphics::Material> m_material;
	Graphics::PrimitiveType m_primitiveType;
};

class LaserGraphic : public Renderable {
public:
	LaserGraphic(Graphics::Renderer *r);

	void SetColor(const Color &c) { m_color = c; }
	void SetSideIntensity(float i) { m_sideIntensity = i; }
	void SetGlowIntensity(float i) { m_glowIntensity = i; }

	virtual void Draw(Graphics::RenderState *rs) override;

private:
	Color m_color;
	float m_sideIntensity;
	float m_glowIntensity;
};

// XXX silly. Just change Model to be a Renderable!
class ModelGraphic : public Renderable {
public:
	ModelGraphic(SceneGraph::Model *model, unsigned int nodeMask);
	virtual void Draw(Graphics::RenderState *rs) override;
	unsigned int GetNodeMask() const { return m_nodeMask; }
	void SetNodeMask(unsigned int i) { m_nodeMask = i; }

private:
	unsigned int m_nodeMask;
	SceneGraph::Model *m_model;
};

#endif // _RENDERABLE_H
