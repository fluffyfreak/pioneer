// Copyright © 2008-2013 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Renderable.h"
#include "Sfx.h"
#include "Projectile.h"
#include "graphics/Drawables.h"
#include "graphics/Renderer.h"
#include "graphics/Material.h"
#include "graphics/VertexArray.h"
#include "scenegraph/Model.h"
#include "scenegraph/Node.h"

using Graphics::Renderer;
using Graphics::VertexArray;
using Graphics::Material;

TriangleGraphic::TriangleGraphic(Renderer *renderer, VertexArray *v, RefCountedPtr<Graphics::Material> m, Graphics::PrimitiveType pt)
: Renderable(renderer)
, m_vertexArray(v)
, m_material(m)
, m_primitiveType(pt)
{
}

VertexArray *TriangleGraphic::GetVertexArray() const
{
	return m_vertexArray.get();
}

void TriangleGraphic::SetVertexArray(VertexArray *v)
{
	m_vertexArray.reset(v);
}

void TriangleGraphic::Draw(Graphics::RenderState *rs)
{
	m_renderer->SetTransform(m_transform);
	m_renderer->DrawTriangles(m_vertexArray.get(), rs, m_material.Get(), m_primitiveType);
}

LaserGraphic::LaserGraphic(Renderer *r)
: Renderable(r)
, m_color(Color::WHITE)
, m_sideIntensity(0.f)
, m_glowIntensity(0.f)
{
}

void LaserGraphic::Draw(Graphics::RenderState *rs)
{
	//Using the shared resources in Projectile
	GetRenderer()->SetTransform(m_transform);
	if (m_sideIntensity > 0.01f) {
		Projectile::GetSideMat()->diffuse = m_color * m_sideIntensity;
		GetRenderer()->DrawTriangles(Projectile::GetSideVerts(), rs, Projectile::GetSideMat());
	}
	if (m_glowIntensity > 0.01f) {
		Projectile::GetGlowMat()->diffuse = m_color * m_glowIntensity;
		GetRenderer()->DrawTriangles(Projectile::GetGlowVerts(), rs, Projectile::GetGlowMat());
	}
}

ModelGraphic::ModelGraphic(SceneGraph::Model *m, unsigned int n)
: Renderable(0)
, m_model(m)
, m_nodeMask(n)
{
}

void ModelGraphic::Draw(Graphics::RenderState *rs)
{
	SceneGraph::RenderData rd;
	rd.nodemask = m_nodeMask;
	m_model->Render(m_transform, &rd);
}