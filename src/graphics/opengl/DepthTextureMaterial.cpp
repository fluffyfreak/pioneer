// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "DepthTextureMaterial.h"
#include "graphics/Material.h"
#include "graphics/Graphics.h"
#include "TextureGL.h"
#include "RendererGL.h"
#include <sstream>
#include "StringF.h"

namespace Graphics {
namespace OGL {

// rendering OF a depth map material
DepthTextureProgram::DepthTextureProgram(const MaterialDescriptor &desc)
{
	m_name = "DepthTexture";
	RendererOGL::CheckErrors();

	LoadShaders(m_name, m_defines);
	InitUniforms();
}

Program *DepthTextureMaterial::CreateProgram(const MaterialDescriptor &desc)
{
	return new DepthTextureProgram(desc);
}

void DepthTextureMaterial::Apply()
{
	OGL::Material::Apply();

	DepthTextureProgram *p = static_cast<DepthTextureProgram*>(m_program);

	p->texture0.Set(this->texture0, 0);
	RendererOGL::CheckErrors();
}

void DepthTextureMaterial::Unapply()
{
	// Might not be necessary to unbind textures, but let's not old graphics code (eg, old-UI)
	if (texture0) {
		static_cast<TextureGL*>(texture0)->Unbind();
	}
}


// render-to-depth map material
DepthRenderProgram::DepthRenderProgram(const MaterialDescriptor &desc)
{
	m_name = "DepthRender";
	RendererOGL::CheckErrors();

	LoadShaders(m_name, m_defines);
	InitUniforms();
}

Program *DepthRenderMaterial::CreateProgram(const MaterialDescriptor &desc)
{
	return new DepthRenderProgram(desc);
}

void DepthRenderMaterial::Apply()
{
	OGL::Material::Apply();

	RendererOGL::CheckErrors();
}

}
}
