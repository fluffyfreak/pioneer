// Copyright © 2008-2014 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "GenGasGiantColourMaterial.h"
#include "graphics/Material.h"
#include "graphics/TextureGL.h"
#include "graphics/Graphics.h"
#include "graphics/RendererGL2.h"
#include <sstream>
#include "StringF.h"
#include "Ship.h"

namespace Graphics {
namespace GL2 {

GenGasGiantColourProgram::GenGasGiantColourProgram(const MaterialDescriptor &desc)
{
	//build some defines
	std::stringstream ss;
	if (desc.textures > 0)
		ss << "#define TEXTURE0\n";
	if (desc.vertexColors)
		ss << "#define VERTEXCOLOR\n";
	if (desc.alphaTest)
		ss << "#define ALPHA_TEST\n";
	
	// No lights
	ss << "#define NUM_LIGHTS 0\n";

	if (desc.specularMap)
		ss << "#define MAP_SPECULAR\n";
	if (desc.glowMap)
		ss << "#define MAP_EMISSIVE\n";
	if (desc.usePatterns)
		ss << "#define MAP_COLOR\n";

	m_name = "GenGasGiantColour";
	m_defines = ss.str();

	LoadShaders(m_name, m_defines);
	InitUniforms();
}

Program *GenGasGiantColourMaterial::CreateProgram(const MaterialDescriptor &desc)
{
	return new GenGasGiantColourProgram(desc);
}

void GenGasGiantColourMaterial::Apply()
{
	GenGasGiantColourProgram *p = static_cast<GenGasGiantColourProgram*>(m_program);
	p->Use();
	p->invLogZfarPlus1.Set(m_renderer->m_invLogZfarPlus1);

	p->diffuse.Set(this->diffuse);

	p->texture0.Set(this->texture0, 0);
	p->texture1.Set(this->texture1, 1);
	p->texture2.Set(this->texture2, 2);
	p->texture3.Set(this->texture3, 3);
	p->texture4.Set(this->texture4, 4);

	p->heatGradient.Set(this->heatGradient, 5);
	if(nullptr!=specialParameter0) {
		HeatGradientParameters_t *pMGP = static_cast<HeatGradientParameters_t*>(specialParameter0);
		p->heatingMatrix.Set(pMGP->heatingMatrix);
		p->heatingNormal.Set(pMGP->heatingNormal);
		p->heatingAmount.Set(pMGP->heatingAmount);
	} else {
		p->heatingMatrix.Set(matrix3x3f::Identity());
		p->heatingNormal.Set(vector3f(0.0f, -1.0f, 0.0f));
		p->heatingAmount.Set(0.0f);
	}
}

void GenGasGiantColourMaterial::Unapply()
{
	// Might not be necessary to unbind textures, but let's not old graphics code (eg, old-UI)
	if (heatGradient) {
		static_cast<TextureGL*>(heatGradient)->Unbind();
		glActiveTexture(GL_TEXTURE4);
	}
	if (texture4) {
		static_cast<TextureGL*>(texture4)->Unbind();
		glActiveTexture(GL_TEXTURE3);
	}
	if (texture3) {
		static_cast<TextureGL*>(texture3)->Unbind();
		glActiveTexture(GL_TEXTURE2);
	}
	if (texture2) {
		static_cast<TextureGL*>(texture2)->Unbind();
		glActiveTexture(GL_TEXTURE1);
	}
	if (texture1) {
		static_cast<TextureGL*>(texture1)->Unbind();
		glActiveTexture(GL_TEXTURE0);
	}
	if (texture0) {
		static_cast<TextureGL*>(texture0)->Unbind();
	}
}

}
}
