// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "TriplanarMaterial.h"
#include "graphics/Material.h"
#include "graphics/Graphics.h"
#include "RendererGL.h"
#include "TextureGL.h"
#include <sstream>
#include "StringF.h"
#include "Ship.h"

namespace Graphics {
namespace OGL {

TriplanarProgram::TriplanarProgram(const MaterialDescriptor &desc, int numLights)
{
	numLights = Clamp(numLights, 1, 4);

	//build some defines
	std::stringstream ss;
	if (desc.textures > 0)
		ss << "#define TEXTURE0\n";
	if (desc.vertexColors)
		ss << "#define VERTEXCOLOR\n";
	//using only one light
	if (desc.lighting && numLights > 0)
		ss << stringf("#define NUM_LIGHTS %0{d}\n", numLights);
	else
		ss << "#define NUM_LIGHTS 0\n";
	if (desc.normalMap && desc.lighting)
		ss << "#define MAP_NORMAL\n";
	if (desc.specularMap)
		ss << "#define MAP_SPECULAR\n";
	if (desc.glowMap)
		ss << "#define MAP_EMISSIVE\n";
	if (desc.ambientMap)
		ss << "#define MAP_AMBIENT\n";
	if (desc.quality & HAS_HEAT_GRADIENT)
		ss << "#define HEAT_COLOURING\n";

	m_name = "triplanar";
	m_defines = ss.str();

	LoadShaders(m_name, m_defines);
	InitUniforms();
}

void TriplanarProgram::InitUniforms()
{
	Program::InitUniforms();

	scale.Init("scale", m_program);
}

TriplanarMaterial::TriplanarMaterial()
: m_programs()
, m_curNumLights(0)
{
}

Program *TriplanarMaterial::CreateProgram(const MaterialDescriptor &desc)
{
	m_curNumLights = m_renderer->m_numDirLights;
	return new TriplanarProgram(desc, m_curNumLights);
}

void TriplanarMaterial::SetProgram(Program *p)
{
	m_programs[m_curNumLights] = p;
	m_program = p;
}

void TriplanarMaterial::Apply()
{
	//request a new light variation
	if (m_curNumLights != m_renderer->m_numDirLights) {
		m_curNumLights = m_renderer->m_numDirLights;
		if (m_programs[m_curNumLights] == 0) {
			m_descriptor.dirLights = m_curNumLights; //hax
			m_programs[m_curNumLights] = m_renderer->GetOrCreateProgram(this);
		}
		m_program = m_programs[m_curNumLights];
	}

	OGL::Material::Apply();

	TriplanarProgram *p = static_cast<TriplanarProgram*>(m_program);

	p->scale.Set(1.0f/4.0f);

	p->diffuse.Set(this->diffuse);
	p->emission.Set(this->emissive);
	p->specular.Set(this->specular);
	p->shininess.Set(float(this->shininess));
	p->sceneAmbient.Set(m_renderer->GetAmbientColor());
	
	//Light uniform parameters
	for (Uint32 i = 0; i<m_renderer->GetNumLights(); i++) {
		const Light& Light = m_renderer->GetLight(i);
		p->lights[i].diffuse.Set(Light.GetDiffuse());
		p->lights[i].specular.Set(Light.GetSpecular());
		const vector3f pos = Light.GetPosition();
		p->lights[i].position.Set(pos.x, pos.y, pos.z, (Light.GetType() == Light::LIGHT_DIRECTIONAL ? 0.f : 1.f));
	}
	RendererOGL::CheckErrors();

	p->texture0.Set(this->texture0, 0);
	p->texture1.Set(this->texture1, 1);
	p->texture2.Set(this->texture2, 2);
	p->texture3.Set(this->texture3, 3);
	p->texture4.Set(this->texture4, 4);
	p->texture5.Set(this->texture5, 5);

	p->heatGradient.Set(this->heatGradient, 6);
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

void TriplanarMaterial::Unapply()
{
	// Might not be necessary to unbind textures, but let's not old graphics code (eg, old-UI)
	if (heatGradient) {
		static_cast<TextureGL*>(heatGradient)->Unbind();
		glActiveTexture(GL_TEXTURE5);
	}
	if (texture5) {
		static_cast<TextureGL*>(texture5)->Unbind();
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
