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
#include "galaxy/StarSystem.h"

namespace Graphics {
namespace GL2 {

GenGasGiantColourProgram::GenGasGiantColourProgram(const MaterialDescriptor &desc)
{
	//build some defines
	std::stringstream ss;
	if (desc.textures > 0)
		ss << "#define TEXTURE0\n";

	switch( desc.effect )
	{
	default:
	case EFFECT_GEN_JUPITER_GASSPHERE_TEXTURE:
		ss << "#define GEN_JUPITER_ESQUE\n";
		break;
	case EFFECT_GEN_NEPTUNE_GASSPHERE_TEXTURE:
		ss << "#define GEN_NEPTUNE_ESQUE\n";
		break;
	case EFFECT_GEN_SATURN_GASSPHERE_TEXTURE:
		ss << "#define GEN_SATURN_ESQUE\n";
		break;
	case EFFECT_GEN_URANUS_GASSPHERE_TEXTURE:
		ss << "#define GEN_URANUS_ESQUE\n";
		break;
	}
	
	// No lights
	ss << "#define NUM_LIGHTS 0\n";

	m_name = "GenGasGiantColour";
	m_defines = ss.str();
	m_libs |= NOISE;	// we need the noise library for generating the colours in the shader

	LoadShaders(m_name, m_defines);
	InitUniforms();
}

void GenGasGiantColourProgram::InitUniforms()
{
	v0.Init("v0", m_program);
	v1.Init("v1", m_program);
	v2.Init("v2", m_program);
	v3.Init("v3", m_program);
	fracStep.Init("fracStep", m_program);

	// XXX omg hacking galore
	time.Init("time", m_program);
	
	octaves.Init("octaves", m_program);
	lacunarity.Init("lacunarity", m_program);
	frequency.Init("frequency", m_program);

	ggdarkColor.Init("ggdarkColor", m_program);
	gglightColor.Init("gglightColor", m_program);
	entropy.Init("entropy", m_program);
	planetEarthRadii.Init("planetEarthRadii", m_program);
	// XXX omg hacking galore
}

Program *GenGasGiantColourMaterial::CreateProgram(const MaterialDescriptor &desc)
{
	return new GenGasGiantColourProgram(desc);
}

void GenGasGiantColourMaterial::Apply()
{
	GenGasGiantColourProgram *p = static_cast<GenGasGiantColourProgram*>(m_program);
	p->Use();

	const Graphics::GenGasGiantColourMaterialParameters params = *static_cast<Graphics::GenGasGiantColourMaterialParameters*>(this->specialParameter0);
	assert(params.v);
	p->v0.Set( params.v[0] );
	p->v1.Set( params.v[1] );
	p->v2.Set( params.v[2] );
	p->v3.Set( params.v[3] );
	p->fracStep.Set( params.fracStep );

	// XXX omg hacking galore
	p->time.Set(params.time);

	assert(params.pTerrain);
	int octaves[10];
	float lacunarity[10];
	float frequency[10];
	for(Uint32 i=0; i<10; i++) {
		octaves[i]		= Clamp(params.pTerrain->GetFracDef(i).octaves, 1, 3);
		lacunarity[i]	= params.pTerrain->GetFracDef(i).lacunarity;
		frequency[i]	= params.pTerrain->GetFracDef(i).frequency;
	}
	p->octaves.Set(octaves, 10);
	p->lacunarity.Set(lacunarity, 10);
	p->frequency.Set(frequency, 10);

	vector3f darkColor[8];
	vector3f lightColor[8];
	for(Uint32 i=0; i<8; i++) {
		const vector3d &dc = params.pTerrain->GetColorValue(Terrain::eGGDARKCOLOR, i);
		const vector3d &lc = params.pTerrain->GetColorValue(Terrain::eGGLIGHTCOLOR, i);
		darkColor[i] = vector3f(dc.x, dc.y, dc.z);
		lightColor[i] = vector3f(lc.x, lc.y, lc.z);
	}
		
	p->ggdarkColor.Set(&darkColor[0], 8);
	p->gglightColor.Set(&lightColor[0], 8);
	p->entropy.Set(float(params.pTerrain->GetEntropy(0)));
	p->planetEarthRadii.Set(float(params.planetRadius / EARTH_RADIUS));
	// XXX omg hacking galore

	p->diffuse.Set(this->diffuse);

	if( this->texture0 )
	{
		p->texture0.Set(this->texture0, 0);
		p->texture1.Set(this->texture1, 1);
		p->texture2.Set(this->texture2, 2);
		p->texture3.Set(this->texture3, 3);
		p->texture4.Set(this->texture4, 4);
	}
}

void GenGasGiantColourMaterial::Unapply()
{
	// Might not be necessary to unbind textures, but let's not old graphics code (eg, old-UI)
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
