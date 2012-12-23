// Copyright © 2008-2012 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "GeoSphereMaterial.h"
#include "GeoSphere.h"
#include "GeoContext.h"
#include "StringF.h"
#include "graphics/Graphics.h"
#include "graphics/TextureGL.h"
#include "graphics/RendererGL2.h"
#include "graphics/GL2/GL2Debug.h"
#include <sstream>

namespace Graphics {
namespace GL2 {

GeoSphereProgram::GeoSphereProgram(const std::string &filename, const std::string &defines)
{
	m_name = filename;
	m_defines = defines;
	LoadShaders(filename, defines);
	InitUniforms();
	Graphics::GL2::CheckGLError();
}

void GeoSphereProgram::InitUniforms()
{
	Program::InitUniforms();
	atmosColor.Init("atmosColor", m_program);
	geosphereAtmosFogDensity.Init("geosphereAtmosFogDensity", m_program);
	geosphereAtmosInvScaleHeight.Init("geosphereAtmosInvScaleHeight", m_program);
	geosphereAtmosTopRad.Init("geosphereAtmosTopRad", m_program);
	geosphereCenter.Init("geosphereCenter", m_program);
	geosphereScale.Init("geosphereScale", m_program);
	geosphereScaledRadius.Init("geosphereScaledRadius", m_program);
}

Program *GeoSphereSurfaceMaterial::CreateProgram(const MaterialDescriptor &desc)
{
	assert(desc.effect == EFFECT_GEOSPHERE_TERRAIN);
	assert(desc.dirLights < 5);
	std::stringstream ss;
	ss << stringf("#define NUM_LIGHTS %0{u}\n", desc.dirLights);
	if (desc.atmosphere)
		ss << "#define ATMOSPHERE\n";
	return new Graphics::GL2::GeoSphereProgram("geosphere_terrain", ss.str());
}

void GeoSphereSurfaceMaterial::Apply()
{
	SetGSUniforms();
}

void GeoSphereSurfaceMaterial::SetGSUniforms()
{
	GeoSphereProgram *p = static_cast<GeoSphereProgram*>(m_program);
	const SystemBody::AtmosphereParameters ap = *static_cast<SystemBody::AtmosphereParameters*>(this->specialParameter0);

	p->Use();
	p->invLogZfarPlus1.Set(m_renderer->m_invLogZfarPlus1);
	p->emission.Set(this->emissive);
	p->sceneAmbient.Set(m_renderer->GetAmbientColor());
	p->atmosColor.Set(ap.atmosCol);
	p->geosphereAtmosFogDensity.Set(ap.atmosDensity);
	p->geosphereAtmosInvScaleHeight.Set(ap.atmosInvScaleHeight);
	p->geosphereAtmosTopRad.Set(ap.atmosRadius);
	p->geosphereCenter.Set(ap.center);
	p->geosphereScaledRadius.Set(ap.planetRadius / ap.scale);
	p->geosphereScale.Set(ap.scale);
}

Program *GeoSphereSkyMaterial::CreateProgram(const MaterialDescriptor &desc)
{
	assert(desc.effect == EFFECT_GEOSPHERE_SKY);
	assert(desc.dirLights > 0 && desc.dirLights < 5);
	std::stringstream ss;
	ss << stringf("#define NUM_LIGHTS %0{u}\n", desc.dirLights);
	ss << "#define ATMOSPHERE\n";
	return new Graphics::GL2::GeoSphereProgram("geosphere_sky", ss.str());
}

void GeoSphereSkyMaterial::Apply()
{
	SetGSUniforms();
}

GeoSphereHeightmapProgram::GeoSphereHeightmapProgram(const std::string &filename, const std::string &defines)
{
	m_name = filename;
	m_defines = defines;
	LoadShaders(filename, defines);
	InitUniforms();
	Graphics::GL2::CheckGLError();
}

void GeoSphereHeightmapProgram::InitUniforms()
{
	Program::InitUniforms();
	AtmosColor.Init("atmosColor", m_program);
	AtmosFogDensity.Init("geosphereAtmosFogDensity", m_program);
	AtmosInvScaleHeight.Init("geosphereAtmosInvScaleHeight", m_program);
	AtmosTopRad.Init("geosphereAtmosTopRad", m_program);
	Center.Init("geosphereCenter", m_program);
	Scale.Init("geosphereScale", m_program);
	ScaledRadius.Init("geosphereScaledRadius", m_program);

	Radius.Init("radius", m_program);
	V0.Init("v0", m_program);
	V1.Init("v1", m_program);
	V2.Init("v2", m_program);
	V3.Init("v3", m_program);
	FracStep.Init("fracStep", m_program);
}

Program *GeoSphereHeightmapMaterial::CreateProgram(const MaterialDescriptor &desc)
{
	assert(desc.effect == EFFECT_GEOSPHERE_HEIGHTMAPTERRAIN);
	assert(desc.dirLights < 5);
	std::stringstream ss;
	ss << stringf("#define NUM_LIGHTS %0{u}\n", desc.dirLights);
	if (desc.atmosphere)
		ss << "#define ATMOSPHERE\n";
	return new Graphics::GL2::GeoSphereProgram("geosphere_heightmapterrain", ss.str());
}

void GeoSphereHeightmapMaterial::Apply()
{
	//XXX replace with actual material parameter
	glMaterialfv (GL_FRONT, GL_EMISSION, &emissive[0]);

	SetGSUniforms();
}

void GeoSphereHeightmapMaterial::SetGSUniforms()
{
	GeoSphereHeightmapProgram *p = static_cast<GeoSphereHeightmapProgram*>(m_program);
	const SystemBody::AtmosphereParameters ap = *static_cast<SystemBody::AtmosphereParameters*>(this->specialParameter0);

	p->Use();
	p->invLogZfarPlus1.Set(m_renderer->m_invLogZfarPlus1);
	p->sceneAmbient.Set(m_renderer->GetAmbientColor());
	p->AtmosColor.Set(ap.atmosCol);
	p->AtmosFogDensity.Set(ap.atmosDensity);
	p->AtmosInvScaleHeight.Set(ap.atmosInvScaleHeight);
	p->AtmosTopRad.Set(ap.atmosRadius);
	p->Center.Set(ap.center);
	p->ScaledRadius.Set(ap.planetRadius / ap.scale);
	p->Scale.Set(ap.scale);

	const GeoPatchParameters gpp = *static_cast<GeoPatchParameters*>(this->specialParameter1);
	p->Radius.Set(ap.planetRadius);
	p->V0.Set(gpp.mV0);
	p->V1.Set(gpp.mV1);
	p->V2.Set(gpp.mV2);
	p->V3.Set(gpp.mV3);
	p->FracStep.Set(1.0f);
}

}
}
