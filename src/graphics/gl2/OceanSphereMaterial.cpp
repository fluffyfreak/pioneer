#include "OceanSphereMaterial.h"
#include "OceanSphere.h"
#include "StringF.h"
#include "graphics/Graphics.h"
#include "graphics/RendererGL2.h"
#include <sstream>

namespace Graphics {
namespace GL2 {

OceanSphereProgram::OceanSphereProgram(const std::string &filename, const std::string &defines)
{
	m_name = filename;
	m_defines = defines;
	LoadShaders(filename, defines);
	InitUniforms();
}

void OceanSphereProgram::InitUniforms()
{
	Program::InitUniforms();
	atmosColor.Init("atmosColor", m_program);
	oceansphereAtmosFogDensity.Init("oceansphereAtmosFogDensity", m_program);
	oceansphereAtmosInvScaleHeight.Init("oceansphereAtmosInvScaleHeight", m_program);
	oceansphereAtmosTopRad.Init("oceansphereAtmosTopRad", m_program);
	oceansphereCenter.Init("oceansphereCenter", m_program);
	oceansphereScale.Init("oceansphereScale", m_program);
	oceansphereScaledRadius.Init("oceansphereScaledRadius", m_program);
}

Program *OceanSphereSurfaceMaterial::CreateProgram(const MaterialDescriptor &desc)
{
	assert(desc.effect == EFFECT_OCEANSPHERE_PATCH);
	assert(desc.dirLights < 5);
	std::stringstream ss;
	ss << stringf("#define NUM_LIGHTS %0{u}\n", desc.dirLights);
	if (desc.atmosphere)
		ss << "#define ATMOSPHERE\n";
	return new Graphics::GL2::OceanSphereProgram("oceansphere_patch", ss.str());
}

void OceanSphereSurfaceMaterial::Apply()
{
	//XXX replace with actual material parameter
	glMaterialfv (GL_FRONT, GL_EMISSION, &emissive[0]);

	SetGSUniforms();
}

void OceanSphereSurfaceMaterial::SetGSUniforms()
{
	OceanSphereProgram *p = static_cast<OceanSphereProgram*>(m_program);
	const SystemBody::AtmosphereParameters ap = *static_cast<SystemBody::AtmosphereParameters*>(this->specialParameter0);

	p->Use();
	p->invLogZfarPlus1.Set(m_renderer->m_invLogZfarPlus1);
	p->sceneAmbient.Set(m_renderer->GetAmbientColor());
	p->atmosColor.Set(ap.atmosCol);
	p->oceansphereAtmosFogDensity.Set(ap.atmosDensity);
	p->oceansphereAtmosInvScaleHeight.Set(ap.atmosInvScaleHeight);
	p->oceansphereAtmosTopRad.Set(ap.atmosRadius);
	p->oceansphereCenter.Set(ap.center);
	p->oceansphereScaledRadius.Set(ap.planetRadius / ap.scale);
	p->oceansphereScale.Set(ap.scale);
}

}
}
