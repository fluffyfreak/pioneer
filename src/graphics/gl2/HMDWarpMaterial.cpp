// Copyright © 2008-2013 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "HMDWarpMaterial.h"
#include "GeoSphere.h"
#include "Camera.h"
#include "StringF.h"
#include "graphics/Graphics.h"
#include "graphics/TextureGL.h"
#include "graphics/RendererGL2.h"
#include "OculusRift.h"
#include <sstream>

namespace Graphics {
namespace GL2 {

HMDWarpProgram::HMDWarpProgram(const std::string &filename, const std::string &defines)
{
	m_name = filename;
	m_defines = defines;
	LoadShaders(filename, defines);
	InitUniforms();
}

void HMDWarpProgram::InitUniforms()
{
	Program::InitUniforms();
	LensCenter.Init("LensCenter", m_program);
	ScreenCenter.Init("ScreenCenter", m_program);
	Scale.Init("Scale", m_program);
	ScaleIn.Init("ScaleIn", m_program);
	HmdWarpParam.Init("HmdWarpParam", m_program);
	ChromAbParam.Init("ChromAbParam", m_program);
	Texm.Init("Texm", m_program);
}

Program *HMDWarpMaterial::CreateProgram(const MaterialDescriptor &desc)
{
	assert((desc.effect == EFFECT_HMDWARP));
	std::stringstream ss;
	if (desc.textures > 0)
		ss << "#define TEXTURE0\n";
	return new Graphics::GL2::HMDWarpProgram("HMDWarp", ss.str());
}
#pragma optimize("",off)
void HMDWarpMaterial::Apply()
{
	SetGSUniforms();
	glPushAttrib(GL_ENABLE_BIT);
	glDisable(GL_CULL_FACE);
}
#pragma optimize("",off)
void HMDWarpMaterial::Unapply()
{
	glPopAttrib();
	if (texture0) {
		static_cast<TextureGL*>(texture0)->Unbind();
	}
	m_program->Unuse();
}
#pragma optimize("",off)
void HMDWarpMaterial::SetGSUniforms()
{
	HMDWarpProgram *p = static_cast<HMDWarpProgram*>(m_program);

	p->Use();
	p->diffuse.Set(this->diffuse);
	p->texture0.Set(this->texture0, 0);
	p->emission.Set(this->emissive);
	p->sceneAmbient.Set(m_renderer->GetAmbientColor());
	const float WindowWidth = Graphics::GetScreenWidth();
	const float WindowHeight = Graphics::GetScreenHeight();

	assert(specialParameter0);
	OculusRiftInterface::Viewport* VP = static_cast<OculusRiftInterface::Viewport*>(specialParameter0);
	const float w = float(VP->w) / float(WindowWidth),
                h = float(VP->h) / float(WindowHeight),
                x = float(VP->x) / float(WindowWidth),
                y = float(VP->y) / float(WindowHeight);

    const float as = float(VP->w) / float(VP->h);

	float XCenterOffset;
	float Scale;
	float K0,K1,K2,K3;
	OculusRiftInterface::GetDistortionValues(XCenterOffset, Scale, K0, K1, K2, K3);
	if(VP->x > 0) // total hack
		XCenterOffset = -XCenterOffset;

    // We are using 1/4 of DistortionCenter offset value here, since it is
    // relative to [-1,1] range that gets mapped to [0, 0.5].
	p->LensCenter.Set(x + (w + XCenterOffset * 0.5f)*0.5f, y + h*0.5f);
	p->ScreenCenter.Set(x + w*0.5f, y + h*0.5f);

    // MA: This is more correct but we would need higher-res texture vertically; we should adopt this
    // once we have asymmetric input texture scale.
    const float scaleFactor = 1.0f / Scale;

    p->Scale.Set(  (w/2) * scaleFactor, (h/2) * scaleFactor * as);
    p->ScaleIn.Set((2/w),               (2/h) / as);

    p->HmdWarpParam.Set(K0, K1, K2, K3);

    /*if (PostProcessShaderRequested == PostProcessShader_DistortionAndChromAb)
    {
        pPostProcessShader->SetUniform4f("ChromAbParam",
                                        Distortion.ChromaticAberration[0], 
                                        Distortion.ChromaticAberration[1],
                                        Distortion.ChromaticAberration[2],
                                        Distortion.ChromaticAberration[3]);
    }*/

	float texf[16]={w  , 0.f, 0.f, x,
                    0.f, h  , 0.f, y,
                    0.f, 0.f, 0.f, 0.f,
					0.f, 0.f, 0.f, 1.f};
    matrix4x4f texm(texf);
    p->Texm.Set(texm,true);
}

}	// namespace GL2
}	// namespace Graphics
