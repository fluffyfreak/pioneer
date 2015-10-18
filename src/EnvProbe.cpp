// Copyright © 2008-2013 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "EnvProbe.h"
#include "Pi.h"
#include "Player.h"
#include "Easing.h"
#include "FileSystem.h"
#include "PngWriter.h"
#include "graphics/Renderer.h"
#include "graphics/RenderTarget.h"
#include "graphics/TextureBuilder.h"
#include "graphics/opengl/TextureGL.h"
#include "graphics/Graphics.h"
#include "StringF.h"
#include <algorithm>

const size_t EnvProbe::NUM_VIEW_DIRECTIONS = 6;

namespace
{
	static const vector3f pos(0.0f, 0.0f, 0.0f);
	static const vector3f target[EnvProbe::NUM_VIEW_DIRECTIONS] = {
		vector3f(  1.0f,  0.0f,  0.0f ),
		vector3f( -1.0f,  0.0f,  0.0f ),
		vector3f(  0.0f,  1.0f,  0.0f ),
		vector3f(  0.0f, -1.0f,  0.0f ),
		vector3f(  0.0f,  0.0f,  1.0f ),
		vector3f(  0.0f,  0.0f, -1.0f )
	};
	static const vector3f upDir[EnvProbe::NUM_VIEW_DIRECTIONS] = {
		vector3f(  0.0f,  1.0f,  0.0f ),
		vector3f(  0.0f,  1.0f,  0.0f ),
		vector3f(  0.0f,  0.0f, -1.0f ),
		vector3f(  0.0f,  0.0f,  1.0f ),
		vector3f(  0.0f,  1.0f,  0.0f ),
		vector3f(  0.0f,  1.0f,  0.0f )
	};
};

EnvProbe::EnvProbe(Graphics::Renderer *r, const int sizeInPixels)
	: m_sizeInPixels(sizeInPixels)
	, m_renderer(r)
{
	//get near & far clipping distances
	float znear, zfar;
	Pi::renderer->GetNearFarRange(znear, zfar);

	// setup camera context and camera(s)
	m_cameraContexts.resize(EnvProbe::NUM_VIEW_DIRECTIONS);
	m_cameras.resize(EnvProbe::NUM_VIEW_DIRECTIONS);
	for (size_t i = 0; i < EnvProbe::NUM_VIEW_DIRECTIONS; i++) {
		RefCountedPtr<CameraContext> pContext( new CameraContext(m_sizeInPixels, m_sizeInPixels, 90.f, znear, zfar) );
		m_cameraContexts[i] = pContext;
		m_cameras[i].reset(new Camera(pContext, Pi::renderer));
	}

	// setup render targets
	for (int i = 0; i<EnvProbe::NUM_VIEW_DIRECTIONS; i++)
	{
		Graphics::TextureDescriptor texDesc(
			Graphics::TEXTURE_RGB_888,
			vector2f(m_sizeInPixels, m_sizeInPixels),
			Graphics::LINEAR_CLAMP, false, false, 0);
		Graphics::Texture* pTexture = Pi::renderer->CreateTexture(texDesc);

		// Complete the RT description so we can request a buffer.
		// NB: we don't want it to create use a texture because we share it with the textured quad created above.
		Graphics::RenderTargetDesc rtDesc(
			m_sizeInPixels,
			m_sizeInPixels,
			Graphics::TEXTURE_NONE,    // don't create a texture
			Graphics::TEXTURE_DEPTH,
			false);

		Graphics::RenderTarget* pRTarget = Pi::renderer->CreateRenderTarget(rtDesc);

		pRTarget->SetColorTexture(pTexture);

		m_renderTargets.push_back(pRTarget);
	}
}

EnvProbe::~EnvProbe()
{
	for (auto rt : m_renderTargets) {
		delete rt;
	}
}

///////////////////////////////////////////////////////////////////////////////
// compute transform axis from object position, target and up direction
///////////////////////////////////////////////////////////////////////////////
void lookAtToAxes(const vector3f& pos, const vector3f& target, const vector3f& upDir, vector3f& left, vector3f& up, vector3f& forward)
{
    // compute the forward vector
    forward = (target - pos).Normalized();

    // compute the left vector
    left = upDir.Cross(forward).Normalized();  // cross product

    // compute the orthonormal up vector
    up = forward.Cross(left).Normalized();     // cross product
}

///////////////////////////////////////////////////////////////////////////////
// compute transform axis from object position, target and up direction
///////////////////////////////////////////////////////////////////////////////
void lookAtMatrix(const vector3f& pos, const vector3f& target, const vector3f& upDir, matrix4x4f& matOut)
{
    vector3f forward, left, up;
	lookAtToAxes(pos, target, upDir, left, up, forward);

	matOut.LoadFrom3x3Matrix( &matrix3x3f::FromVectors(left, up, forward)[0] );
}

void lookAtMatrix(const vector3f& pos, const vector3f& target, const vector3f& upDir, matrix3x3f& matOut)
{
	vector3f forward, left, up;
	lookAtToAxes(pos, target, upDir, left, up, forward);
	matOut = matrix3x3f::FromVectors(left, up, forward);
}

void EnvProbe::BeginRenderTarget(Graphics::RenderTarget* pTarget) 
{
	assert(pTarget);
	const bool bTargetSet = m_renderer->SetRenderTarget(pTarget);
	assert(bTargetSet);
}

void EnvProbe::EndRenderTarget() 
{
	m_renderer->SetRenderTarget(NULL);
}

void EnvProbe::Draw()
{
	m_renderer->SetViewport(0, 0, m_sizeInPixels, m_sizeInPixels);
	m_renderer->SetPerspectiveProjection(90.0f, 1.f, 1.f, 10000.f);
	m_renderer->SetTransform(matrix4x4f::Identity());

	Graphics::RenderStateDesc rsd;
	rsd.depthTest = true;
	rsd.depthWrite = true;
	m_renderer->SetRenderState(m_renderer->CreateRenderState(rsd));

	const Color oldSceneAmbientColor = m_renderer->GetAmbientColor();
	m_renderer->SetAmbientColor(Color(0.f));
	
	Graphics::TextureCubeData tcd;
	memset(&tcd, 0, sizeof(Graphics::TextureCubeData));

	for (int c = 0; c < EnvProbe::NUM_VIEW_DIRECTIONS; c++)
	{
		// current  pointers
		Graphics::RenderTarget *it = m_renderTargets[c];
		CameraContext *pCamContext = m_cameraContexts[c].Get();
		Camera *pCamera = m_cameras[c].get();

		matrix3x3f matOut;
		lookAtMatrix(pos, target[c], upDir[c], matOut);
		matrix3x3d transform;
		matrix3x3ftod(matOut, transform);

		pCamContext->SetPosition(Pi::player->GetPosition());
		pCamContext->SetOrient(Pi::player->GetOrient() * transform);

		// apply camera settings
		pCamContext->ApplyDrawTransforms(m_renderer);
		pCamContext->BeginFrame();
		pCamera->Update();

		// render state
		Graphics::Renderer::StateTicket ticket(m_renderer);
		m_renderer->PushMatrix();
		
		// render each face of the cubemap
		BeginRenderTarget(it);
		{
			// don't ever draw the player when drawing from the players location.
			pCamera->Draw(Pi::player);
		}
		EndRenderTarget();

		Graphics::Texture* pTex = it->GetColorTexture();

		switch (c)
		{
		case 0: tcd.posX = pTex; break;
		case 1: tcd.negX = pTex; break;
		case 2: tcd.posY = pTex; break;
		case 3: tcd.negY = pTex; break;
		case 4: tcd.posZ = pTex; break;
		case 5: tcd.negZ = pTex; break;
		default:
			assert(false);
			break;
		}
#if 0
		Graphics::TextureGL* pGL = static_cast<Graphics::TextureGL*>(it->GetColorTexture());
		// pad rows to 4 bytes, which is the default row alignment for OpenGL
		const int stride = (3* m_sizeInPixels + 3) & ~3;

		std::vector<Uint8> pixel_data(stride * m_sizeInPixels);
		pGL->Bind();
		glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_UNSIGNED_BYTE, &pixel_data[0]);
		pGL->Unbind();
		glFinish();

		const std::string dir = "cubemaps";
		FileSystem::userFiles.MakeDirectory(dir);
		const std::string destFile = stringf("cube_face_%0{d}.png", c);
		const std::string fname = FileSystem::JoinPathBelow(dir, destFile);

		write_png(FileSystem::userFiles, fname, &pixel_data[0], m_sizeInPixels, m_sizeInPixels, stride, 3);

		printf("cubemap %s saved\n", fname.c_str());
#endif

		m_renderer->PopMatrix();
		pCamContext->EndFrame();
	}

	m_renderer->SetAmbientColor(oldSceneAmbientColor);
}
