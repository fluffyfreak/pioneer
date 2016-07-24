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

#define DUMP_TO_TEXTURE 0

#if DUMP_TO_TEXTURE
#include "FileSystem.h"
#include "PngWriter.h"
#include "graphics/opengl/TextureGL.h"
void textureDump(const char* destFile, const int width, const int height, const Color* buf)
{
	PROFILE_SCOPED()
	const std::string dir = "cubemaps";
	FileSystem::userFiles.MakeDirectory(dir);
	const std::string fname = FileSystem::JoinPathBelow(dir, destFile);

	// pad rows to 4 bytes, which is the default row alignment for OpenGL
	//const int stride = (3*width + 3) & ~3;
	const int stride = width * 4;

	write_png(FileSystem::userFiles, fname, &buf[0].r, width, height, stride, 4);

	printf("texture %s saved\n", fname.c_str());
}
#endif

const size_t EnvProbe::NUM_VIEW_DIRECTIONS = 6;

namespace
{
	static const vector3f pos(0.0f, 0.0f, 0.0f);
	static const vector3f target[EnvProbe::NUM_VIEW_DIRECTIONS] = {
		vector3f(  1.0f,  0.0f,  0.0f ),
		vector3f( -1.0f,  0.0f,  0.0f ),
		vector3f(  0.0f, -1.0f,  0.0f ),
		vector3f(  0.0f,  1.0f,  0.0f ),
		vector3f(  0.0f,  0.0f,  1.0f ),
		vector3f(  0.0f,  0.0f, -1.0f )
	};
	static const vector3f upDir[EnvProbe::NUM_VIEW_DIRECTIONS] = {
		vector3f(  0.0f, -1.0f,  0.0f ),
		vector3f(  0.0f, -1.0f,  0.0f ),
		vector3f(  0.0f,  0.0f, -1.0f ),
		vector3f(  0.0f,  0.0f,  1.0f ),
		vector3f(  0.0f, -1.0f,  0.0f ),
		vector3f(  0.0f, -1.0f,  0.0f )
	};
	static const Color4f colorDir[EnvProbe::NUM_VIEW_DIRECTIONS] = {
		Color4f(  Color4f::RED ),
		Color4f(  Color4f::RED * 0.25f ),
		Color4f(  Color4f::GREEN ),
		Color4f(  Color4f::GREEN * 0.25f ),
		Color4f(  Color4f::BLUE ),
		Color4f(  Color4f::BLUE * 0.25f )
	};
};

EnvProbe::EnvProbe(Graphics::Renderer *r, const int sizeInPixels)
	: m_sizeInPixels(sizeInPixels)
	, m_renderer(r)
{
	PROFILE_SCOPED()
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

	Graphics::TextureDescriptor texDesc(
		Graphics::TEXTURE_RGBA_8888,
		vector2f(m_sizeInPixels, m_sizeInPixels),
		Graphics::LINEAR_CLAMP, false, false, false, 0, Graphics::TEXTURE_CUBE_MAP);
	m_cubemap.Reset( Pi::renderer->CreateTexture(texDesc) );

	// setup render target
	// Complete the RT description so we can request a buffer.
	// NB: we don't want it to create use a texture because we share it with the textured quad created above.
	Graphics::RenderTargetDesc rtDesc(
		m_sizeInPixels,
		m_sizeInPixels,
		Graphics::TEXTURE_NONE,    // don't create a texture
		Graphics::TEXTURE_DEPTH,
		false);

	m_renderTarget.reset( Pi::renderer->CreateRenderTarget(rtDesc) );
}

EnvProbe::~EnvProbe()
{
}

///////////////////////////////////////////////////////////////////////////////
// compute transform axis from object position, target and up direction
///////////////////////////////////////////////////////////////////////////////
void lookAtToAxes(const vector3f& pos, const vector3f& target, const vector3f& upDir, vector3f& left, vector3f& up, vector3f& forward)
{
	PROFILE_SCOPED()
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
/*void lookAtMatrix(const vector3f& pos, const vector3f& target, const vector3f& upDir, matrix4x4f& matOut)
{
	PROFILE_SCOPED()
    vector3f forward, left, up;
	lookAtToAxes(pos, target, upDir, left, up, forward);

	matOut.LoadFrom3x3Matrix( &matrix3x3f::FromVectors(left, up, forward)[0] );
}*/

void lookAtMatrix(const vector3f& pos, const vector3f& target, const vector3f& upDir, matrix3x3f& matOut)
{
	PROFILE_SCOPED()
	vector3f forward, left, up;
	lookAtToAxes(pos, target, upDir, left, up, forward);
	matOut = matrix3x3f::FromVectors(left, up, forward);
}

void EnvProbe::BeginRenderTarget(Graphics::RenderTarget* pTarget) 
{
	PROFILE_SCOPED()
	assert(pTarget);
	const bool bTargetSet = m_renderer->SetRenderTarget(pTarget);
	assert(bTargetSet);
}

void EnvProbe::EndRenderTarget() 
{
	PROFILE_SCOPED()
	m_renderer->SetRenderTarget(NULL);
}

static bool s_hacky_clear = false;
void EnvProbe::Draw(Body *pBody)
{
	PROFILE_SCOPED()
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
		m_cameraContexts[c]->SetFrame(pBody->GetFrame());
		CameraContext *pCamContext = m_cameraContexts[c].Get();
		Camera *pCamera = m_cameras[c].get();

		matrix3x3f matOut;
		lookAtMatrix(pos, target[c], upDir[c], matOut);
		matrix3x3d transform;
		matrix3x3ftod(matOut, transform);

		pCamContext->SetPosition(pBody->GetPosition());
		pCamContext->SetOrient(pBody->GetOrient() * transform);

		// apply camera settings
		pCamContext->ApplyDrawTransforms(m_renderer);
		pCamContext->BeginFrame();
		pCamera->Update();

		// render state
		Graphics::Renderer::StateTicket ticket(m_renderer);
		m_renderer->PushMatrix();

		m_renderTarget->SetCubeFaceTexture(c, m_cubemap.Get());
		
		// render each face of the cubemap
		BeginRenderTarget(m_renderTarget.get());
		{
			// don't ever draw the body itself.
			if(s_hacky_clear) {
				Pi::renderer->SetClearColor(colorDir[c]);
				Pi::renderer->ClearScreen();
			} else {
				pCamera->Draw(pBody);
			}
		}
		EndRenderTarget();

		m_renderer->PopMatrix();
		pCamContext->EndFrame();
	}

#if DUMP_TO_TEXTURE
	static int s_count = 0;
	Graphics::TextureGL* pGLTex = static_cast<Graphics::TextureGL*>(m_renderTarget->GetColorTexture());
	for (int c = 0; c < EnvProbe::NUM_VIEW_DIRECTIONS; c++)
	{
		// current  pointers
		std::unique_ptr<Color, FreeDeleter> buffer(static_cast<Color*>(malloc(m_sizeInPixels*m_sizeInPixels*4)));
		
		pGLTex->Bind();
		glGetTexImage(GL_TEXTURE_CUBE_MAP_POSITIVE_X + c, 0, GL_RGBA, GL_UNSIGNED_BYTE, buffer.get());
		pGLTex->Unbind();

		const std::string destFile = stringf("cube_face_%0{d}_%1{d}.png", s_count, c);
		textureDump(destFile.c_str(), m_sizeInPixels, m_sizeInPixels, buffer.get());
	}
	++s_count;
#endif

	m_renderer->SetAmbientColor(oldSceneAmbientColor);
}
