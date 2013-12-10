// Copyright © 2008-2013 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "BackgroundGen.h"
#include "Pi.h"
#include "Easing.h"
#include "FileSystem.h"
#include "PngWriter.h"
#include "graphics/Renderer.h"
#include "graphics/RenderTarget.h"
#include "graphics/TextureBuilder.h"
#include "graphics/TextureGL.h"
#include "graphics/Graphics.h"
#include "StringF.h"
#include <algorithm>

namespace
{
	static const vector3f pos(0.0f, 0.0f, 0.0f);
	static const vector3f target[6] = {
		vector3f(  1.0f,  0.0f,  0.0f ),
		vector3f( -1.0f,  0.0f,  0.0f ),
		vector3f(  0.0f,  1.0f,  0.0f ),
		vector3f(  0.0f, -1.0f,  0.0f ),
		vector3f(  0.0f,  0.0f,  1.0f ),
		vector3f(  0.0f,  0.0f, -1.0f )
	};
	static const vector3f upDir[6] = {
		vector3f(  0.0f,  1.0f,  0.0f ),
		vector3f(  0.0f,  1.0f,  0.0f ),
		vector3f(  0.0f,  0.0f, -1.0f ),
		vector3f(  0.0f,  0.0f,  1.0f ),
		vector3f(  0.0f,  1.0f,  0.0f ),
		vector3f(  0.0f,  1.0f,  0.0f )
	};
};

BackgroundGen::BackgroundGen(Graphics::Renderer *r, int width, int height)
: m_width(width)
, m_height(height)
, m_aspectRatio(float(width)/float(height))
, m_renderer(r)
{
	m_background.reset(new Background::Container(r, UNIVERSE_SEED));
}

BackgroundGen::~BackgroundGen()
{
}

///////////////////////////////////////////////////////////////////////////////
// compute transform axis from object position, target and up direction
///////////////////////////////////////////////////////////////////////////////
void lookAtToAxes(const vector3f& pos, const vector3f& target, const vector3f& upDir,
                  vector3f& left, vector3f& up, vector3f& forward)
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
void lookAtMatrix(const vector3f& pos, const vector3f& target, const vector3f& upDir,
                  matrix4x4f& matOut)
{
    vector3f forward, left, up;
	lookAtToAxes(pos, target, upDir, left, up, forward);

	matOut.LoadFrom3x3Matrix( &matrix3x3f::FromVectors(left, up, forward)[0] );
}

static void BeginRenderTarget(Graphics::Renderer *r, Graphics::RenderTarget* pTarget) {
	assert(pTarget);
	const bool bTargetSet = r->SetRenderTarget(pTarget);
	assert(bTargetSet);
}

static void EndRenderTarget(Graphics::Renderer *r) {
	r->SetRenderTarget(NULL);
}

typedef std::pair<Graphics::RenderTarget*, Graphics::Texture*> RTTexPair;

void BackgroundGen::Draw()
{
	// setup render targets
	std::vector<RTTexPair> RTtargets;
	for( int i=0; i<6; i++ )
	{
		Graphics::TextureDescriptor texDesc(
			Graphics::TEXTURE_RGB_888,
			vector2f(m_width, m_height),
			Graphics::LINEAR_CLAMP, false, false, 0);
		Graphics::Texture* pTexture = Pi::renderer->CreateTexture(texDesc);

		// Complete the RT description so we can request a buffer.
		// NB: we don't want it to create use a texture because we share it with the textured quad created above.
		Graphics::RenderTargetDesc rtDesc(
			m_width,
			m_height,
			Graphics::TEXTURE_NONE,    // don't create a texture
			Graphics::TEXTURE_DEPTH,
			false);
		
		Graphics::RenderTarget* pRTarget = Pi::renderer->CreateRenderTarget(rtDesc);

		pRTarget->SetColorTexture(pTexture);

		RTtargets.push_back(std::make_pair(pRTarget, pTexture));
	}

	m_renderer->SetViewport(0, 0, m_width, m_height);
	m_renderer->SetPerspectiveProjection(90.0f, m_aspectRatio, 1.f, 10000.f);
	m_renderer->SetTransform(matrix4x4f::Identity());

	m_renderer->SetDepthTest(true);
	m_renderer->SetDepthWrite(true);

	glPushAttrib(GL_ALL_ATTRIB_BITS & (~GL_POINT_BIT));

	const Color oldSceneAmbientColor = m_renderer->GetAmbientColor();
	m_renderer->SetAmbientColor(Color(0.f));
	
    matrix4x4f matOut;
	Graphics::TextureCubeData tcd;

	int i=0;
	for (auto it = RTtargets.begin(), itEnd = RTtargets.end(); it != itEnd; ++it)
	{
		glPushMatrix();
		lookAtMatrix(pos, target[i], upDir[i], matOut);
		matrix4x4d transform;
		matrix4x4ftod(matOut, transform);
		
		// render each face of the cubemap
		BeginRenderTarget(m_renderer, (*it).first);
		{
			m_background->SetDrawFlags( Background::Container::DRAW_STARS | Background::Container::DRAW_MILKY );
			m_background->Draw(m_renderer, transform);
		}
		EndRenderTarget(m_renderer);

		Graphics::Texture* pGL = (*it).first->GetColorTexture();

		switch (i)
		{
		case 0: tcd.posX = pGL; break;
		case 1: tcd.negX = pGL; break;
		case 2: tcd.posY = pGL; break;
		case 3: tcd.negY = pGL; break;
		case 4: tcd.posZ = pGL; break;
		case 5: tcd.negZ = pGL; break;
		default:
			assert(false);
			break;
		}
#if 0
		Graphics::TextureGL* pGL = static_cast<Graphics::TextureGL*>((*it).first->GetColorTexture());
		// pad rows to 4 bytes, which is the default row alignment for OpenGL
		const int stride = (3*m_width + 3) & ~3;

		std::vector<Uint8> pixel_data(stride * m_height);
		pGL->Bind();
		glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_UNSIGNED_BYTE, &pixel_data[0]);
		pGL->Unbind();
		glFinish();

		const std::string dir = "cubemaps";
		FileSystem::userFiles.MakeDirectory(dir);
		const std::string destFile = stringf("cube_face_%0{d}.png", i);
		const std::string fname = FileSystem::JoinPathBelow(dir, destFile);

		write_png(FileSystem::userFiles, fname, &pixel_data[0], m_width, m_height, stride, 3);

		printf("cubemap %s saved\n", fname.c_str());
#endif
		glPopMatrix();
		++i;
	}

	m_renderer->SetAmbientColor(oldSceneAmbientColor);

	glPopAttrib();
}
