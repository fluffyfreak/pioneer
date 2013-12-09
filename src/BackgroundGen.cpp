// Copyright © 2008-2013 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "BackgroundGen.h"
#include "Pi.h"
#include "Easing.h"
#include "graphics/Renderer.h"
#include "graphics/RenderTarget.h"
#include "graphics/TextureBuilder.h"
#include "graphics/Graphics.h"
#include <algorithm>

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

static void BeginRenderTarget(Graphics::Renderer *r, Graphics::RenderTarget* pTarget) {
	assert(pTarget);
	const bool bTargetSet = r->SetRenderTarget(pTarget);
	assert(bTargetSet);
}

static void EndRenderTarget(Graphics::Renderer *r) {
	r->SetRenderTarget(NULL);
}

typedef std::pair<Graphics::RenderTarget*, Graphics::Texture*> RTTexPair;

void BackgroundGen::Draw(float _time)
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

	for (auto it = RTtargets.begin(), itEnd = RTtargets.end(); it != itEnd; ++it)
	{
		// render each face of the cubemap
		BeginRenderTarget(m_renderer, (*it).first);
		{
			matrix4x4d brot = matrix4x4d::RotateXMatrix(-0.25*_time) * matrix4x4d::RotateZMatrix(0.6);
			m_background->SetDrawFlags( Background::Container::DRAW_STARS | Background::Container::DRAW_MILKY );
			m_background->Draw(m_renderer, brot);
		}
		EndRenderTarget(m_renderer);
	}

	m_renderer->SetAmbientColor(oldSceneAmbientColor);

	glPopAttrib();
}
