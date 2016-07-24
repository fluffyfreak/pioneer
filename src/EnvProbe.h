// Copyright © 2008-2013 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _ENVPROBE_H
#define _ENVPROBE_H

#include "Camera.h"
#include "CameraController.h"
#include "Background.h"

class EnvProbe {
public:
	EnvProbe(Graphics::Renderer *r, const int sizeInPixels);
	~EnvProbe();
	void Draw(Body*);

	static const size_t NUM_VIEW_DIRECTIONS;

	Graphics::Texture* GetCubemap() const { return m_cubemap.Get(); }

private:
	void BeginRenderTarget(Graphics::RenderTarget* pTarget);
	void EndRenderTarget();

	int m_sizeInPixels;
	Graphics::Renderer *m_renderer;
	RefCountedPtr<Graphics::Texture> m_cubemap;
	std::unique_ptr<Graphics::RenderTarget> m_renderTarget;
	std::vector<RefCountedPtr<CameraContext>> m_cameraContexts;
	std::vector<std::unique_ptr<Camera>> m_cameras;
};

#endif
