// Copyright © 2008-2024 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "VkRenderer.h"

namespace Graphics {

	static Renderer *CreateRenderer(const Settings &vs)
	{
		Graphics::Settings gs;
		return new VkRenderer(nullptr, gs);
	}

	void VkRenderer::RegisterRenderer()
	{
		Graphics::RegisterRenderer(Graphics::RENDERER_VULKAN, CreateRenderer);
	}

} // namespace Graphics
