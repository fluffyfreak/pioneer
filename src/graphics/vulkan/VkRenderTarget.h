// Copyright © 2008-2024 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#pragma once

#ifndef _VKRENDERTARGET_H
#define _VKRENDERTARGET_H

#include "graphics/RenderTarget.h"

namespace Graphics {

	class VkRenderer;

	namespace Vulkan {

		class RenderTarget : public Graphics::RenderTarget {
		public:
			virtual Texture *GetColorTexture() const { return m_colorTexture.Get(); }
			virtual Texture *GetDepthTexture() const { return m_depthTexture.Get(); }
			virtual void SetCubeFaceTexture(const Uint32 face, Texture *t) final { m_colorTexture.Reset(t); }
			virtual void SetColorTexture(Texture *t) { m_colorTexture.Reset(t); }
			virtual void SetDepthTexture(Texture *t) { m_depthTexture.Reset(t); }

		protected:
			friend class Graphics::VkRenderer;
			RenderTarget(const RenderTargetDesc &d) :
				Graphics::RenderTarget(d) {}

		private:
			RefCountedPtr<Texture> m_colorTexture;
			RefCountedPtr<Texture> m_depthTexture;
		};

	} // namespace Vulkan

} // namespace Graphics

#endif
