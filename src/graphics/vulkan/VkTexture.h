// Copyright © 2008-2024 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#pragma once

#ifndef _VKTEXTURE_H
#define _VKTEXTURE_H

#include "graphics/Texture.h"

namespace Graphics {

	namespace Vulkan {

		class VkTexture : public Texture {
		public:
			VkTexture(const TextureDescriptor &descriptor) :
				Texture(descriptor) {}

			void Update(const void *data, const vector2f &pos, const vector3f &dataSize, TextureFormat format, const unsigned int numMips) final {}
			void Update(const TextureCubeData &data, const vector3f &dataSize, TextureFormat format, const unsigned int numMips) final {}
			void Update(const vecDataPtr &data, const vector3f &dataSize, const TextureFormat format, const unsigned int numMips) final {}

			void Bind() final {}
			void Unbind() final {}

			void SetSampleMode(TextureSampleMode) final {}
			void BuildMipmaps(const uint32_t) final {}
			uint32_t GetTextureID() const final { return 0U; }
			uint32_t GetTextureMemSize() const final { return 0U; }

		private:
		};

	} // namespace Vulkan

} // namespace Graphics

#endif
