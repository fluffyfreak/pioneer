// Copyright © 2008-2024 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#pragma once

#ifndef _VKMATERIAL_H
#define _VKMATERIAL_H

#include "graphics/Material.h"
#include "graphics/RenderState.h"
#include "graphics/Renderer.h"

namespace Graphics {

	class VkRenderer;

	namespace Vulkan {

		class Program;

		class Material : public Graphics::Material {
		public:
			Material(RenderStateDesc rsd) :
				rsd(rsd) {}

			// Create an appropriate program for this material.
			virtual Program *CreateProgram(const MaterialDescriptor &) { return nullptr; }
			virtual bool IsProgramLoaded() const override final { return false; }
			virtual void SetProgram(Program *p) {}

			virtual bool SetTexture(size_t name, Texture *tex) override final { return false; }
			virtual bool SetBuffer(size_t hash, BufferBinding<Graphics::UniformBuffer> uboBinding) override final { return false; }
			virtual bool SetBufferDynamic(size_t name, void *data, size_t size) override { return false; }

			virtual bool SetPushConstant(size_t name, int i) override final { return false; }
			virtual bool SetPushConstant(size_t name, float f) override final { return false; }
			virtual bool SetPushConstant(size_t name, vector3f v3) override final { return false; }
			virtual bool SetPushConstant(size_t name, vector3f v4, float f4) override final { return false; }
			virtual bool SetPushConstant(size_t name, Color c) override final { return false; }
			virtual bool SetPushConstant(size_t name, matrix3x3f mat3) override final { return false; }
			virtual bool SetPushConstant(size_t name, matrix4x4f mat4) override final { return false; }

			RenderStateDesc rsd; // here to ensure validation works correctly
		};
	} // namespace Vulkan
} // namespace Graphics

#endif
