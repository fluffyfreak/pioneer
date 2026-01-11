// Copyright © 2008-2024 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#pragma once

#include "graphics/RenderState.h"
#include "graphics/Renderer.h"
#include "graphics/Types.h"
#include "graphics/UniformBuffer.h"

#include "graphics/VertexArray.h"
#include "graphics/VertexBuffer.h"
#include "graphics/vulkan/VkMaterial.h"
#include "graphics/vulkan/VkRenderTarget.h"
#include "graphics/vulkan/VkTexture.h"
#include "graphics/vulkan/VkUniformBuffer.h"
#include "graphics/vulkan/VkVertexBuffer.h"

#include "RefCounted.h"
#include <stack>
#include <unordered_map>

#include <vulkan/vulkan.hpp>

namespace Graphics {

	class Texture;
	struct Settings;

	namespace Vulkan {
		class CachedVertexBuffer;
		class CommandList;
		class InstanceBuffer;
		class IndexBuffer;
		class Material;
		class MeshObject;
		class RenderState;
		class RenderStateCache;
		class RenderTarget;
		class Shader;
		class UniformBuffer;
		class UniformLinearBuffer;
		class VertexBuffer;
	} // namespace Vulkan

	class VkRenderer final : public Renderer {
	public:
		static void RegisterRenderer();

		VkRenderer(SDL_Window *window, const Graphics::Settings &vs);
		~VkRenderer() final;

		void InitVulkan();

		const char *GetName() const final { return "Vulkan renderer"; }
		RendererType GetRendererType() const final { return RENDERER_VULKAN; }

		bool SupportsInstancing() final { return false; }
		int GetMaximumNumberAASamples() const final { return 0; }
		bool GetNearFarRange(float &near_, float &far_) const final { return true; }

		void SetVSyncEnabled(bool) override {}

		bool BeginFrame() final { return true; }
		bool EndFrame() final { return true; }
		bool SwapBuffers() final { return true; }

		RenderTarget *GetRenderTarget() final { return m_rt; }
		bool SetRenderTarget(RenderTarget *rt) final
		{
			m_rt = rt;
			return true;
		}
		bool SetScissor(ViewportExtents ext) final { return true; }

		void CopyRenderTarget(RenderTarget *, RenderTarget *, ViewportExtents, ViewportExtents, bool) final {}
		void ResolveRenderTarget(RenderTarget *, RenderTarget *, ViewportExtents) final {}

		bool ClearScreen(const Color &, bool) final { return true; }
		bool ClearDepthBuffer() final { return true; }

		bool SetViewport(ViewportExtents v) final { return true; }
		ViewportExtents GetViewport() const final { return {}; }

		bool SetTransform(const matrix4x4f &m) final { return true; }
		matrix4x4f GetTransform() const final { return matrix4x4f::Identity; }
		bool SetPerspectiveProjection(float fov, float aspect, float near_, float far_) final { return true; }
		bool SetOrthographicProjection(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax) final { return true; }
		bool SetProjection(const matrix4x4f &m) final { return true; }
		matrix4x4f GetProjection() const final { return matrix4x4f::Identity; }

		bool SetWireFrameMode(bool enabled) final { return true; }

		bool SetLightIntensity(Uint32, const float *) final { return true; }
		bool SetLights(Uint32 numlights, const Light *l) final { return true; }
		Uint32 GetNumLights() const final { return 1; }
		bool SetAmbientColor(const Color &c) final { return true; }

		bool FlushCommandBuffers() final { return true; }

		bool DrawBuffer(const VertexArray *, Material *) final { return true; }
		bool DrawBufferDynamic(VertexBuffer *, uint32_t, IndexBuffer *, uint32_t, uint32_t, Material *) final { return true; }
		bool DrawMesh(MeshObject *, Material *) final { return true; }

		void Draw(Span<VertexBuffer *const>, IndexBuffer *, Material *m, uint32_t, uint32_t) final {}

		Material *CreateMaterial(const std::string &s, const MaterialDescriptor &d, const RenderStateDesc &rsd, const VertexFormatDesc &vfd) final { return new Graphics::Vulkan::Material(rsd); }
		Material *CloneMaterial(const Material *m, const MaterialDescriptor &d, const RenderStateDesc &rsd, const VertexFormatDesc &vfd) final { return new Graphics::Vulkan::Material(rsd); }
		Texture *CreateTexture(const TextureDescriptor &d) final { return new Graphics::Vulkan::VkTexture(d); }
		RenderTarget *CreateRenderTarget(const RenderTargetDesc &d) final { return new Graphics::Vulkan::RenderTarget(d); }
		VertexBuffer *CreateVertexBuffer(BufferUsage usage, uint32_t numVertices, uint32_t stride) final { return new Graphics::Vulkan::VertexBuffer(usage, numVertices, stride); }
		IndexBuffer *CreateIndexBuffer(Uint32 size, BufferUsage bu, IndexBufferSize el) final { return new Graphics::Vulkan::IndexBuffer(size, bu, el); }
		UniformBuffer *CreateUniformBuffer(Uint32 size, BufferUsage bu) final { return new Graphics::Vulkan::UniformBuffer(size, bu); }
		MeshObject *CreateMeshObject(const VertexFormatDesc &desc, VertexBuffer *v, IndexBuffer *i) final { return new Graphics::Vulkan::MeshObject(desc, static_cast<Vulkan::VertexBuffer *>(v), static_cast<Vulkan::IndexBuffer *>(i)); }
		MeshObject *CreateMeshObjectFromArray(const VertexArray *v, IndexBuffer *i = nullptr, BufferUsage usage = BUFFER_USAGE_STATIC) final
		{
			VertexFormatDesc desc = VertexFormatDesc::FromAttribSet(v->GetAttributeSet());
			Graphics::VertexBuffer *vertexBuffer = CreateVertexBuffer(usage, v->GetNumVerts(), desc.bindings[0].stride);
			v->Populate(vertexBuffer);

			return CreateMeshObject(desc, vertexBuffer, i);
		}

		const RenderStateDesc &GetMaterialRenderState(const Graphics::Material *m) final { return static_cast<const Vulkan::Material *>(m)->rsd; }

		bool ReloadShaders() final { return true; }

	protected:
		void PushState() final {}
		void PopState() final {}

	private:
		Graphics::RenderTarget *m_rt;

		VkInstance m_vkInst;
		VkDevice m_device;
		VkQueue m_graphicsQueue;
		VkQueue m_presentQueue;
	};

} // namespace Graphics
