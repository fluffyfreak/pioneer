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
		virtual ~VkRenderer() override final;

		void InitVulkan();

		virtual const char *GetName() const override final { return "Vulkan renderer"; }
		virtual RendererType GetRendererType() const override final { return RENDERER_VULKAN; }

		virtual bool SupportsInstancing() override final { return false; }
		virtual int GetMaximumNumberAASamples() const override final { return 0; }
		virtual bool GetNearFarRange(float &near_, float &far_) const override final { return true; }

		virtual void SetVSyncEnabled(bool) override {}

		virtual bool BeginFrame() override final { return true; }
		virtual bool EndFrame() override final { return true; }
		virtual bool SwapBuffers() override final { return true; }

		virtual RenderTarget *GetRenderTarget() override final { return m_rt; }
		virtual bool SetRenderTarget(RenderTarget *rt) override final
		{
			m_rt = rt;
			return true;
		}
		virtual bool SetScissor(ViewportExtents ext) override final { return true; }

		virtual void CopyRenderTarget(RenderTarget *, RenderTarget *, ViewportExtents, ViewportExtents, bool) override final {}
		virtual void ResolveRenderTarget(RenderTarget *, RenderTarget *, ViewportExtents) override final {}

		virtual bool ClearScreen(const Color &, bool) override final { return true; }
		virtual bool ClearDepthBuffer() override final { return true; }

		virtual bool SetViewport(ViewportExtents v) override final { return true; }
		virtual ViewportExtents GetViewport() const override final { return {}; }

		virtual bool SetTransform(const matrix4x4f &m) override final { return true; }
		virtual matrix4x4f GetTransform() const override final { return matrix4x4f::Identity(); }
		virtual bool SetPerspectiveProjection(float fov, float aspect, float near_, float far_) override final { return true; }
		virtual bool SetOrthographicProjection(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax) override final { return true; }
		virtual bool SetProjection(const matrix4x4f &m) override final { return true; }
		virtual matrix4x4f GetProjection() const override final { return matrix4x4f::Identity(); }

		virtual bool SetWireFrameMode(bool enabled) override final { return true; }

		virtual bool SetLightIntensity(Uint32, const float *) override final { return true; }
		virtual bool SetLights(Uint32 numlights, const Light *l) override final { return true; }
		virtual Uint32 GetNumLights() const override final { return 1; }
		virtual bool SetAmbientColor(const Color &c) override final { return true; }

		virtual bool FlushCommandBuffers() override final { return true; }

		virtual bool DrawBuffer(const VertexArray *, Material *) override final { return true; }
		virtual bool DrawBufferDynamic(VertexBuffer *, uint32_t, IndexBuffer *, uint32_t, uint32_t, Material *) override final { return true; }
		virtual bool DrawMesh(MeshObject *, Material *) override final { return true; }
		virtual bool DrawMeshInstanced(MeshObject *, Material *, InstanceBuffer *) override final { return true; }

		virtual Material *CreateMaterial(const std::string &s, const MaterialDescriptor &d, const RenderStateDesc &rsd) override final { return new Graphics::Vulkan::Material(rsd); }
		virtual Material *CloneMaterial(const Material *m, const MaterialDescriptor &d, const RenderStateDesc &rsd) override final { return new Graphics::Vulkan::Material(rsd); }
		virtual Texture *CreateTexture(const TextureDescriptor &d) override final { return new Graphics::Vulkan::VkTexture(d); }
		virtual RenderTarget *CreateRenderTarget(const RenderTargetDesc &d) override final { return new Graphics::Vulkan::RenderTarget(d); }
		virtual VertexBuffer *CreateVertexBuffer(const VertexBufferDesc &d) override final { return new Graphics::Vulkan::VertexBuffer(d); }
		virtual IndexBuffer *CreateIndexBuffer(Uint32 size, BufferUsage bu, IndexBufferSize el) override final { return new Graphics::Vulkan::IndexBuffer(size, bu, el); }
		virtual InstanceBuffer *CreateInstanceBuffer(Uint32 size, BufferUsage bu) override final { return new Graphics::Vulkan::InstanceBuffer(size, bu); }
		virtual UniformBuffer *CreateUniformBuffer(Uint32 size, BufferUsage bu) override final { return new Graphics::Vulkan::UniformBuffer(size, bu); }
		virtual MeshObject *CreateMeshObject(VertexBuffer *v, IndexBuffer *i) override final { return new Graphics::Vulkan::MeshObject(static_cast<Vulkan::VertexBuffer *>(v), static_cast<Vulkan::IndexBuffer *>(i)); }
		virtual MeshObject *CreateMeshObjectFromArray(const VertexArray *v, IndexBuffer *i = nullptr, BufferUsage = BUFFER_USAGE_STATIC) override final
		{
			auto desc = Graphics::VertexBufferDesc::FromAttribSet(v->GetAttributeSet());
			desc.numVertices = v->GetNumVerts();
			return new Graphics::Vulkan::MeshObject(static_cast<Vulkan::VertexBuffer *>(CreateVertexBuffer(desc)), static_cast<Vulkan::IndexBuffer *>(i));
		}

		virtual const RenderStateDesc &GetMaterialRenderState(const Graphics::Material *m) override final { return static_cast<const Vulkan::Material *>(m)->rsd; }

		virtual bool ReloadShaders() override final { return true; }

	protected:
		virtual void PushState() override final {}
		virtual void PopState() override final {}

	private:
		const matrix4x4f m_identity;
		Graphics::RenderTarget *m_rt;

		VkInstance m_vkInst;
		VkDevice m_device;
		VkQueue m_graphicsQueue;
		VkQueue m_presentQueue;
	};

} // namespace Graphics
