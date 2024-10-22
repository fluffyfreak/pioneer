// Copyright © 2008-2024 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "VkRenderer.h"

#include "MathUtil.h"
#include "RefCounted.h"
#include "SDL_stdinc.h"
#include "StringF.h"

#include "graphics/Graphics.h"
#include "graphics/Light.h"
#include "graphics/Material.h"
#include "graphics/RenderState.h"
#include "graphics/Texture.h"
#include "graphics/TextureBuilder.h"
#include "graphics/Types.h"
#include "graphics/VertexArray.h"
#include "graphics/VertexBuffer.h"

#include "profiler/Profiler.h"

#include "core/Log.h"

#include <SDL.h>
#include <SDL_vulkan.h>

namespace Graphics {

	constexpr int MAX_FRAMES_IN_FLIGHT = 2;

	std::string VkErrorString(VkResult errorCode)
	{
		switch (errorCode) {
#define STR(r) \
	case VK_##r: return #r
			STR(NOT_READY);
			STR(TIMEOUT);
			STR(EVENT_SET);
			STR(EVENT_RESET);
			STR(INCOMPLETE);
			STR(ERROR_OUT_OF_HOST_MEMORY);
			STR(ERROR_OUT_OF_DEVICE_MEMORY);
			STR(ERROR_INITIALIZATION_FAILED);
			STR(ERROR_DEVICE_LOST);
			STR(ERROR_MEMORY_MAP_FAILED);
			STR(ERROR_LAYER_NOT_PRESENT);
			STR(ERROR_EXTENSION_NOT_PRESENT);
			STR(ERROR_FEATURE_NOT_PRESENT);
			STR(ERROR_INCOMPATIBLE_DRIVER);
			STR(ERROR_TOO_MANY_OBJECTS);
			STR(ERROR_FORMAT_NOT_SUPPORTED);
			STR(ERROR_SURFACE_LOST_KHR);
			STR(ERROR_NATIVE_WINDOW_IN_USE_KHR);
			STR(SUBOPTIMAL_KHR);
			STR(ERROR_OUT_OF_DATE_KHR);
			STR(ERROR_INCOMPATIBLE_DISPLAY_KHR);
			STR(ERROR_VALIDATION_FAILED_EXT);
			STR(ERROR_INVALID_SHADER_NV);
			STR(ERROR_INCOMPATIBLE_SHADER_BINARY_EXT);
#undef STR
		default:
			return "UNKNOWN_ERROR";
		}
	}

	std::string VkPhysicalDeviceTypeString(VkPhysicalDeviceType type)
	{
		switch (type) {
#define STR(r) \
	case VK_PHYSICAL_DEVICE_TYPE_##r: return #r
			STR(OTHER);
			STR(INTEGRATED_GPU);
			STR(DISCRETE_GPU);
			STR(VIRTUAL_GPU);
			STR(CPU);
#undef STR
		default: return "UNKNOWN_DEVICE_TYPE";
		}
	}

	static bool LoadLibraryAndCreateWindow(const char *name, const Graphics::Settings &vs, SDL_Window *&m_window)
	{
		PROFILE_SCOPED()

		if (0 != SDL_Vulkan_LoadLibrary(nullptr))
		{
			Error("SDL_Vulkan_LoadLibrary Failed");
			return false;
		}

		Uint32 winFlags = SDL_WINDOW_VULKAN;

		winFlags |= (vs.hidden ? SDL_WINDOW_HIDDEN : SDL_WINDOW_SHOWN);
		if (!vs.hidden && vs.fullscreen) // TODO: support for borderless fullscreen and changing m_window size
			winFlags |= SDL_WINDOW_FULLSCREEN;

		if (vs.canBeResized)
			winFlags |= SDL_WINDOW_RESIZABLE;

		m_window = SDL_CreateWindow(name, SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, vs.width, vs.height, winFlags);
		if (!m_window)
			return false;

		return true;
	}

	static Renderer *CreateRenderer(const Settings &vs)
	{
		PROFILE_SCOPED()
		assert(vs.rendererType == Graphics::RendererType::RENDERER_VULKAN);

		const std::string name("Pioneer");
		SDL_Window *m_window = nullptr;

		bool ok = LoadLibraryAndCreateWindow(name.c_str(), vs, m_window);
		if (!ok) {
			Error("Failed to set video mode: %s", SDL_GetError());
			return nullptr;
		}

		SDLSurfacePtr iconSurface = LoadSurfaceFromFile(vs.iconFile);
		if (iconSurface)
			SDL_SetWindowIcon(m_window, iconSurface.Get());

		SDL_SetWindowTitle(m_window, vs.title);
		SDL_ShowCursor(0);

		//SDL_GL_SetSwapInterval((vs.vsync != 0) ? -1 : 0);

		return new VkRenderer(m_window, vs);
	}

	void VkRenderer::RegisterRenderer()
	{
		Graphics::RegisterRenderer(Graphics::RENDERER_VULKAN, CreateRenderer);
	}

	VkRenderer::VkRenderer(SDL_Window *window, const Graphics::Settings &vs) :
		Renderer(window, vs.width, vs.height),
		m_rt(nullptr)
	{
		PROFILE_SCOPED()

		const char *sdl_driver = SDL_GetCurrentVideoDriver();
		Log::Info("SDL video driver used: {}", sdl_driver);

		InitVulkan();
	}

	VkRenderer::~VkRenderer()
	{
		for (size_t i = 0; i < m_swapchain.image_count; i++) {
			m_dispTable.destroySemaphore(m_finishedSemaphore[i], nullptr);
		}
		for (size_t i = 0; i < MAX_FRAMES_IN_FLIGHT; i++) {
			m_dispTable.destroySemaphore(m_availableSemaphores[i], nullptr);
			m_dispTable.destroyFence(m_inFlightFences[i], nullptr);
		}

		m_dispTable.destroyCommandPool(m_commandPool, nullptr);

		for (auto framebuffer : m_framebuffers) {
			m_dispTable.destroyFramebuffer(framebuffer, nullptr);
		}

		m_dispTable.destroyPipeline(m_graphicsPipeline, nullptr);
		m_dispTable.destroyPipelineLayout(m_pipelineLayout, nullptr);
		m_dispTable.destroyRenderPass(m_renderPass, nullptr);

		m_swapchain.destroy_image_views(m_swapchainImageViews);

		vkb::destroy_swapchain(m_swapchain);
		vkb::destroy_device(m_device);
		vkb::destroy_surface(m_instance, m_surface);
		vkb::destroy_instance(m_instance);

		SDL_DestroyWindow(m_window);

		SDL_Vulkan_UnloadLibrary();
	}

	int VkRenderer::InitVulkan()
	{
		if (0 != DeviceInitialization()) return -1;
		if (0 != CreateSwapchain()) return -1;
		if (0 != GetQueues()) return -1;
		if (0 != CreateRenderPass()) return -1;
		if (0 != CreateGraphicsPipeline()) return -1;
		if (0 != CreateFramebuffers()) return -1;
		if (0 != CreateCommandPool()) return -1;
		if (0 != CreateCommandBuffers()) return -1;
		if (0 != CreateSyncObjects()) return -1;

		return 0;
	}

	int VkRenderer::DeviceInitialization()
	{
		vkb::InstanceBuilder instance_builder;
		auto instance_ret = instance_builder.use_default_debug_messenger().request_validation_layers().build();
		if (!instance_ret) {
			Log::Info(instance_ret.error().message().c_str());
			return -1;
		}
		m_instance = instance_ret.value();

		m_instDispTable = m_instance.make_table();

		if (SDL_FALSE == SDL_Vulkan_CreateSurface(m_window, m_instance, &m_surface)) {
			Log::Error("SDL_Vulkan_CreateSurface failed");
		}

		vkb::PhysicalDeviceSelector phys_device_selector(m_instance);

		auto phys_device_ret = phys_device_selector.set_surface(m_surface).select();
		if (!phys_device_ret) {
			Log::Info(phys_device_ret.error().message().c_str());
			if (phys_device_ret.error() == vkb::PhysicalDeviceError::no_suitable_device) {
				const auto& detailed_reasons = phys_device_ret.detailed_failure_reasons();
				if (!detailed_reasons.empty()) {
					Log::Error("GPU Selection failure reasons:\n");
					for (const std::string& reason : detailed_reasons) {
						Log::Error((reason + "\n").c_str());
					}
				}
			}
			return -1;
		}
		vkb::PhysicalDevice physical_device = phys_device_ret.value();

		vkb::DeviceBuilder device_builder{ physical_device };

		auto device_ret = device_builder.build();
		if (!device_ret) {
			Log::Error(device_ret.error().message().c_str());
			return -1;
		}
		m_device = device_ret.value();

		m_dispTable = m_device.make_table();

		return 0;
	}

	int VkRenderer::CreateSwapchain()
	{
		vkb::SwapchainBuilder swapchain_builder{ m_device };
		auto swap_ret = swapchain_builder.set_old_swapchain(m_swapchain).build();
		if (!swap_ret) {
			Log::Error(fmt::format("{} {}", swap_ret.error().message(), VkErrorString(swap_ret.vk_result()) ).c_str());
			return -1;
		}
		vkb::destroy_swapchain(m_swapchain);
		m_swapchain = swap_ret.value();

		return 0;
	}

	int VkRenderer::GetQueues()
	{
		auto gq = m_device.get_queue(vkb::QueueType::graphics);
		if (!gq.has_value()) {
			Log::Error(fmt::format("failed to get graphics queue: {}", gq.error().message()).c_str());
			return -1;
		}
		m_graphicsQueue = gq.value();

		auto pq = m_device.get_queue(vkb::QueueType::present);
		if (!pq.has_value()) {
			Log::Error(fmt::format("failed to get present queue: {}", pq.error().message()).c_str());
			return -1;
		}
		m_presentQueue = pq.value();
		return 0;
	}

	int VkRenderer::CreateRenderPass()
	{
		VkAttachmentDescription color_attachment = {};
		color_attachment.format = m_swapchain.image_format;
		color_attachment.samples = VK_SAMPLE_COUNT_1_BIT;
		color_attachment.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
		color_attachment.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
		color_attachment.stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
		color_attachment.stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
		color_attachment.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
		color_attachment.finalLayout = VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;

		VkAttachmentReference color_attachment_ref = {};
		color_attachment_ref.attachment = 0;
		color_attachment_ref.layout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;

		VkSubpassDescription subpass = {};
		subpass.pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS;
		subpass.colorAttachmentCount = 1;
		subpass.pColorAttachments = &color_attachment_ref;

		VkSubpassDependency dependency = {};
		dependency.srcSubpass = VK_SUBPASS_EXTERNAL;
		dependency.dstSubpass = 0;
		dependency.srcStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
		dependency.srcAccessMask = 0;
		dependency.dstStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
		dependency.dstAccessMask = VK_ACCESS_COLOR_ATTACHMENT_READ_BIT | VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;

		VkRenderPassCreateInfo render_pass_info = {};
		render_pass_info.sType = VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO;
		render_pass_info.attachmentCount = 1;
		render_pass_info.pAttachments = &color_attachment;
		render_pass_info.subpassCount = 1;
		render_pass_info.pSubpasses = &subpass;
		render_pass_info.dependencyCount = 1;
		render_pass_info.pDependencies = &dependency;

		if (m_dispTable.createRenderPass(&render_pass_info, nullptr, &m_renderPass) != VK_SUCCESS) {
			Log::Error("failed to create render pass");
			return -1; // failed to create render pass!
		}
		return 0;
	}

	int VkRenderer::CreateGraphicsPipeline()
	{
		/*auto vert_code = readFile(std::string(EXAMPLE_SOURCE_DIRECTORY) + "/example/shaders/triangle.vert.spv");
		auto frag_code = readFile(std::string(EXAMPLE_SOURCE_DIRECTORY) + "/example/shaders/triangle.frag.spv");

		VkShaderModule vert_module = createShaderModule(init, vert_code);
		VkShaderModule frag_module = createShaderModule(init, frag_code);
		if (vert_module == VK_NULL_HANDLE || frag_module == VK_NULL_HANDLE) {
			sLog::Error("failed to create shader module");
			return -1; // failed to create shader modules
		}

		VkPipelineShaderStageCreateInfo vert_stage_info = {};
		vert_stage_info.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
		vert_stage_info.stage = VK_SHADER_STAGE_VERTEX_BIT;
		vert_stage_info.module = vert_module;
		vert_stage_info.pName = "main";

		VkPipelineShaderStageCreateInfo frag_stage_info = {};
		frag_stage_info.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
		frag_stage_info.stage = VK_SHADER_STAGE_FRAGMENT_BIT;
		frag_stage_info.module = frag_module;
		frag_stage_info.pName = "main";

		VkPipelineShaderStageCreateInfo shader_stages[] = { vert_stage_info, frag_stage_info };

		VkPipelineVertexInputStateCreateInfo vertex_input_info = {};
		vertex_input_info.sType = VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO;
		vertex_input_info.vertexBindingDescriptionCount = 0;
		vertex_input_info.vertexAttributeDescriptionCount = 0;

		VkPipelineInputAssemblyStateCreateInfo input_assembly = {};
		input_assembly.sType = VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO;
		input_assembly.topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;
		input_assembly.primitiveRestartEnable = VK_FALSE;

		VkViewport viewport = {};
		viewport.x = 0.0f;
		viewport.y = 0.0f;
		viewport.width = (float)m_swapchain.extent.width;
		viewport.height = (float)m_swapchain.extent.height;
		viewport.minDepth = 0.0f;
		viewport.maxDepth = 1.0f;

		VkRect2D scissor = {};
		scissor.offset = { 0, 0 };
		scissor.extent = m_swapchain.extent;

		VkPipelineViewportStateCreateInfo viewport_state = {};
		viewport_state.sType = VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO;
		viewport_state.viewportCount = 1;
		viewport_state.pViewports = &viewport;
		viewport_state.scissorCount = 1;
		viewport_state.pScissors = &scissor;

		VkPipelineRasterizationStateCreateInfo rasterizer = {};
		rasterizer.sType = VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO;
		rasterizer.depthClampEnable = VK_FALSE;
		rasterizer.rasterizerDiscardEnable = VK_FALSE;
		rasterizer.polygonMode = VK_POLYGON_MODE_FILL;
		rasterizer.lineWidth = 1.0f;
		rasterizer.cullMode = VK_CULL_MODE_BACK_BIT;
		rasterizer.frontFace = VK_FRONT_FACE_CLOCKWISE;
		rasterizer.depthBiasEnable = VK_FALSE;

		VkPipelineMultisampleStateCreateInfo multisampling = {};
		multisampling.sType = VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO;
		multisampling.sampleShadingEnable = VK_FALSE;
		multisampling.rasterizationSamples = VK_SAMPLE_COUNT_1_BIT;

		VkPipelineColorBlendAttachmentState colorBlendAttachment = {};
		colorBlendAttachment.colorWriteMask =
			VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT | VK_COLOR_COMPONENT_B_BIT | VK_COLOR_COMPONENT_A_BIT;
		colorBlendAttachment.blendEnable = VK_FALSE;

		VkPipelineColorBlendStateCreateInfo color_blending = {};
		color_blending.sType = VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO;
		color_blending.logicOpEnable = VK_FALSE;
		color_blending.logicOp = VK_LOGIC_OP_COPY;
		color_blending.attachmentCount = 1;
		color_blending.pAttachments = &colorBlendAttachment;
		color_blending.blendConstants[0] = 0.0f;
		color_blending.blendConstants[1] = 0.0f;
		color_blending.blendConstants[2] = 0.0f;
		color_blending.blendConstants[3] = 0.0f;

		VkPipelineLayoutCreateInfo pipeline_layout_info = {};
		pipeline_layout_info.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
		pipeline_layout_info.setLayoutCount = 0;
		pipeline_layout_info.pushConstantRangeCount = 0;

		if (m_dispTable.createPipelineLayout(&pipeline_layout_info, nullptr, &m_pipelineLayout) != VK_SUCCESS) {
			Log::Error("failed to create pipeline layout");
			return -1; // failed to create pipeline layout
		}

		std::vector<VkDynamicState> dynamic_states = { VK_DYNAMIC_STATE_VIEWPORT, VK_DYNAMIC_STATE_SCISSOR };

		VkPipelineDynamicStateCreateInfo dynamic_info = {};
		dynamic_info.sType = VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO;
		dynamic_info.dynamicStateCount = static_cast<uint32_t>(dynamic_states.size());
		dynamic_info.pDynamicStates = dynamic_states.data();

		VkGraphicsPipelineCreateInfo pipeline_info = {};
		pipeline_info.sType = VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO;
		pipeline_info.stageCount = 2;
		pipeline_info.pStages = shader_stages;
		pipeline_info.pVertexInputState = &vertex_input_info;
		pipeline_info.pInputAssemblyState = &input_assembly;
		pipeline_info.pViewportState = &viewport_state;
		pipeline_info.pRasterizationState = &rasterizer;
		pipeline_info.pMultisampleState = &multisampling;
		pipeline_info.pColorBlendState = &color_blending;
		pipeline_info.pDynamicState = &dynamic_info;
		pipeline_info.layout = m_pipelineLayout;
		pipeline_info.renderPass = m_renderPass;
		pipeline_info.subpass = 0;
		pipeline_info.basePipelineHandle = VK_NULL_HANDLE;

		if (m_dispTable.createGraphicsPipelines(VK_NULL_HANDLE, 1, &pipeline_info, nullptr, &m_graphicsPipeline) != VK_SUCCESS) {
			Log::Error("failed to create pipline");
			return -1; // failed to create graphics pipeline
		}

		m_dispTable.destroyShaderModule(frag_module, nullptr);
		m_dispTable.destroyShaderModule(vert_module, nullptr);*/
		return 0;
	}

	int VkRenderer::CreateFramebuffers()
	{
		m_swapchainImages = m_swapchain.get_images().value();
		m_swapchainImageViews = m_swapchain.get_image_views().value();

		m_framebuffers.resize(m_swapchainImageViews.size());

		for (size_t i = 0; i < m_swapchainImageViews.size(); i++) {
			VkImageView attachments[] = { m_swapchainImageViews[i] };

			VkFramebufferCreateInfo framebuffer_info = {};
			framebuffer_info.sType = VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO;
			framebuffer_info.renderPass = m_renderPass;
			framebuffer_info.attachmentCount = 1;
			framebuffer_info.pAttachments = attachments;
			framebuffer_info.width = m_swapchain.extent.width;
			framebuffer_info.height = m_swapchain.extent.height;
			framebuffer_info.layers = 1;

			if (m_dispTable.createFramebuffer(&framebuffer_info, nullptr, &m_framebuffers[i]) != VK_SUCCESS) {
				return -1; // failed to create framebuffer
			}
		}
		return 0;
	}

	int VkRenderer::CreateCommandPool()
	{
		VkCommandPoolCreateInfo pool_info = {};
		pool_info.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
		pool_info.queueFamilyIndex = m_device.get_queue_index(vkb::QueueType::graphics).value();

		if (m_dispTable.createCommandPool(&pool_info, nullptr, &m_commandPool) != VK_SUCCESS) {
			Log::Error("failed to create command pool");
			return -1; // failed to create command pool
		}
		return 0;
	}

	int VkRenderer::CreateCommandBuffers()
	{
		m_commandBuffers.resize(m_framebuffers.size());

		VkCommandBufferAllocateInfo allocInfo = {};
		allocInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
		allocInfo.commandPool = m_commandPool;
		allocInfo.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
		allocInfo.commandBufferCount = (uint32_t)m_commandBuffers.size();

		if (m_dispTable.allocateCommandBuffers(&allocInfo, m_commandBuffers.data()) != VK_SUCCESS) {
			return -1; // failed to allocate command buffers;
		}

		for (size_t i = 0; i < m_commandBuffers.size(); i++) {
			VkCommandBufferBeginInfo begin_info = {};
			begin_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;

			if (m_dispTable.beginCommandBuffer(m_commandBuffers[i], &begin_info) != VK_SUCCESS) {
				return -1; // failed to begin recording command buffer
			}

			VkRenderPassBeginInfo render_pass_info = {};
			render_pass_info.sType = VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO;
			render_pass_info.renderPass = m_renderPass;
			render_pass_info.framebuffer = m_framebuffers[i];
			render_pass_info.renderArea.offset = { 0, 0 };
			render_pass_info.renderArea.extent = m_swapchain.extent;
			VkClearValue clearColor{ { { 0.0f, 0.0f, 0.0f, 1.0f } } };
			render_pass_info.clearValueCount = 1;
			render_pass_info.pClearValues = &clearColor;

			VkViewport viewport = {};
			viewport.x = 0.0f;
			viewport.y = 0.0f;
			viewport.width = (float)m_swapchain.extent.width;
			viewport.height = (float)m_swapchain.extent.height;
			viewport.minDepth = 0.0f;
			viewport.maxDepth = 1.0f;

			VkRect2D scissor = {};
			scissor.offset = { 0, 0 };
			scissor.extent = m_swapchain.extent;

			m_dispTable.cmdSetViewport(m_commandBuffers[i], 0, 1, &viewport);
			m_dispTable.cmdSetScissor(m_commandBuffers[i], 0, 1, &scissor);

			m_dispTable.cmdBeginRenderPass(m_commandBuffers[i], &render_pass_info, VK_SUBPASS_CONTENTS_INLINE);

			m_dispTable.cmdBindPipeline(m_commandBuffers[i], VK_PIPELINE_BIND_POINT_GRAPHICS, m_graphicsPipeline);

			m_dispTable.cmdDraw(m_commandBuffers[i], 3, 1, 0, 0);

			m_dispTable.cmdEndRenderPass(m_commandBuffers[i]);

			if (m_dispTable.endCommandBuffer(m_commandBuffers[i]) != VK_SUCCESS) {
				Log::Error("failed to record command buffer");
				return -1; // failed to record command buffer!
			}
		}
		return 0;
	}

	int VkRenderer::CreateSyncObjects()
	{
		m_availableSemaphores.resize(MAX_FRAMES_IN_FLIGHT);
		m_finishedSemaphore.resize(m_swapchain.image_count);
		m_inFlightFences.resize(MAX_FRAMES_IN_FLIGHT);
		m_imageInFlight.resize(m_swapchain.image_count, VK_NULL_HANDLE);

		VkSemaphoreCreateInfo semaphore_info = {};
		semaphore_info.sType = VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO;

		VkFenceCreateInfo fence_info = {};
		fence_info.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
		fence_info.flags = VK_FENCE_CREATE_SIGNALED_BIT;

		for (size_t i = 0; i < m_swapchain.image_count; i++) {
			if (m_dispTable.createSemaphore(&semaphore_info, nullptr, &m_finishedSemaphore[i]) != VK_SUCCESS) {
				Log::Error("failed to create sync objects");
				return -1; // failed to create synchronization objects for a frame
			}
		}

		for (size_t i = 0; i < MAX_FRAMES_IN_FLIGHT; i++) {
			if (m_dispTable.createSemaphore(&semaphore_info, nullptr, &m_availableSemaphores[i]) != VK_SUCCESS ||
				m_dispTable.createFence(&fence_info, nullptr, &m_inFlightFences[i]) != VK_SUCCESS) {
				Log::Error("failed to create sync objects");
				return -1; // failed to create synchronization objects for a frame
			}
		}
		return 0;
	}

} // namespace Graphics
