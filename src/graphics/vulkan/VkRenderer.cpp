// Copyright © 2008-2024 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "VkRenderer.h"

#include "MathUtil.h"
#include "RefCounted.h"
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

#include <vulkan/vulkan.hpp>

#include "vk-bootstrap/VkBootstrap.h"

#include <SDL.h>
#include <SDL_vulkan.h>

namespace Graphics {

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
		SDL_GLContext glContext = nullptr;

		bool ok = LoadLibraryAndCreateWindow(name.c_str(), vs, m_window);
		if (!ok) {
			Error("Failed to set video mode: %s", SDL_GetError());
			return nullptr;
		}

		SDLSurfacePtr surface = LoadSurfaceFromFile(vs.iconFile);
		if (surface)
			SDL_SetWindowIcon(m_window, surface.Get());

		SDL_SetWindowTitle(m_window, vs.title);
		SDL_ShowCursor(0);

		//SDL_GL_SetSwapInterval((vs.vsync != 0) ? -1 : 0);

		return new VkRenderer(m_window, vs);
	}

	void VkRenderer::RegisterRenderer()
	{
		Graphics::RegisterRenderer(Graphics::RENDERER_VULKAN, CreateRenderer);
	}

	VkRenderer::VkRenderer(SDL_Window *m_window, const Graphics::Settings &vs) :
		Renderer(m_window, vs.width, vs.height),
		m_identity(matrix4x4f::Identity()),
		m_rt(nullptr)
	{
		PROFILE_SCOPED()

		const char *sdl_driver = SDL_GetCurrentVideoDriver();
		Log::Info("SDL video driver used: {}", sdl_driver);

		InitVulkan();
	}

	VkRenderer::~VkRenderer()
	{
		vkDestroyDevice(m_device, nullptr);
		vkDestroyInstance(m_vkInst, nullptr);
		SDL_Vulkan_UnloadLibrary();
	}

	void VkRenderer::InitVulkan()
	{
		VkResult result;
		uint32_t extensionCount;
		const char **extensionNames = 0;
		SDL_Vulkan_GetInstanceExtensions(m_window, &extensionCount, nullptr);
		extensionNames = new const char *[extensionCount];
		SDL_Vulkan_GetInstanceExtensions(m_window, &extensionCount, extensionNames);
		const VkInstanceCreateInfo instInfo = {
			VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO, // sType
			nullptr,								// pNext
			0,										// flags
			nullptr,								// pApplicationInfo
			0,										// enabledLayerCount
			nullptr,								// ppEnabledLayerNames
			extensionCount,							// enabledExtensionCount
			extensionNames,							// ppEnabledExtensionNames
		};
		result = vkCreateInstance(&instInfo, nullptr, &m_vkInst);
		if (result != VK_SUCCESS) {
			Error("Unable to create Vulkan instance : %s", VkErrorString(result).c_str());
			return;
		}

		uint32_t physicalDeviceCount;
		result = vkEnumeratePhysicalDevices(m_vkInst, &physicalDeviceCount, nullptr);
		if (result != VK_SUCCESS || physicalDeviceCount == 0) {
			Error("No device with Vulkan support found : %s", VkErrorString(result).c_str());
			return;
		}

		std::vector<VkPhysicalDevice> physicalDevices(physicalDeviceCount);
		result = vkEnumeratePhysicalDevices(m_vkInst, &physicalDeviceCount, physicalDevices.data());
		if (result != VK_SUCCESS) {
			Error("Could not enumerate physical devices : %s", VkErrorString(result).c_str());
			return;
		}

		// pick the best device... or just the first one for now
		VkPhysicalDevice physicalDevice = physicalDevices[0];

		uint32_t queueFamilyCount;
		vkGetPhysicalDeviceQueueFamilyProperties(physicalDevice, &queueFamilyCount, nullptr);
		std::vector<VkQueueFamilyProperties> queueFamilies(queueFamilyCount);
		vkGetPhysicalDeviceQueueFamilyProperties(physicalDevice, &queueFamilyCount, queueFamilies.data());

		VkSurfaceKHR surface;
		SDL_Vulkan_CreateSurface(m_window, m_vkInst, &surface);

		uint32_t graphicsQueueIndex = UINT32_MAX;
		uint32_t presentQueueIndex = UINT32_MAX;
		VkBool32 support;
		uint32_t i = 0;
		for (const VkQueueFamilyProperties &queueFamily : queueFamilies) {
			if (graphicsQueueIndex == UINT32_MAX && queueFamily.queueCount > 0 && queueFamily.queueFlags & VK_QUEUE_GRAPHICS_BIT)
				graphicsQueueIndex = i;
			if (presentQueueIndex == UINT32_MAX) {
				vkGetPhysicalDeviceSurfaceSupportKHR(physicalDevice, i, surface, &support);
				if (support)
					presentQueueIndex = i;
			}
			++i;
		}

		float queuePriority = 1.0f;
		VkDeviceQueueCreateInfo queueInfo = {
			VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO, // sType
			nullptr,									// pNext
			0,											// flags
			graphicsQueueIndex,							// graphicsQueueIndex
			1,											// queueCount
			&queuePriority,								// pQueuePriorities
		};

		VkPhysicalDeviceFeatures deviceFeatures = {};
		const char *deviceExtensionNames[] = { VK_KHR_SWAPCHAIN_EXTENSION_NAME };
		VkDeviceCreateInfo createInfo = {
			VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO, // sType
			nullptr,							  // pNext
			0,									  // flags
			1,									  // queueCreateInfoCount
			&queueInfo,							  // pQueueCreateInfos
			0,									  // enabledLayerCount
			nullptr,							  // ppEnabledLayerNames
			1,									  // enabledExtensionCount
			deviceExtensionNames,				  // ppEnabledExtensionNames
			&deviceFeatures,					  // pEnabledFeatures
		};
		result = vkCreateDevice(physicalDevice, &createInfo, nullptr, &m_device);
		if (result != VK_SUCCESS) {
			Error("vkCreateDevice failed to create logical device : %s", VkErrorString(result).c_str());
			return;
		}

		vkGetDeviceQueue(m_device, graphicsQueueIndex, 0, &m_graphicsQueue);

		vkGetDeviceQueue(m_device, presentQueueIndex, 0, &m_presentQueue);
	}

} // namespace Graphics
