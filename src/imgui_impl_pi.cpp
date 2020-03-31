// ImGui SDL2 binding with OpenGL3
// In this binding, ImTextureID is used to store an OpenGL 'GLuint' texture identifier. Read the FAQ about ImTextureID in imgui.cpp.

// You can copy and use unmodified imgui_impl_* files in your project. See main.cpp for an example of using this.
// If you use this binding you'll need to call 4 functions: ImGui_ImplXXXX_Init(), ImGui_ImplXXXX_NewFrame(), ImGui::Render() and ImGui_ImplXXXX_Shutdown().
// If you are new to ImGui, see examples/README.txt and documentation at the top of imgui.cpp.
// https://github.com/ocornut/imgui


#include "imgui_impl_pi.h"
#include "Pi.h"

#include "graphics/Material.h"
#include "graphics/Types.h"
#include "graphics/VertexBuffer.h"

// SDL,GL3W
#include <SDL.h>
#include <SDL_syswm.h>
#define GL_FUNC_ADD 0x8006

// Data
static double       g_Time = 0.0f;
static bool         g_MousePressed[3] = { false, false, false };
static float        g_MouseWheel = 0.0f;
static uint32_t     g_FontTexture = 0;

static Graphics::RenderState *pRS = nullptr;
static RefCountedPtr<Graphics::Texture> fontTexture;
static RefCountedPtr<Graphics::Material> rMaterial;
static RefCountedPtr<Graphics::VertexBuffer> vertexBuffer;
static RefCountedPtr<Graphics::IndexBuffer> indexBuffer;

// Pioneer ImGui Render function.
// (this used to be set in io.RenderDrawListsFn and called by ImGui::Render(), but you can now call this directly from your main loop)

// This is the main rendering function that you have to implement and provide to ImGui (via setting up 'RenderDrawListsFn' in the ImGuiIO structure)
// If text or lines are blurry when integrating ImGui in your engine:
// - in your Render function, try translating your projection matrix by (0.5f,0.5f) or (0.375f,0.375f)
void ImGui_ImplPi_RenderDrawData(ImDrawData *draw_data)
{
	// Avoid rendering when minimized, scale coordinates for retina displays (screen coordinates != framebuffer coordinates)
	ImGuiIO& io = ImGui::GetIO();
	int fb_width = (int)(io.DisplaySize.x * io.DisplayFramebufferScale.x);
	int fb_height = (int)(io.DisplaySize.y * io.DisplayFramebufferScale.y);
	if ( fb_width == 0 || fb_height == 0 )
		return;
	draw_data->ScaleClipRects(io.DisplayFramebufferScale);

	//glEnable(GL_BLEND);
	//glBlendEquation(GL_FUNC_ADD);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glDisable(GL_CULL_FACE);
	//glDisable(GL_DEPTH_TEST);
	//glEnable(GL_SCISSOR_TEST);
	//glActiveTexture(GL_TEXTURE0);

	// Setup orthographic projection matrix
	const float ortho_projection[4][4] =
	{
		{ 2.0f / io.DisplaySize.x, 0.0f,                   0.0f, 0.0f },
		{ 0.0f,                  2.0f / -io.DisplaySize.y, 0.0f, 0.0f },
		{ 0.0f,                  0.0f,                  -1.0f, 0.0f },
		{-1.0f,                  1.0f,                   0.0f, 1.0f },
	};
	Pi::renderer->SetViewport(0, 0, fb_width, fb_height);
	Pi::renderer->SetTransform(matrix4x4f::Identity());
	Pi::renderer->SetProjection(matrix4x4f(&ortho_projection[0][0]));

	//ImDrawVert *vtxptr = vertexBuffer->Map<ImDrawVert>(Graphics::BUFFER_MAP_WRITE);
	//assert(vertexBuffer->GetDesc().stride == sizeof(ImDrawVert));

	//glUseProgram(g_ShaderHandle);
	//glUniform1i(g_AttribLocationTex, 0);
	//glUniformMatrix4fv(g_AttribLocationProjMtx, 1, GL_FALSE, &ortho_projection[0][0]);
	//glBindVertexArray(g_VaoHandle);

	for ( int n = 0; n < draw_data->CmdListsCount; n++ )
	{
		const ImDrawList* cmd_list = draw_data->CmdLists[n];
		const ImDrawIdx* idx_buffer_offset = 0;

		vertexBuffer->BufferData(cmd_list->VtxBuffer.Size * sizeof(ImDrawVert), cmd_list->VtxBuffer.Data);
		vertexBuffer->SetVertexCount(cmd_list->VtxBuffer.Size);

		indexBuffer->BufferData(cmd_list->IdxBuffer.Size * sizeof(ImDrawIdx), cmd_list->IdxBuffer.Data);
		indexBuffer->SetIndexCount(cmd_list->IdxBuffer.Size);

		for ( int cmd_i = 0; cmd_i < cmd_list->CmdBuffer.Size; cmd_i++ )
		{
			const ImDrawCmd* pcmd = &cmd_list->CmdBuffer[cmd_i];
			if ( pcmd->UserCallback )
			{
				pcmd->UserCallback(cmd_list, pcmd);
			}
			else
			{
				//glBindTexture(GL_TEXTURE_2D, (GLuint)(intptr_t)pcmd->TextureId);
				//glScissor((int)pcmd->ClipRect.x, (int)(fb_height - pcmd->ClipRect.w), (int)(pcmd->ClipRect.z - pcmd->ClipRect.x), (int)(pcmd->ClipRect.w - pcmd->ClipRect.y));
				Pi::renderer->SetScissor(true,vector2f((int)pcmd->ClipRect.x, (int)(fb_height - pcmd->ClipRect.w)), vector2f((int)(pcmd->ClipRect.z - pcmd->ClipRect.x), (int)(pcmd->ClipRect.w - pcmd->ClipRect.y)));
				//glDrawElements(GL_TRIANGLES, (GLsizei)pcmd->ElemCount, sizeof(ImDrawIdx) == 2 ? GL_UNSIGNED_SHORT : GL_UNSIGNED_INT, idx_buffer_offset);
				Pi::renderer->DrawBufferIndexed(vertexBuffer.Get(), indexBuffer.Get(), pRS, rMaterial.Get());
			}
			idx_buffer_offset += pcmd->ElemCount;
		}
	}

	// Restore modified state
	Pi::renderer->SetScissor(false);
}

static const char* ImGui_ImplPi_GetClipboardText(void *userData)
{
	return SDL_GetClipboardText();
}

static void ImGui_ImplPi_SetClipboardText(void *userData, const char *text)
{
	SDL_SetClipboardText(text);
}

bool ImGui_ImplPi_ProcessEvent(SDL_Event* event)
{
	ImGuiIO& io = ImGui::GetIO();
	switch ( event->type )
	{
	case SDL_MOUSEWHEEL:
	{
		if ( event->wheel.y > 0 )
			g_MouseWheel = 1;
		if ( event->wheel.y < 0 )
			g_MouseWheel = -1;
		return true;
	}
	case SDL_MOUSEBUTTONDOWN:
	{
		if ( event->button.button == SDL_BUTTON_LEFT ) g_MousePressed[0] = true;
		if ( event->button.button == SDL_BUTTON_RIGHT ) g_MousePressed[1] = true;
		if ( event->button.button == SDL_BUTTON_MIDDLE ) g_MousePressed[2] = true;
		return true;
	}
	case SDL_TEXTINPUT:
	{
		io.AddInputCharactersUTF8(event->text.text);
		return true;
	}
	case SDL_KEYDOWN:
	case SDL_KEYUP:
	{
		int key = event->key.keysym.sym;
		if ( key & SDLK_SCANCODE_MASK )
			key = (key & ~SDLK_SCANCODE_MASK) | 0x100;
		io.KeysDown[key] = (event->type == SDL_KEYDOWN);
		io.KeyShift = ((SDL_GetModState() & KMOD_SHIFT) != 0);
		io.KeyCtrl = ((SDL_GetModState() & KMOD_CTRL) != 0);
		io.KeyAlt = ((SDL_GetModState() & KMOD_ALT) != 0);
		io.KeySuper = ((SDL_GetModState() & KMOD_GUI) != 0);
		return true;
	}
	}
	return false;
}

void ImGui_ImplPi_CreateFontsTexture()
{
	// Build texture atlas
	ImGuiIO& io = ImGui::GetIO();
	unsigned char* pixels;
	int width, height;
	io.Fonts->GetTexDataAsRGBA32(&pixels, &width, &height);   // Load as RGBA 32-bits for OpenGL3 demo because it is more likely to be compatible with user's existing shader.


	const vector2f texSize(1.0f, 1.0f);
	const vector3f dataSize(width, height, 0.0f);
	Graphics::TextureDescriptor texDesc(
			Graphics::TEXTURE_RGBA_8888,
			dataSize, texSize, Graphics::LINEAR_CLAMP,
			false, false, false, 0, Graphics::TEXTURE_2D);
	fontTexture.Reset(Pi::renderer->CreateTexture(texDesc));
	// Upload texture to graphics system
	//GLint last_texture;
	//glGetIntegerv(GL_TEXTURE_BINDING_2D, &last_texture);
	//glGenTextures(1, &g_FontTexture);
	//glBindTexture(GL_TEXTURE_2D, g_FontTexture);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	//glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
	g_FontTexture = fontTexture->GetTextureID();

	// Store our identifier
	io.Fonts->TexID = (void *)(intptr_t)g_FontTexture;
}

bool ImGui_ImplPi_CreateDeviceObjects()
{
	if (vertexBuffer != nullptr)
		return true;

	// Setup render state: alpha-blending enabled, no face culling, no depth testing, scissor enabled
	Graphics::RenderStateDesc rsDesc;
	rsDesc.blendMode = Graphics::BLEND_ADDITIVE;
	rsDesc.cullMode = Graphics::CULL_NONE;
	rsDesc.depthTest = false;
	pRS = Pi::renderer->CreateRenderState(rsDesc);

//	glGenBuffers(1, &g_VboHandle);
//	glGenBuffers(1, &g_ElementsHandle);
//
//	glGenVertexArrays(1, &g_VaoHandle);
//	glBindVertexArray(g_VaoHandle);
//	glBindBuffer(GL_ARRAY_BUFFER, g_VboHandle);
//	glEnableVertexAttribArray(g_AttribLocationPosition);
//	glEnableVertexAttribArray(g_AttribLocationUV);
//	glEnableVertexAttribArray(g_AttribLocationColor);
//
//#define OFFSETOF(TYPE, ELEMENT) ((size_t)&(((TYPE *)0)->ELEMENT))
//	glVertexAttribPointer(g_AttribLocationPosition, 2, GL_FLOAT, GL_FALSE, sizeof(ImDrawVert), (GLvoid*)OFFSETOF(ImDrawVert, pos));
//	glVertexAttribPointer(g_AttribLocationUV, 2, GL_FLOAT, GL_FALSE, sizeof(ImDrawVert), (GLvoid*)OFFSETOF(ImDrawVert, uv));
//	glVertexAttribPointer(g_AttribLocationColor, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(ImDrawVert), (GLvoid*)OFFSETOF(ImDrawVert, col));
//#undef OFFSETOF
	Graphics::VertexBufferDesc vbDesc;
	vbDesc.attrib[0].semantic = Graphics::ATTRIB_POSITION;
	vbDesc.attrib[0].format   = Graphics::ATTRIB_FORMAT_FLOAT2;	// NB: this might go really wrong
	vbDesc.attrib[1].semantic = Graphics::ATTRIB_UV0;
	vbDesc.attrib[1].format   = Graphics::ATTRIB_FORMAT_FLOAT2;
	vbDesc.attrib[2].semantic = Graphics::ATTRIB_DIFFUSE;
	vbDesc.attrib[2].format   = Graphics::ATTRIB_FORMAT_UBYTE4;
	vbDesc.numVertices = 65535;//box->GetNumVerts();	// this is baaaad
	vbDesc.usage = Graphics::BUFFER_USAGE_DYNAMIC; // need one for GL_STREAM_DRAW type?
	vertexBuffer.Reset(Pi::renderer->CreateVertexBuffer(vbDesc));

	//create buffer & copy
	indexBuffer.Reset(Pi::renderer->CreateIndexBuffer(65535, Graphics::BUFFER_USAGE_DYNAMIC));


	Graphics::MaterialDescriptor desc;
	desc.effect = Graphics::EFFECT_IMGUI;
	desc.textures = 1;
	rMaterial.Reset(Pi::renderer->CreateMaterial(desc));
	rMaterial->texture0 = nullptr;

	ImGui_ImplPi_CreateFontsTexture();
	rMaterial->texture0 = fontTexture.Get();

	return true;
}

void    ImGui_ImplPi_InvalidateDeviceObjects()
{
	vertexBuffer.Reset();
	rMaterial.Reset();

	if ( g_FontTexture )
	{
		fontTexture.Reset();
		ImGui::GetIO().Fonts->TexID = 0;
		g_FontTexture = 0;
	}
}

bool    ImGui_ImplPi_Init(SDL_Window* window)
{
	ImGuiIO& io = ImGui::GetIO();
	io.KeyMap[ImGuiKey_Tab] = SDLK_TAB;                     // Keyboard mapping. ImGui will use those indices to peek into the io.KeyDown[] array.
	io.KeyMap[ImGuiKey_LeftArrow] = SDL_SCANCODE_LEFT;
	io.KeyMap[ImGuiKey_RightArrow] = SDL_SCANCODE_RIGHT;
	io.KeyMap[ImGuiKey_UpArrow] = SDL_SCANCODE_UP;
	io.KeyMap[ImGuiKey_DownArrow] = SDL_SCANCODE_DOWN;
	io.KeyMap[ImGuiKey_PageUp] = SDL_SCANCODE_PAGEUP;
	io.KeyMap[ImGuiKey_PageDown] = SDL_SCANCODE_PAGEDOWN;
	io.KeyMap[ImGuiKey_Home] = SDL_SCANCODE_HOME;
	io.KeyMap[ImGuiKey_End] = SDL_SCANCODE_END;
	io.KeyMap[ImGuiKey_Delete] = SDLK_DELETE;
	io.KeyMap[ImGuiKey_Backspace] = SDLK_BACKSPACE;
	io.KeyMap[ImGuiKey_Enter] = SDLK_RETURN;
	io.KeyMap[ImGuiKey_Escape] = SDLK_ESCAPE;
	io.KeyMap[ImGuiKey_A] = SDLK_a;
	io.KeyMap[ImGuiKey_C] = SDLK_c;
	io.KeyMap[ImGuiKey_V] = SDLK_v;
	io.KeyMap[ImGuiKey_X] = SDLK_x;
	io.KeyMap[ImGuiKey_Y] = SDLK_y;
	io.KeyMap[ImGuiKey_Z] = SDLK_z;

	io.SetClipboardTextFn = ImGui_ImplPi_SetClipboardText;
	io.GetClipboardTextFn = ImGui_ImplPi_GetClipboardText;

#ifdef _WIN32
	SDL_SysWMinfo wmInfo;
	SDL_VERSION(&wmInfo.version);
	SDL_GetWindowWMInfo(window, &wmInfo);
	io.ImeWindowHandle = wmInfo.info.win.window;
#else
	(void)window;
#endif

	return true;
}

void ImGui_ImplPi_Shutdown()
{
	ImGui_ImplPi_InvalidateDeviceObjects();
}

void ImGui_ImplPi_NewFrame(SDL_Window* window)
{
	if ( !g_FontTexture )
		ImGui_ImplPi_CreateDeviceObjects();

	ImGuiIO& io = ImGui::GetIO();

	// Setup display size (every frame to accommodate for window resizing)
	int w, h;
	int display_w, display_h;
	SDL_GetWindowSize(window, &w, &h);
	SDL_GL_GetDrawableSize(window, &display_w, &display_h);
	io.DisplaySize = ImVec2((float)w, (float)h);
	io.DisplayFramebufferScale = ImVec2(w > 0 ? ((float)display_w / w) : 0, h > 0 ? ((float)display_h / h) : 0);

	// Setup time step
	Uint32	time = SDL_GetTicks();
	double current_time = time / 1000.0;
	io.DeltaTime = g_Time > 0.0 ? (float)(current_time - g_Time) : (float)(1.0f / 60.0f);

	if ( io.DeltaTime <= 0.0 ) {
		// Delta-T should *never* be negative, but there seem to be bugs in SDL_GetTicks
		// (or the underlying platform implementation) which can lead to this happening.
		// In that case, just assume 60 FPS.
		io.DeltaTime = 1.0 / 60.0;
	}

	g_Time = current_time;

	// Setup inputs
	// (we already got mouse wheel, keyboard keys & characters from SDL_PollEvent())
	int mx, my;
	Uint32 mouseMask = SDL_GetMouseState(&mx, &my);
	if ( SDL_GetWindowFlags(window) & SDL_WINDOW_MOUSE_FOCUS )
		io.MousePos = ImVec2((float)mx, (float)my);   // Mouse position, in pixels (set to -1,-1 if no mouse / on another screen, etc.)
	else
		io.MousePos = ImVec2(-1, -1);

	io.MouseDown[0] = g_MousePressed[0] || (mouseMask & SDL_BUTTON(SDL_BUTTON_LEFT)) != 0;		// If a mouse press event came, always pass it as "mouse held this frame", so we don't miss click-release events that are shorter than 1 frame.
	io.MouseDown[1] = g_MousePressed[1] || (mouseMask & SDL_BUTTON(SDL_BUTTON_RIGHT)) != 0;
	io.MouseDown[2] = g_MousePressed[2] || (mouseMask & SDL_BUTTON(SDL_BUTTON_MIDDLE)) != 0;
	g_MousePressed[0] = g_MousePressed[1] = g_MousePressed[2] = false;

	io.MouseWheel = g_MouseWheel;
	g_MouseWheel = 0.0f;

	// Hide OS mouse cursor if ImGui is drawing it
	SDL_ShowCursor(io.MouseDrawCursor ? 0 : 1);
}
