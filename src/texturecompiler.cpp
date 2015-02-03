// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "libs.h"
#include "utils.h"
#include <cstdio>
#include <cstdlib>

#include "scenegraph/SceneGraph.h"

#include "FileSystem.h"
#include "GameConfig.h"
#include "graphics/dummy/RendererDummy.h"
#include "graphics/opengl/RendererGL.h"
#include "graphics/Graphics.h"
#include "graphics/Light.h"
#include "graphics/Renderer.h"
#include "graphics/Texture.h"
#include "graphics/TextureBuilder.h"
#include "graphics/Drawables.h"
#include "graphics/VertexArray.h"
#include "scenegraph/DumpVisitor.h"
#include "scenegraph/FindNodeVisitor.h"
#include "scenegraph/BinaryConverter.h"
#include "OS.h"
#include "StringF.h"
#include "ModManager.h"
#include <sstream>

std::unique_ptr<GameConfig> s_config;
std::unique_ptr<Graphics::Renderer> s_renderer;

static const std::string s_dummyPath("");

void SetupRenderer()
{
	s_config.reset(new GameConfig);

	OS::RedirectStdio();

	//init components
	FileSystem::userFiles.MakeDirectory(""); // ensure the config directory exists
	if (SDL_Init(SDL_INIT_VIDEO) < 0)
		Error("SDL initialization failed: %s\n", SDL_GetError());

	ModManager::Init();

	Graphics::RendererDummy::RegisterRenderer();
	Graphics::RendererOGL::RegisterRenderer();

	//video
	Graphics::Settings videoSettings = {};
	videoSettings.rendererType = Graphics::RENDERER_OPENGL;
	videoSettings.width = s_config->Int("ScrWidth");
	videoSettings.height = s_config->Int("ScrHeight");
	videoSettings.fullscreen = false;
	videoSettings.hidden = true;
	videoSettings.requestedSamples = s_config->Int("AntiAliasingMode");
	videoSettings.vsync = false;
	videoSettings.useTextureCompression = true;
	videoSettings.iconFile = OS::GetIconFilename();
	videoSettings.title = "Texture Compiler";
	s_renderer.reset(Graphics::Init(videoSettings));
}

/*vector2f encode(const vector3f& n)
{
	const vector2f nxy(n.x, n.y);
	vector2f enc = nxy.Normalized() * (sqrt(-n.z*0.5 + 0.5));
	enc = (enc * 0.5f) + vector2f(0.5f);
	return enc;
}

vector3f decode(const vector4f& enc)
{
	const vector4f enc2 = enc * vector4f(2, 2, 0, 0);
	vector4f nn = enc2 + vector4f(-1, -1, 1, -1);
	const float l = nn.xyz().Dot(-vector3f(nn.x, nn.y, nn.w));
	nn.z = l;
	nn.xy() *= sqrt(l);
	return (nn.xyz() * 2.0f) + vector3f(0, 0, -1);

	//vector4f nn = enc*vector4f(2,2,0,0) + vector4f(-1,-1,1,-1);
	//const float l = dot(nn.xyz,-nn.xyw);
	//nn.z = l;
	//nn.xy *= sqrt(l);
	//return nn.xyz * 2 + vector3f(0,0,-1);
}*/

vector2f encode(const vector3f& n)
{
	const float scale = 1.7777f;
	vector2f enc = vector2f(n.x, n.y) / (n.z + 1.0f);
	enc /= scale;
	enc = (enc * 0.5f) + 0.5f;
	return enc;
}

vector3f decode(const vector4f& enc)
{
	const float scale = 1.7777f;
	const vector3f nn = enc.xyz() * vector3f(2.0f * scale, 2.0f * scale, 0.0f) + vector3f(-scale, -scale, 1.0f);
	const float g = 2.0f / nn.Dot(nn);
	vector3f n;
	n.x = g * nn.x;
	n.y = g * nn.y;
	n.z = g - 1.0f;
	return n;
}

void RunCompiler(const std::string &modelName, const std::string &filepath)
{
	Profiler::Timer timer;
	timer.Start();
	Output("\n---\nStarting compiler for (%s)\n", modelName.c_str());

	//load the image, encode it, then save the encoded version
	std::unique_ptr<vector3f> imageData;
	int Width = 0;
	int Height = 0;
	try {
		SDLSurfacePtr im = LoadSurfaceFromFile(modelName);
		if (im) {
			// now that we have our raw image loaded
			// allocate the space for our processed representation
			imageData.reset(new vector3f[(im->w * im->h)]);
			vector3f* pimg = imageData.get();

			// lock the image once so we can read from it
			SDL_LockSurface(im.Get());

			// setup our map dimensions for later
			Width = im->w;
			Height = im->w;

			// copy every pixel value from the red channel (image is greyscale, channel is irrelevant)
			for (int x = 0; x<im->w; x++) {
				for (int y = 0; y<im->h; y++) {
					const unsigned char v0 = static_cast<unsigned char*>(im->pixels)[((x * 3) + 0) + (y*im->pitch)];
					const unsigned char v1 = static_cast<unsigned char*>(im->pixels)[((x * 3) + 1) + (y*im->pitch)];
					const unsigned char v2 = static_cast<unsigned char*>(im->pixels)[((x * 3) + 2) + (y*im->pitch)];
					pimg[x + y*Width] = vector3f(v0 / 255.0f, v1 / 255.0f, v2 / 255.0f);
				}
			}

			// unlock the surface and then release it
			SDL_UnlockSurface(im.Get());
			/*if (im) {
				SDL_FreeSurface(im.Get());
			}*/
		} else {
			Output("Failed to load image.\n");
		}
	} catch (...) {
		//minimal error handling, this is not expected to happen since we got this far.
		return;
	}

	try {
		const int xmax = Width;
		const int ymax = Height;
		for (int x = 0; x < xmax; x++) {
			for (int y = 0; y < ymax; y++) {
				const vector3f val = imageData.get()[(x + (y * xmax))].Normalized();
				const vector2f enc = encode(val);
				const vector3f dec = decode(vector4f(enc.x, enc.y, 0.0f, 0.0f));
				const vector3f decN = dec.Normalized();
				assert(dec == val);
			}
		}
		//const std::string DataPath = FileSystem::NormalisePath(filepath.substr(0, filepath.size()-6));
		//SceneGraph::BinaryConverter bc(s_renderer.get());
		//bc.Save(modelName, DataPath, model.get(), false);
	} catch (const CouldNotOpenFileException&) {
	} catch (const CouldNotWriteToFileException&) {
	}

	timer.Stop();
	Output("Compiling \"%s\" took: %lf\n", modelName.c_str(), timer.millicycles());
}


enum RunMode {
	MODE_MODELCOMPILER=0,
	MODE_MODELBATCHEXPORT,
	MODE_VERSION,
	MODE_USAGE,
	MODE_USAGE_ERROR
};

int main(int argc, char** argv)
{
#ifdef PIONEER_PROFILER
	Profiler::detect( argc, argv );
#endif

	RunMode mode = MODE_USAGE_ERROR;

	if (argc > 1) {
		const char switchchar = argv[1][0];
		const std::string modeopt(std::string(argv[1]).substr(1));

		if (!(switchchar == '-' || switchchar == '/')) {
			mode = MODE_USAGE_ERROR;
		} else if (modeopt == "compile" || modeopt == "c") {
			mode = MODE_MODELCOMPILER;
		} else if (modeopt == "version" || modeopt == "v") {
			mode = MODE_VERSION;
		} else if (modeopt == "help" || modeopt == "h" || modeopt == "?") {
			mode = MODE_USAGE;
		}
	}
	
	// Init here since we'll need it for both batch and RunCompiler modes.
	FileSystem::Init();

	// what mode are we in?
	switch (mode) {
		case MODE_MODELCOMPILER: {
			std::string modelName;
			std::string filePath;
			if (argc > 2) {
				filePath = modelName = argv[2];
				SetupRenderer();
				RunCompiler(modelName, filePath);
			}
			break;
		}

		case MODE_VERSION: {
			std::string version(PIONEER_VERSION);
			if (strlen(PIONEER_EXTRAVERSION)) version += " (" PIONEER_EXTRAVERSION ")";
			Output("texturecompiler %s\n", version.c_str());
			break;
		}

		case MODE_USAGE_ERROR:
			Output("texturecompiler: unknown mode %s\n", argv[1]);
			// fall through

		case MODE_USAGE:
			Output(
				"usage: texturecompiler [mode] [options...]\n"
				"available modes:\n"
				"    -compile          [-c ...]          texture compiler\n"
				"    -version          [-v]              show version\n"
				"    -help             [-h,-?]           this help\n"
			);
			break;
	}

	return 0;
}
