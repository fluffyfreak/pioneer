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
#include "PngWriter.h"
#include <sstream>

std::unique_ptr<GameConfig> s_config;
//std::unique_ptr<Graphics::Renderer> s_renderer;

static const std::string s_dummyPath("");

void SetupBasics()
{
	s_config.reset(new GameConfig);

	OS::RedirectStdio();

	//init components
	FileSystem::userFiles.MakeDirectory(""); // ensure the config directory exists
	if (SDL_Init(SDL_INIT_VIDEO) < 0)
		Error("SDL initialization failed: %s\n", SDL_GetError());

	ModManager::Init();
}

/*void SetupRenderer()
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
}*/

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

struct TImgData {
	TImgData() : data(nullptr), w(-1), h(-1) {}
	~TImgData() {}
	vector3f* data;
	int w;
	int h;
};

bool LoadImageData(const std::string &imageName, TImgData &imgData)
{
	Profiler::Timer timer;
	timer.Start();
	Output("\n---\nLoading image (%s)\n", imageName.c_str());
	try {
		SDLSurfacePtr im = LoadSurfaceFromFile(imageName);
		if (im) {
			// now that we have our raw image loaded
			// allocate the space for our processed representation
			imgData.data = new vector3f[(im->w * im->h)];
			vector3f* pimg = imgData.data;

			// lock the image once so we can read from it
			SDL_LockSurface(im.Get());

			// setup our map dimensions for later
			imgData.w = im->w;
			imgData.h = im->h;

			// copy every pixel value from the red channel (image is greyscale, channel is irrelevant)
			for (int y = 0; y<im->h; y++) {
				for (int x = 0; x<im->w; x++) {
					const unsigned char v0 = static_cast<unsigned char*>(im->pixels)[((x * 3) + 0) + (y*im->pitch)];
					const unsigned char v1 = static_cast<unsigned char*>(im->pixels)[((x * 3) + 1) + (y*im->pitch)];
					const unsigned char v2 = static_cast<unsigned char*>(im->pixels)[((x * 3) + 2) + (y*im->pitch)];
					pimg[x + (y*im->w)] = vector3f(v0 / 255.0f, v1 / 255.0f, v2 / 255.0f);
				}
			}

			// unlock the surface and then release it
			SDL_UnlockSurface(im.Get());
		} else {
			Output("Failed to load image.\n");
		}
	} catch (...) {
		//minimal error handling, this is not expected to happen since we got this far.
		timer.Stop();
		return false;
	}

	timer.Stop();
	Output("Loading \"%s\" took: %lf\n", imageName.c_str(), timer.millicycles());

	return true;
}

void EncodeNormalData(vector2f *outNormal, const TImgData imgData) {
	const int xmax = imgData.w;
	const int ymax = imgData.h;
	for (int y = 0; y < ymax; y++) {
		for (int x = 0; x < xmax; x++) {
			const vector3f val = imgData.data[(x + (y * xmax))].Normalized();
			const vector2f enc = encode(val);
			//const vector3f dec = decode(vector4f(enc.x, enc.y, 0.0f, 0.0f));
			//const vector3f decN = dec.Normalized();
			outNormal[(x + (y * xmax))] = enc;
		}
	}
}

void AverageRGBData(float *outIntensity, const TImgData imgData) {
	const int xmax = imgData.w;
	const int ymax = imgData.h;
	for (int y = 0; y < ymax; y++) {
		for (int x = 0; x < xmax; x++) {
			const vector3f val = imgData.data[(x + (y * xmax))];
			// average the RGB value to get a greyscale "Intensity" like value
			outIntensity[(x + (y * xmax))] = (val.x + val.y + val.z) / 3.0f;
		}
	}
}

void PackRGBandA(vector4f *outRGBAf, const vector3f *rgb, const float *a, const int w, const int h)
{
	for(int x=0; x<w; x++) {
		for(int y=0; y<h; y++) {
			const int idx(x + (y * w));
			outRGBAf[idx] = vector4f(rgb[idx], a[idx]);
		}
	}
}

void PackRGandBandA(vector4f *outRGBAf, const vector2f *rg, const float *b, const float *a, const int w, const int h)
{
	for(int x=0; x<w; x++) {
		for(int y=0; y<h; y++) {
			const int idx(x + (y * w));
			outRGBAf[idx] = vector4f(rg[idx].x, rg[idx].y, b[idx], a[idx]);
		}
	}
}

void RunCompiler(const std::string &diffuseName, const std::string &normalName, const std::string &specularName, const std::string &AOName, const std::string &prepend)
{
	Profiler::Timer timer;
	timer.Start();
	Output("\n---\nStarting texture compiler for (%s, %s, %s, %s)\n", diffuseName.c_str(), normalName.c_str(), specularName.c_str(), AOName.c_str());

	//load the images, encode & pack them, save the resulting versions
	TImgData diffuseImg;
	TImgData normalImg;
	TImgData specularImg;
	TImgData aoImg;
	const bool resDif = LoadImageData(diffuseName, diffuseImg);
	const bool resNor = LoadImageData(normalName, normalImg);
	const bool resSpe = LoadImageData(specularName, specularImg);
	const bool resAmb = LoadImageData(AOName, aoImg);

	const bool equalDims = (diffuseImg.w && normalImg.w && specularImg.w && aoImg.w)
		&& (diffuseImg.h && normalImg.h && specularImg.h && aoImg.h);
	
	if(resDif && resNor && resSpe && resAmb && equalDims) {
		// get the greyscale of the diffuse RGB data
		std::unique_ptr<float[]> avgDiffuseData( new float[(diffuseImg.w * diffuseImg.h)] );
		AverageRGBData(avgDiffuseData.get(), diffuseImg);

		// Encode the normal data
		std::unique_ptr<vector2f[]> encNormalData( new vector2f[(normalImg.w * normalImg.h)] );
		EncodeNormalData(encNormalData.get(), normalImg);

		// get the greyscale of the specular RGB data
		std::unique_ptr<float[]> avgSpecData( new float[(specularImg.w * specularImg.h)] );
		AverageRGBData(avgSpecData.get(), specularImg);

		// get the greyscale of the ambient RGB data
		std::unique_ptr<float[]> avgAOData( new float[(aoImg.w * aoImg.h)] );
		AverageRGBData(avgAOData.get(), aoImg);

		// pack diffuseImg (RGB) with avgDiffuseData (A)
		std::unique_ptr<vector4f[]> packedDiffuse( new vector4f[(diffuseImg.w * diffuseImg.h)] );
		PackRGBandA( packedDiffuse.get(), diffuseImg.data, avgDiffuseData.get(), diffuseImg.w, diffuseImg.h);
		
		// pack encNormalData (RG) with avgSpecData (B) with avgAOData (A)
		std::unique_ptr<vector4f[]> packedNSAO( new vector4f[(normalImg.w * normalImg.h)] );
		PackRGandBandA( packedNSAO.get(), encNormalData.get(), avgSpecData.get(), avgAOData.get(), normalImg.w, normalImg.h);

		// save out the data
		static const std::string dir("textures-compiled");
		FileSystem::userFiles.MakeDirectory(dir);
		{
			// copy data into pixel friendly format / values
			const int w = diffuseImg.w;
			const int h = diffuseImg.h;
			std::unique_ptr<Uint8> pixels(new Uint8[(w * h) * 4]);
			// pad rows to 4 bytes, which is the default row alignment for OpenGL
			const int stride = (4*w + 3) & ~3;
			for(int y=0; y<h; y++) {
				for(int x=0; x<w; x++) {
					const int idx((x*4) + (y * stride));
					const vector4f &v4 = packedDiffuse.get()[x + (((h-1)-y) * w)];
					pixels.get()[idx+0] = Uint8(v4.x * 255.0f);
					pixels.get()[idx+1] = Uint8(v4.y * 255.0f);
					pixels.get()[idx+2] = Uint8(v4.z * 255.0f);
					pixels.get()[idx+3] = Uint8(v4.w * 255.0f);
				}
			}
			// write to png file
			const std::string fname = FileSystem::JoinPathBelow(dir, prepend+"diffusePacked.png");
			write_png(FileSystem::userFiles, fname, pixels.get(), w, h, stride, 4);
		}
		
		{
			// copy data into pixel friendly format / values
			const int w = normalImg.w;
			const int h = normalImg.h;
			std::unique_ptr<Uint8> pixels(new Uint8[(w * h) * 4]);
			// pad rows to 4 bytes, which is the default row alignment for OpenGL
			const int stride = (4*w + 3) & ~3;
			for(int y=0; y<h; y++) {
				for(int x=0; x<w; x++) {
					const int idx((x*4) + (y * stride));
					const vector4f &v4 = packedNSAO.get()[x + (((h-1)-y) * w)];
					pixels.get()[idx+0] = Uint8(v4.x * 255.0f);
					pixels.get()[idx+1] = Uint8(v4.y * 255.0f);
					pixels.get()[idx+2] = Uint8(v4.z * 255.0f);
					pixels.get()[idx+3] = Uint8(v4.w * 255.0f);
				}
			}
			// write to png file
			const std::string fname = FileSystem::JoinPathBelow(dir, prepend+"nsAOPacked.png");
			write_png(FileSystem::userFiles, fname, pixels.get(), w, h, stride, 4);
		}
	} 

	timer.Stop();
	Output("Compiling took: %lf\n", timer.millicycles());
}

void DisplayUsage() {
	// -compile textures/asteroid/Stones-Diffuse.png textures/asteroid/Stones-Normal.png textures/asteroid/Stones-Specular.png textures/asteroid/Stones-AO.png Stones-
	// -compile textures/asteroid/Pumice-Diffuse.png textures/asteroid/Pumice-Normal.png textures/asteroid/Pumice-Specular.png textures/asteroid/Pumice-AO.png Pumice-
	Output(
		"usage: texturecompiler [mode] [options...]\n"
		"available modes:\n"
		"    -compile          [-c ...]          texture compiler\n"
		"    -version          [-v]              show version\n"
		"    -help             [-h,-?]           this help\n"
	);
}

enum RunMode {
	MODE_TEXTURECOMPILER=0,
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
			mode = MODE_TEXTURECOMPILER;
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
		case MODE_TEXTURECOMPILER: {
			std::string diffuseName;
			std::string normalName;
			std::string specularName;
			std::string AOName;
			std::string prepend;
			if (argc > 6) {
				SetupBasics();
				//SetupRenderer();
				diffuseName = argv[2];
				normalName = argv[3];
				specularName = argv[4];
				AOName = argv[5];
				prepend = argv[6];
				RunCompiler(diffuseName, normalName, specularName, AOName, prepend);
			} else {
				DisplayUsage();
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
			DisplayUsage();
			break;
	}

	return 0;
}
