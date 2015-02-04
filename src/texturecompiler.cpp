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

namespace encodings 
{
	class EncDec {
	public:
		virtual vector2f encode(const vector3f&) = 0;
		virtual vector3f decode(const vector4f&) = 0;
	};

	class SpheremapTransform : public EncDec {
	public:
		virtual vector2f encode(const vector3f& n) override
		{
			const vector2f nxy(n.x, n.y);
			vector2f enc = nxy.Normalized() * (sqrt(-n.z*0.5 + 0.5));
			enc = (enc * 0.5f) + vector2f(0.5f);
			return enc;
		}

		virtual vector3f decode(const vector4f& enc) override
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
		}
	};

	class LambertAzimuthalEqualAreaProjection : public EncDec {
	public:
		virtual vector2f encode(const vector3f& n) override
		{
			const float f = sqrt(8.0f * n.z + 8.0f);
			return n.xy() / f + 0.5f;
		}

		virtual vector3f decode(const vector4f& enc) override
		{
			vector2f fenc = (enc.xy() * 4.0f) - 2.0f;
			float f = fenc.Dot(fenc);
			float g = sqrt(1.0f - f / 4.0f);
			vector3f n;
			n.xy() = fenc*g;
			n.z = 1.0f - f / 2.0f;
			return n;
		}
	};

	class SphericalCoordinates : public EncDec {
	public:
		#define kPI 3.1415926536f
		virtual vector2f encode(const vector3f& n) override
		{
			// convert to spherical coordinates, adjust range to 0..1 from -1..1
			return (vector2f(atan2(n.y, n.x) / kPI, n.z) + vector2f(1.0f)) * vector2f(0.5f);
		}

		virtual vector3f decode(const vector4f& enc) override
		{
			vector2f ang = (enc.xy() * 2.0f) - 1.0f;
			vector2f scth;
			sincos(ang.x * kPI, scth.x, scth.y);
			vector2f scphi = vector2f(sqrt(1.0f - ang.y*ang.y), ang.y);
			return vector3f(scth.y*scphi.x, scth.x*scphi.x, scphi.y);
		}
	private:
		void sincos(const float x, float &sinval, float &cosval)
		{
			sinval = sin(x);
			cosval = sqrt(1.0 - sinval * sinval);
		}
	};

	class StereographicProjection : public EncDec {
	public:
		virtual vector2f encode(const vector3f& n) override
		{
			const float scale = 1.7777f;
			vector2f enc = vector2f(n.x, n.y) / (n.z + 1.0f);
			enc /= scale;
			enc = (enc * 0.5f) + 0.5f;
			return enc;
		}

		virtual vector3f decode(const vector4f& enc) override
		{
			const float scale = 1.7777f;
			const vector3f nn = vector3f(enc.x, enc.y, 0.0f) * vector3f(2.0f * scale, 2.0f * scale, 0.0f) + vector3f(-scale, -scale, 1.0f);
			const float g = 2.0f / nn.Dot(nn);
			vector3f n;
			n.x = g * nn.x;
			n.y = g * nn.y;
			n.z = g - 1.0f;
			return n;
		}
	};
}
using namespace encodings;

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

template <typename T>
struct TImgData {
	TImgData() : data(nullptr), w(-1), h(-1) {}
	~TImgData() {}

	bool LoadImageData(const std::string &imageName, FileSystem::FileSource &fs = FileSystem::gameDataFiles)
	{
		Profiler::Timer timer;
		timer.Start();
		Output("\n---\nLoading image (%s)\n", imageName.c_str());
		try {
			SDLSurfacePtr im = LoadSurfaceFromFile(imageName, fs);
			if (im) {
				// now that we have our raw image loaded
				// allocate the space for our processed representation
				data = new T[(im->w * im->h)];
				T* pimg = data;

				// lock the image once so we can read from it
				SDL_LockSurface(im.Get());

				// setup our map dimensions for later
				w = im->w;
				h = im->h;
				const int bpp = (im->pitch / im->w);

				// copy every pixel value from the red channel (image is greyscale, channel is irrelevant)
				for (int y = 0; y<im->h; y++) {
					for (int x = 0; x<im->w; x++) {
						for(int b = 0; b<bpp; b++) {
							const unsigned char v = static_cast<unsigned char*>(im->pixels)[((x * bpp) + b) + (y*im->pitch)];
							pimg[x + (y*im->w)][b] = v / 255.0f;
						}
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
	T* data;
	int w;
	int h;
};
typedef TImgData<vector3f> TImgData3f;
typedef TImgData<vector4f> TImgData4f;

void EncodeNormalData(vector2f *outNormal, const TImgData3f imgData) {
	SphericalCoordinates encoder;
	const int xmax = imgData.w;
	const int ymax = imgData.h;
	for (int y = 0; y < ymax; y++) {
		for (int x = 0; x < xmax; x++) {
			const vector3f val = imgData.data[(x + (y * xmax))].Normalized();
			const vector2f enc = encoder.encode(val);
			outNormal[(x + (y * xmax))] = enc;
		}
	}
}

void DecodeNormalData(vector3f *outNormal, const TImgData4f imgData) {
	SphericalCoordinates encoder;
	const int xmax = imgData.w;
	const int ymax = imgData.h;
	for (int y = 0; y < ymax; y++) {
		for (int x = 0; x < xmax; x++) {
			const vector4f val = imgData.data[(x + (y * xmax))];
			const vector3f dec = encoder.decode(val);
			outNormal[(x + (y * xmax))] = dec.Normalized();
		}
	}
}

void AverageRGBData(float *outIntensity, const TImgData3f imgData) {
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
			outRGBAf[idx] = vector4f(rg[idx].x, rg[idx].y,  b[idx],  a[idx]);
		}
	}
}

static const std::string dir("textures-compiled");
void SaveData(const std::string &filename, vector4f *packed, const int w, const int h)
{
	FileSystem::userFiles.MakeDirectory(dir);
	// copy data into pixel friendly format / values
	std::unique_ptr<Uint8> pixels(new Uint8[(w * h) * 4]);
	// pad rows to 4 bytes, which is the default row alignment for OpenGL
	const int stride = (4 * w + 3) & ~3;
	for (int y = 0; y<h; y++) {
		for (int x = 0; x<w; x++) {
			const int idx((x * 4) + (y * stride));
			const vector4f &v4 = packed[x + (((h - 1) - y) * w)];
			assert(0.0f <= v4.x && 0.0f <= v4.y);
			assert(1.0f >= v4.x && 1.0f >= v4.y);
			pixels.get()[idx + 0] = Uint8(v4.x * 255.0f);
			pixels.get()[idx + 1] = Uint8(v4.y * 255.0f);
			pixels.get()[idx + 2] = Uint8(v4.z * 255.0f);
			pixels.get()[idx + 3] = Uint8(v4.w * 255.0f);
		}
	}
	// write to png file
	const std::string fname = FileSystem::JoinPathBelow(dir, filename);
	write_png(FileSystem::userFiles, fname, pixels.get(), w, h, stride, 4);
}

void RunCompiler(const std::string &srcFolder, const std::string &diffuseName, const std::string &normalName, const std::string &specularName, const std::string &AOName, const std::string &prepend)
{
	Profiler::Timer timer;
	timer.Start();
	Output("\n---\nStarting texture compiler for (%s, %s, %s, %s)\n", diffuseName.c_str(), normalName.c_str(), specularName.c_str(), AOName.c_str());

	//load the images, encode & pack them, save the resulting versions
	TImgData3f diffuseImg;
	TImgData3f normalImg;
	TImgData3f specularImg;
	TImgData3f aoImg;
	const bool resDif = diffuseImg.LoadImageData(FileSystem::JoinPathBelow(srcFolder,diffuseName));
	const bool resNor = normalImg.LoadImageData(FileSystem::JoinPathBelow(srcFolder,normalName));
	const bool resSpe = specularImg.LoadImageData(FileSystem::JoinPathBelow(srcFolder,specularName));
	const bool resAmb = aoImg.LoadImageData(FileSystem::JoinPathBelow(srcFolder,AOName));

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
		SaveData(prepend + "diffusePacked.png", packedDiffuse.get(), diffuseImg.w, diffuseImg.h);
		SaveData(prepend + "nsAOPacked.png", packedNSAO.get(), normalImg.w, normalImg.h);
#if 0 
		// test normal map decoding
		TImgData4f nsAO;
		const std::string fname = FileSystem::JoinPathBelow(dir, prepend + "nsAOPacked.png");
		nsAO.LoadImageData(fname, FileSystem::userFiles);
		std::unique_ptr<vector3f[]> decNormalData( new vector3f[(nsAO.w * nsAO.h)] );
		DecodeNormalData(decNormalData.get(), nsAO);
		for (int y = 0; y < normalImg.h; y++) {
			for (int x = 0; x < normalImg.w; x++) {
				const vector3f val = normalImg.data[(x + (y * normalImg.w))].Normalized();
				const vector3f dec = decNormalData[(x + (y * normalImg.w))];
				const float diff = val.Dot(dec);
				if( diff < 0.9f ) {
					Output("Big diff (%f)\n", diff);
				}
			}
		}
#endif
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
			std::string srcFolder;
			std::string diffuseName;
			std::string normalName;
			std::string specularName;
			std::string AOName;
			std::string prepend;
			if (argc > 7) {
				SetupBasics();
				//SetupRenderer();
				srcFolder = argv[2];
				diffuseName = argv[3];
				normalName = argv[4];
				specularName = argv[5];
				AOName = argv[6];
				prepend = argv[7];
				RunCompiler(srcFolder, diffuseName, normalName, specularName, AOName, prepend);
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
