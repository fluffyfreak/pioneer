// Copyright © 2008-2017 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "libs.h"
#include "utils.h"
#include <cstdio>
#include <cstdlib>

#include "scenegraph/SceneGraph.h"

#include "FileSystem.h"
#include "StringRange.h"
#include "GameConfig.h"
#include "JobQueue.h"
#include "graphics/dummy/RendererDummy.h"
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
std::unique_ptr<AsyncJobQueue> asyncJobQueue;

static const std::string s_dummyPath("");

// fwd decl'
void RunCompiler(const std::string &modelName, const std::string &filepath, const bool bInPlace);

// ********************************************************************************
// functions
// ********************************************************************************
void Setup()
{
	PROFILE_SCOPED()
	s_config.reset(new GameConfig);

	OS::RedirectStdio();

	//init components
	FileSystem::userFiles.MakeDirectory(""); // ensure the config directory exists

	ModManager::Init();
}

struct HYGField
{
	HYGField() { memset(this, 0, sizeof(HYGField)); }
	int id ; // : The database primary key.
	int hip ; // : The star's ID in the Hipparcos catalog, if known.
	int hd ; // : The star's ID in the Henry Draper catalog, if known.
	int hr ; // : The star's ID in the Harvard Revised catalog, which is the same as its number in the Yale Bright Star Catalog.
	int gl ; // : The star's ID in the third edition of the Gliese Catalog of Nearby Stars.
	std::string bf ; // : The Bayer / Flamsteed designation, primarily from the Fifth Edition of the Yale Bright Star Catalog.This is a combination of the two designations.The Flamsteed number, if present, is given first; then a three - letter abbreviation for the Bayer Greek letter; the Bayer superscript number, if present; and finally, the three - letter constellation abbreviation.Thus Alpha Andromedae has the field value "21Alp And", and Kappa1 Sculptoris (no Flamsteed number) has "Kap1Scl".
	float ra, dec; // : The star's right ascension and declination, for epoch and equinox 2000.0.
	std::string proper ; // : A common name for the star, such as "Barnard's Star" or "Sirius".I have taken these names primarily from the Hipparcos project's web site, which lists representative names for the 150 brightest stars and many of the 150 closest stars. I have added a few names to this list. Most of the additions are designations from catalogs mostly now forgotten (e.g., Lalande, Groombridge, and Gould ["G."]) except for certain nearby stars which are still best known by these designations.
	float dist ; // : The star's distance in parsecs, the most common unit in astrometry. To convert parsecs to light years, multiply by 3.262. A value >= 10000000 indicates missing or dubious (e.g., negative) parallax data in Hipparcos.
	float pmra, pmdec ; // : The star's proper motion in right ascension and declination, in milliarcseconds per year.
	float rv ; // : The star's radial velocity in km/sec, where known.
	float mag ; // : The star's apparent visual magnitude.
	float absmag ; // : The star's absolute visual magnitude (its apparent magnitude from a distance of 10 parsecs).
	float spect ; // : The star's spectral type, if known.
	float ci ; // : The star's color index (blue magnitude - visual magnitude), where known.
	float x, y, z ; // : The Cartesian coordinates of the star, in a system based on the equatorial coordinates as seen from Earth. + X is in the direction of the vernal equinox (at epoch 2000), +Z towards the north celestial pole, and +Y in the direction of R.A. 6 hours, declination 0 degrees.
	float vx, vy, vz ; // : The Cartesian velocity components of the star, in the same coordinate system described immediately above.They are determined from the proper motion and the radial velocity (when known).The velocity unit is parsecs per year; these are small values (around 1 millionth of a parsec per year), but they enormously simplify calculations using parsecs as base units for celestial mapping.
	float rarad, decrad, pmrarad, prdecrad; // : The positions in radians, and proper motions in radians per year.
	float bayer ; // : The Bayer designation as a distinct value
	float flam ; // : The Flamsteed number as a distinct value
	std::string con ; // : The standard constellation abbreviation
	int comp, comp_primary, base ; // : Identifies a star in a multiple star system.comp = ID of companion star, comp_primary = ID of primary star for this component, and base = catalog ID or name for this multi - star system.Currently only used for Gliese stars.
	float lum ; // : Star's luminosity as a multiple of Solar luminosity.
	std::string var ; // : Star's standard variable star designation, when known.
	float var_min, var_max ; // : Star's approximate magnitude range, for variables. This value is based on the Hp magnitudes for the range in the original Hipparcos catalog, adjusted to the V magnitude scale to match the "mag" field.
};
#pragma optimize("",off)
std::string ReplaceAll(std::string str, const std::string& from, const std::string& to) {
	size_t start_pos = 0;
	while((start_pos = str.find(from, start_pos)) != std::string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos += to.length()-1; // Handles case where 'to' is a substring of 'from'
	}
	return str;
}

static inline size_t SplitSpec(const std::string &spec, std::vector<std::string> &output)
{
	static const std::string delim(",");

	std::string specClone(ReplaceAll(spec, ",,", ",null,"));
	specClone = ReplaceAll(specClone, ",\n", ",null\n");
	if (specClone.back() == ',')
		specClone += "null";

	size_t i = 0, start = 0, end = 0;
	while (end != std::string::npos) {
		// get to the first non-delim char
		start = specClone.find_first_not_of(delim, end);

		// read the end, no more to do
		if (start == std::string::npos)
			break;

		// find the end - next delim or end of string
		end = specClone.find_first_of(delim, start);

		// extract the fragment and remember it
		output[i++] = (specClone.substr(start, (end == std::string::npos) ? std::string::npos : end - start).c_str());
	}

	return i;
}

void RunCompiler(const std::string &modelName, const std::string &filepath)
{
	PROFILE_SCOPED()
	Profiler::Timer timer;
	timer.Start();
	Output("\n---\nStarting compiler for (%s)\n", modelName.c_str());

	std::vector<HYGField> HYGFields;
	HYGFields.reserve(120000);

	FileSystem::userFiles.MakeDirectory("starcharts");

	try
	{
		auto file = FileSystem::userFiles.ReadFile("starcharts/" + filepath);
		if (!file)
			throw CouldNotOpenFileException();

		// parse csv line by line
		StringRange buffer = file->AsStringRange();
		buffer = buffer.StripUTF8BOM();

		while (!buffer.Empty())
		{
			StringRange line = buffer.ReadLine().StripSpace();

			// if the line is a comment, skip it
			if (line.Empty() || (line[0] == '#'))
				continue;

			//
			std::vector<std::string> split(37);
			size_t numFound = SplitSpec(line.ToString(), split);
			printf("%zu\n", numFound);

			//id,hip,hd,hr,gl,bf,proper,ra,dec,dist,pmra,pmdec,rv,mag,absmag,spect,ci,x,y,z,vx,vy,vz,rarad,decrad,pmrarad,pmdecrad,bayer,flam,con,comp,comp_primary,base,lum,var,var_min,var_max

			HYGField hf;
			hf.id = atoi(split[0].c_str());
			hf.x = atoi(split[17].c_str());
			hf.y = atoi(split[18].c_str());
			hf.z = atoi(split[19].c_str());
			HYGFields.push_back(hf);
		}
	}
	catch (const CouldNotOpenFileException&)
	{
		// error
	}
	catch (const CouldNotWriteToFileException&)
	{
		// error
	}

	timer.Stop();
	Output("Compiling \"%s\" took: %lf\n", modelName.c_str(), timer.millicycles());
}
#pragma optimize("",on)

// ********************************************************************************
// functions
// ********************************************************************************
enum RunMode {
	MODE_CHARTCOMPILER=0,
	MODE_VERSION,
	MODE_USAGE,
	MODE_USAGE_ERROR
};

int main(int argc, char** argv)
{
#ifdef PIONEER_PROFILER
	Profiler::detect( argc, argv );
#endif

	RunMode mode = MODE_CHARTCOMPILER;

	if (argc > 1) {
		const char switchchar = argv[1][0];
		if (!(switchchar == '-' || switchchar == '/')) {
			mode = MODE_USAGE_ERROR;
			goto start;
		}

		const std::string modeopt(std::string(argv[1]).substr(1));

		if (modeopt == "filename" || modeopt == "f") {
			mode = MODE_CHARTCOMPILER;
			goto start;
		}

		if (modeopt == "version" || modeopt == "v") {
			mode = MODE_VERSION;
			goto start;
		}

		if (modeopt == "help" || modeopt == "h" || modeopt == "?") {
			mode = MODE_USAGE;
			goto start;
		}

		mode = MODE_USAGE_ERROR;
	}

start:

	// Init here since we'll need it for both batch and RunCompiler modes.
	FileSystem::Init();
	FileSystem::userFiles.MakeDirectory(""); // ensure the config directory exists
#ifdef PIONEER_PROFILER
	FileSystem::userFiles.MakeDirectory("profiler");
	const std::string profilerPath = FileSystem::JoinPathBelow(FileSystem::userFiles.GetRoot(), "profiler");
#endif

	// what mode are we in?
	switch (mode) {
		case MODE_CHARTCOMPILER: {
			std::string modelName;
			std::string filePath;
			if (argc > 2) {
				filePath = modelName = argv[2];
				// determine if we're meant to be writing these in the source directory
				Setup();
				RunCompiler(modelName, filePath);
			}
			break;
		}

		case MODE_VERSION: {
			std::string version(PIONEER_VERSION);
			if (strlen(PIONEER_EXTRAVERSION)) version += " (" PIONEER_EXTRAVERSION ")";
			Output("starcharts %s\n", version.c_str());
			break;
		}

		case MODE_USAGE_ERROR:
			Output("starcharts: unknown mode %s\n", argv[1]);
			// fall through

		case MODE_USAGE:
			Output(
				"usage: starcharts [mode] [options...]\n"
				"available modes:\n"
				"    -compile          [-c ...]          star chart compiler\n"
				"    -version          [-v]              show version\n"
				"    -help             [-h,-?]           this help\n"
			);
			break;
	}

#ifdef PIONEER_PROFILER
	Profiler::dumphtml(profilerPath.c_str());
#endif

	return 0;
}
