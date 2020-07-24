// Copyright © 2008-2020 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "buildopts.h"
#include "core/Log.h"
#include "libs.h"
#include "utils.h"
#include <cstdio>
#include <cstdlib>

#include "scenegraph/SceneGraph.h"

#include "FileSystem.h"
#include "GameConfig.h"
#include "GameSaveError.h"
#include "JobQueue.h"
#include "ModManager.h"
#include "StringF.h"
#include "core/OS.h"
#include "graphics/Drawables.h"
#include "graphics/Graphics.h"
#include "graphics/Light.h"
#include "graphics/Renderer.h"
#include "graphics/Texture.h"
#include "graphics/TextureBuilder.h"
#include "graphics/VertexArray.h"
#include "graphics/dummy/RendererDummy.h"
#include "scenegraph/BinaryConverter.h"
#include "scenegraph/DumpVisitor.h"
#include "scenegraph/FindNodeVisitor.h"
#include <sstream>

std::unique_ptr<GameConfig> s_config;
std::unique_ptr<Graphics::Renderer> s_renderer;

//#define USES_THREADS
#ifdef USES_THREADS
std::unique_ptr<AsyncJobQueue> asyncJobQueue;
#endif

static const std::string s_dummyPath("");

// fwd decl'
void RunCompiler(const std::string &modelName, const std::string &filepath, const bool bInPlace);
void CompileModel(const bool isInPlace, const std::string &modelName);
void BatchExport(const bool isInPlace);

// ********************************************************************************
// Overloaded PureJob class to handle compiling each model
// ********************************************************************************
class CompileJob : public Job {
public:
	CompileJob(){};
	CompileJob(const std::string &name, const std::string &path, const bool inPlace) :
		m_name(name),
		m_path(path),
		m_inPlace(inPlace) {}

	virtual void OnRun() override final { RunCompiler(m_name, m_path, m_inPlace); } // RUNS IN ANOTHER THREAD!! MUST BE THREAD SAFE!
	virtual void OnFinish() override final {}
	virtual void OnCancel() override final {}

protected:
	std::string m_name;
	std::string m_path;
	bool m_inPlace;
};

// ********************************************************************************
// functions
// ********************************************************************************
void SetupRenderer()
{
	PROFILE_SCOPED()
	s_config.reset(new GameConfig);

	Log::GetLog()->SetLogFile("modelcompiler.log");

	//init components
	FileSystem::userFiles.MakeDirectory(""); // ensure the config directory exists
	static const Uint32 sdl_init_nothing = 0;
	if (SDL_Init(sdl_init_nothing) < 0)
		Error("SDL initialization failed: %s\n", SDL_GetError());

	ModManager::Init();

	Graphics::RendererDummy::RegisterRenderer();

	//video
	Graphics::Settings videoSettings = {};
	videoSettings.rendererType = Graphics::RENDERER_DUMMY;
	videoSettings.width = s_config->Int("ScrWidth");
	videoSettings.height = s_config->Int("ScrHeight");
	videoSettings.fullscreen = false;
	videoSettings.hidden = true;
	videoSettings.requestedSamples = s_config->Int("AntiAliasingMode");
	videoSettings.vsync = false;
	videoSettings.useTextureCompression = true;
	videoSettings.useAnisotropicFiltering = true;
	videoSettings.iconFile = OS::GetIconFilename();
	videoSettings.title = "Model Compiler";
	s_renderer.reset(Graphics::Init(videoSettings));

#ifdef USES_THREADS
	// get threads up
	Uint32 numThreads = s_config->Int("WorkerThreads");
	const int numCores = OS::GetNumCores();
	assert(numCores > 0);
	if (numThreads == 0)
		numThreads = std::max(Uint32(numCores), 1U); // this is a tool, we can use all of the cores for processing unlike Pioneer
	asyncJobQueue.reset(new AsyncJobQueue(numThreads));
	Output("started %d worker threads\n", numThreads);
#endif
}

void RunCompiler(const std::string &modelName, const std::string &filepath, const bool bInPlace)
{
	PROFILE_SCOPED()
	Profiler::Timer timer;
	timer.Start();
	Output("\n---\nStarting compiler for (%s)\n", modelName.c_str());

	//load the current model in a pristine state (no navlights, shields...)
	//and then save it into binary
	std::unique_ptr<SceneGraph::Model> model;
	try {
		SceneGraph::Loader ld(s_renderer.get(), true, false);
		model.reset(ld.LoadModel(modelName));
		//dump warnings
		for (std::vector<std::string>::const_iterator it = ld.GetLogMessages().begin();
			 it != ld.GetLogMessages().end(); ++it) {
			Output("%s\n", (*it).c_str());
		}
	} catch (...) {
		//minimal error handling, this is not expected to happen since we got this far.
		return;
	}

	try {
		const std::string DataPath = FileSystem::NormalisePath(filepath.substr(0, filepath.size() - 6));
		SceneGraph::BinaryConverter bc(s_renderer.get());
		bc.Save(modelName, DataPath, model.get(), bInPlace);
	} catch (const CouldNotOpenFileException &) {
	} catch (const CouldNotWriteToFileException &) {
	}

	timer.Stop();
	Output("Compiling \"%s\" took: %lf\n", modelName.c_str(), timer.millicycles());
}

void CompileModel(const bool isInPlace, const std::string &modelName)
{
	std::string filePath;
	filePath = modelName; // if it's isInPlace then this will get overwritten below.

	// determine if we're meant to be writing these in the source directory
	if (isInPlace) {
		// find all of the models
		FileSystem::FileSource &fileSource = FileSystem::gameDataFiles;
		for (FileSystem::FileEnumerator files(fileSource, "models", FileSystem::FileEnumerator::Recurse); !files.Finished(); files.Next()) {
			const FileSystem::FileInfo &info = files.Current();
			const std::string &fpath = info.GetPath();

			//check it's the expected type
			if (info.IsFile()) {
				if (ends_with_ci(fpath, ".model")) { // store the path for ".model" files
					const std::string shortname(info.GetName().substr(0, info.GetName().size() - 6));
					if (shortname == modelName) {
						filePath = fpath;
						break;
					}
				}
			}
		}
	}
	SetupRenderer();
	RunCompiler(modelName, filePath, isInPlace);
}

void BatchExport(const bool isInPlace)
{
	// find all of the models
	std::vector<std::pair<std::string, std::string>> list_model;
	FileSystem::FileSource &fileSource = FileSystem::gameDataFiles;
	for (FileSystem::FileEnumerator files(fileSource, "models", FileSystem::FileEnumerator::Recurse); !files.Finished(); files.Next()) {
		const FileSystem::FileInfo &info = files.Current();
		const std::string &fpath = info.GetPath();

		//check it's the expected type
		if (info.IsFile()) {
			if (ends_with_ci(fpath, ".model")) { // store the path for ".model" files
				list_model.push_back(std::make_pair(info.GetName().substr(0, info.GetName().size() - 6), fpath));
			}
		}
	}

	SetupRenderer();

#ifdef USES_THREADS
	std::deque<Job::Handle> handles;
	for (auto &modelName : list_model) {
		handles.push_back(asyncJobQueue->Queue(new CompileJob(modelName.first, modelName.second, isInPlace)));
	}

	while (true) {
		asyncJobQueue->FinishJobs();
		bool hasJobs = false;
		for (auto &handle : handles)
			hasJobs |= handle.HasJob();

		if (!hasJobs)
			break;
	}
#else
	for (auto &modelName : list_model) {
		RunCompiler(modelName.first, modelName.second, isInPlace);
	}
#endif
}

// ********************************************************************************
// functions
// ********************************************************************************
enum class RunMode {
	MODE_MODELCOMPILER = 0,
	MODE_MODELBATCHEXPORT,
	MODE_USAGE_ERROR
};

static FileSystem::FileSourceFS customDataDir(".");

extern "C" int main(int argc, char **argv)
{
#ifdef PIONEER_PROFILER
	Profiler::detect(argc, argv);
#endif

	bool showUsage = false;
	bool showVersion = false;
	bool isInPlace = false;
	bool hasCustomDir = false;
	bool hasError = false;
	RunMode mode = RunMode::MODE_MODELCOMPILER;
	std::string modelName; // used by MODE_MODELCOMPILER
	std::string modeString;

	if (argc > 1) {
		// start at 1 as the first argument is just the path+filename to the exe being run
		for (int argi = 1; argi < argc && !hasError; argi++) {
			const char switchchar = argv[argi][0];
			if (!(switchchar == '-' || switchchar == '/')) {
				mode = RunMode::MODE_USAGE_ERROR;
			}

			const std::string argopt(std::string(argv[argi]).substr(1));

			if (argopt == "inplace") {
				isInPlace = true;
			} else if (argopt == "compile" || argopt == "c" && argc > argi + 1) {
				++argi; // grab the next arg
				const std::string modelArg(argv[argi]);
				if (!modelArg.empty() && modelArg[0] != '-' && modelArg[0] != '/') {
					mode = RunMode::MODE_MODELCOMPILER;
					modeString = argopt;
					modelName = modelArg;
				} else {
					mode = RunMode::MODE_USAGE_ERROR;
					modeString = argopt;
					hasError = true;
				}
			} else if (argopt == "batch" || argopt == "b") {
				mode = RunMode::MODE_MODELBATCHEXPORT;
				modeString = argopt;
			} else if (argopt == "adddir" && argc > argi + 1) {
				++argi; // grab the next arg
				const std::string dirArg(argv[argi]);
				if (!dirArg.empty() && dirArg[0] != '-' && dirArg[0] != '/') {
					customDataDir = FileSystem::FileSourceFS(dirArg);
					hasCustomDir = true;
					isInPlace = true;
				} else {
					Output("modelcompiler: adddir passed empty dir\n");
					hasError = true;
				}
			} else if (argopt == "version" || argopt == "v") {
				showVersion = true;
			} else if (argopt == "help" || argopt == "h" || argopt == "?") {
				showUsage = true;
			} else {
				mode = RunMode::MODE_USAGE_ERROR;
				modeString = argopt;
			}
		}
	} else {
		showVersion = true;
		mode = RunMode::MODE_USAGE_ERROR;
	}

	// Init here since we'll need it for both batch and RunCompiler modes.
	FileSystem::Init();
	FileSystem::userFiles.MakeDirectory(""); // ensure the config directory exists
	// setup a custom source directory if one hs been added
	if (isInPlace && hasCustomDir) {
		FileSystem::gameDataFiles.AppendSource(&customDataDir);
	}
#ifdef PIONEER_PROFILER
	FileSystem::userFiles.MakeDirectory("profiler");
	const std::string profilerPath = FileSystem::JoinPathBelow(FileSystem::userFiles.GetRoot(), "profiler");
#endif

	// what mode are we in?
	switch (mode) {
	case RunMode::MODE_MODELCOMPILER: {
		CompileModel(isInPlace, modelName);
		break;
	}

	case RunMode::MODE_MODELBATCHEXPORT:
		BatchExport(isInPlace);
		break;

	case RunMode::MODE_USAGE_ERROR:
		Output("modelcompiler: unknown mode %s\n", modeString.c_str());
		showUsage = true;
		break;
	}

	if (showVersion) {
		std::string version(PIONEER_VERSION);
		if (strlen(PIONEER_EXTRAVERSION)) version += " (" PIONEER_EXTRAVERSION ")";
		Output("modelcompiler ver: %s, SGM ver: %u\n", version.c_str(), SceneGraph::BinaryConverter::GetSGMVersion());
	}

	if (showUsage) {
		Output(
			"usage: modelcompiler [mode] [options...]\n"
			"available modes:\n"
			"    -compile          [-c <filename>]       model compiler\n"
			"    -batch            [-b]                  batch mode\n"
			"    -inplace          [-inplace]            output into the source folder, default is users /home/Pioneer directory\n"
			"    -adddir           [-adddir <dir>]       add a custom source directory"
			"    -version          [-v]                  show version\n"
			"    -help             [-h,-?]               this help\n");
	}

#ifdef PIONEER_PROFILER
	Profiler::dumphtml(profilerPath.c_str());
#endif

	Graphics::Uninit();
	SDL_Quit();
	FileSystem::Uninit();
#ifdef USES_THREADS
	asyncJobQueue.reset();
#endif
	//exit(0);

	return 0;
}
