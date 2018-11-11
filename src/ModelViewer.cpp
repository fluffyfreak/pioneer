// Copyright © 2008-2018 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "ModelViewer.h"
#include "FileSystem.h"
#include "graphics/opengl/RendererGL.h"
#include "graphics/Graphics.h"
#include "graphics/Light.h"
#include "graphics/TextureBuilder.h"
#include "graphics/Drawables.h"
#include "graphics/VertexArray.h"
#include "scenegraph/DumpVisitor.h"
#include "scenegraph/FindNodeVisitor.h"
#include "scenegraph/BinaryConverter.h"
#include "scenegraph/ModelSkin.h"
#include "OS.h"
#include "Pi.h"
#include "StringF.h"
#include "ModManager.h"
#include "GameSaveError.h"
#include "SH.h"
#include <sstream>

//default options
ModelViewer::Options::Options()
: attachGuns(false)
, showTags(false)
, showDockingLocators(false)
, showCollMesh(false)
, showAabb(false)
, showShields(false)
, showGrid(false)
, showLandingPad(false)
, showUI(true)
, wireframe(false)
, mouselookEnabled(false)
, gridInterval(10.f)
, lightPreset(0)
, orthoView(false)
{
}

//some utility functions
namespace {
	//azimuth/elevation in degrees to a dir vector
	vector3f az_el_to_dir(float yaw, float pitch) {
		//0,0 points to "right" (1,0,0)
		vector3f v;
		v.x = cos(DEG2RAD(yaw)) * cos(DEG2RAD(pitch));
		v.y = sin(DEG2RAD(pitch));
		v.z = sin(DEG2RAD(yaw)) * cos(DEG2RAD(pitch));
		return v;
	}

	//extract color from RGB sliders
	Color get_slider_color(UI::Slider *r, UI::Slider *g, UI::Slider *b) {
		return Color(r->GetValue() * 255.f, g->GetValue() * 255.f, b->GetValue() * 255.f);
	}

	float get_thrust(const UI::Slider *s) {
		return 1.f - (2.f * s->GetValue());
	}

	//add a horizontal button/label pair to a box
	void add_pair(UI::Context *c, UI::Box *box, UI::Widget *widget, const std::string &label) {
		box->PackEnd(c->HBox(5)->PackEnd(UI::WidgetSet( widget, c->Label(label) )));
	}

	void collect_decals(std::vector<std::string> &list)
	{
		const std::string basepath("textures/decals");
		FileSystem::FileSource &fileSource = FileSystem::gameDataFiles;
		for (FileSystem::FileEnumerator files(fileSource, basepath); !files.Finished(); files.Next())
		{
			const FileSystem::FileInfo &info = files.Current();
			const std::string &fpath = info.GetPath();

			//check it's the expected type
			if (info.IsFile() && ends_with_ci(fpath, ".dds")) {
				list.push_back(info.GetName().substr(0, info.GetName().size()-4));
			}
		}
	}

	float zoom_distance(const float base_distance, const float zoom)
	{
		return base_distance * powf(2.0f, zoom);
	}
}

ModelViewer::ModelViewer(Graphics::Renderer *r, LuaManager *lm)
: m_done(false)
, m_screenshotQueued(false)
, m_shieldIsHit(false)
, m_settingColourSliders(false)
, m_shieldHitPan(-1.48f)
, m_frameTime(0.0)
, m_renderer(r)
, m_decalTexture(0)
, m_rotX(0), m_rotY(0), m_zoom(0)
, m_baseDistance(100.0f)
, m_rng(time(0))
, m_currentAnimation(0)
, m_model(0)
, m_modelName("")
{
	OS::RedirectStdio();
	m_ui.Reset(new UI::Context(lm, r, Graphics::GetScreenWidth(), Graphics::GetScreenHeight()));
	m_ui->SetMousePointer("icons/cursors/mouse_cursor_2.png", UI::Point(15, 8));

	m_log = m_ui->MultiLineText("");
	m_log->SetFont(UI::Widget::FONT_SMALLEST);

	m_logScroller.Reset(m_ui->Scroller());
	m_logScroller->SetInnerWidget(m_ui->ColorBackground(Color(0x0,0x0,0x0,0x40))->SetInnerWidget(m_log));

	std::fill(m_mouseButton, m_mouseButton + COUNTOF(m_mouseButton), false);
	std::fill(m_mouseMotion, m_mouseMotion + 2, 0);

	//some widgets
	animSlider = 0;

	onModelChanged.connect(sigc::mem_fun(*this, &ModelViewer::SetupUI));

	//for grid, background
	Graphics::RenderStateDesc rsd;
	rsd.depthWrite = false;
	rsd.cullMode = Graphics::CULL_NONE;
	m_bgState = m_renderer->CreateRenderState(rsd);
}

ModelViewer::~ModelViewer()
{
	ClearModel();
}

void ModelViewer::Run(const std::string &modelName)
{
	std::unique_ptr<GameConfig> config(new GameConfig);

	//init components
	FileSystem::Init();
	FileSystem::userFiles.MakeDirectory(""); // ensure the config directory exists
	if (SDL_Init(SDL_INIT_VIDEO) < 0)
		Error("SDL initialization failed: %s\n", SDL_GetError());

	Lua::Init();

	ModManager::Init();

	Graphics::RendererOGL::RegisterRenderer();

	// determine what renderer we should use, default to Opengl 3.x
	const std::string rendererName = config->String("RendererName", Graphics::RendererNameFromType(Graphics::RENDERER_OPENGL_3x));
	Graphics::RendererType rType = Graphics::RENDERER_OPENGL_3x;
	//if(rendererName == Graphics::RendererNameFromType(Graphics::RENDERER_OPENGL_3x))
	//{
	//	rType = Graphics::RENDERER_OPENGL_3x;
	//}

	//video
	Graphics::Settings videoSettings = {};
	videoSettings.rendererType = rType;
	videoSettings.width = config->Int("ScrWidth");
	videoSettings.height = config->Int("ScrHeight");
	videoSettings.fullscreen = (config->Int("StartFullscreen") != 0);
	videoSettings.hidden = false;
	videoSettings.requestedSamples = config->Int("AntiAliasingMode");
	videoSettings.vsync = (config->Int("VSync") != 0);
	videoSettings.useTextureCompression = (config->Int("UseTextureCompression") != 0);
	videoSettings.useAnisotropicFiltering = (config->Int("UseAnisotropicFiltering") != 0);
	videoSettings.iconFile = OS::GetIconFilename();
	videoSettings.title = "Model viewer";
	Graphics::Renderer *renderer = Graphics::Init(videoSettings);

	NavLights::Init(renderer);
	Shields::Init(renderer);

	//run main loop until quit
	ModelViewer *viewer = new ModelViewer(renderer, Lua::manager);
	viewer->SetModel(modelName);
	viewer->ResetCamera();
	viewer->MainLoop();

	//uninit components
	delete viewer;
	Lua::Uninit();
	delete renderer;
	Shields::Uninit();
	NavLights::Uninit();
	Graphics::Uninit();
	FileSystem::Uninit();
	SDL_Quit();
}

bool ModelViewer::OnPickModel(UI::List *list)
{
	m_requestedModelName = list->GetSelectedOption();
	return true;
}

bool ModelViewer::OnQuit()
{
	m_done = true;
	return true;
}

bool ModelViewer::OnReloadModel(UI::Widget *w)
{
	//camera is not reset, it would be annoying when
	//tweaking materials
	SetModel(m_modelName);
	return true;
}

bool ModelViewer::OnToggleCollMesh(UI::CheckBox *w)
{
	m_options.showDockingLocators = !m_options.showDockingLocators;
	m_options.showCollMesh = !m_options.showCollMesh;
	m_options.showAabb = m_options.showCollMesh;
	return m_options.showCollMesh;
}

bool ModelViewer::OnToggleShowShields(UI::CheckBox *w)
{
	m_options.showShields = !m_options.showShields;
	return m_options.showShields;
}

bool ModelViewer::OnToggleGrid(UI::Widget *)
{
	if (!m_options.showGrid) {
		m_options.showGrid = true;
		m_options.gridInterval = 1.0f;
	}
	else {
		m_options.gridInterval = powf(10, ceilf(log10f(m_options.gridInterval))+1);
		if (m_options.gridInterval >= 10000.0f) {
			m_options.showGrid = false;
			m_options.gridInterval = 0.0f;
		}
	}
	AddLog(m_options.showGrid
		? stringf("Grid: %0{d}", int(m_options.gridInterval))
		: "Grid: off");
	return m_options.showGrid;
}

bool ModelViewer::OnToggleGuns(UI::CheckBox *w)
{
	if (!m_gunModel) {
		CreateTestResources();
	}

	if (!m_gunModel) {
		AddLog("test_gun.model not available");
		return false;
	}

	m_options.attachGuns = !m_options.attachGuns;
	SceneGraph::Model::TVecMT tags;
	m_model->FindTagsByStartOfName("tag_gun_", tags);
	if (tags.empty()) {
		AddLog("Missing tags \"tag_gun_XXX\" in model");
		return false;
	}
	if (m_options.attachGuns) {
		for (auto tag : tags) {
			tag->AddChild(new SceneGraph::ModelNode(m_gunModel.get()));
		}
	} else { //detach
		//we know there's nothing else
		for (auto tag : tags) {
			tag->RemoveChildAt(0);
		}
	}
	return true;
}

bool ModelViewer::OnRandomColor(UI::Widget*)
{
	if (!m_model || !m_model->SupportsPatterns()) return false;

	SceneGraph::ModelSkin skin;
	skin.SetRandomColors(m_rng);
	skin.Apply(m_model);

	// We need this flag setting so that we don't override what we're changing in OnModelColorsChanged
	m_settingColourSliders = true;
	const std::vector<Color> &colors = skin.GetColors();
	for(unsigned int i=0; i<3; i++) {
		for(unsigned int j=0; j<3; j++) {
			// use ToColor4f to get the colours in 0..1 range required
			if( colorSliders[(i*3)+j] )
				colorSliders[(i*3)+j]->SetValue(colors[i].ToColor4f()[j]);
		}
	}
	m_settingColourSliders = false;

	return true;
}

namespace
{
	using namespace SceneGraph;
	class FindStaticGeometryVisitor : public NodeVisitor {
	public:

		FindStaticGeometryVisitor() {};
		virtual void ApplyNode(Node &n)
		{
			if(n.GetTypeName() == "StaticGeometry") {
				auto sg = dynamic_cast<StaticGeometry*>(&n);
				m_results.push_back(sg);
			}

			n.Traverse(*this);
		}

		const std::vector<StaticGeometry*> &GetResults() { return m_results; }

	private:
		std::vector<StaticGeometry*> m_results;
	};
	class FindLODVisitor : public NodeVisitor {
	public:

		FindLODVisitor() {};
		virtual void ApplyNode(Node &n)
		{
			if (n.GetTypeName() == "LOD") {
				auto sg = dynamic_cast<LOD*>(&n);
				m_results.push_back(sg);
			}

			n.Traverse(*this);
		}

		const std::vector<LOD*> &GetResults() { return m_results; }

	private:
		std::vector<LOD*> m_results;
	};

	struct Sampler
	{
		//Cartesian direction
		vector3f direction;

		//Values of each SH function at this point
		std::vector<double> shValues;
	};

	bool GenerateSamples(const int numSamples, const int sqrtNumSamples, const int numBands, const int numFunctions, std::vector<Sampler> &samples)
	{
		//Create space for the SH values in each sample
		samples.resize(numSamples);
		for (int i = 0; i < numSamples; ++i)
		{
			samples[i].shValues.resize(numFunctions);
		}

		int index = 0;
		for (int i = 0; i < sqrtNumSamples; ++i)
		{
			for (int j = 0; j < sqrtNumSamples; ++j)
			{
				// Generate the position of this sample in [0, 1)x[0, 1)
				const float x = (i + ((float)rand() / RAND_MAX)) / sqrtNumSamples;
				const float y = (j + ((float)rand() / RAND_MAX)) / sqrtNumSamples;

				// Convert to spherical polars
				const float theta = 2.0*acos(sqrt(1.0 - x));
				const float phi = 2.0*M_PI*y;

				// Convert to cartesians
				samples[index].direction.x = float(sin(theta)*cos(phi)),
				samples[index].direction.y = float(sin(theta)*sin(phi)),
				samples[index].direction.z = float(cos(theta));

				//Compute SH coefficients for this sample
				for (int l = 0; l < numBands; ++l)
				{
					for (int m = -l; m <= l; ++m)
					{
						int index2 = l * (l + 1) + m;

						samples[index].shValues[index2] = SH::SampleSHFunc(l, m, theta, phi);
					}
				}

				++index;
			}
		}

		return true;
	}

	//double Light(const double theta, const double phi)
	//{
	//	return (theta < M_PI / 6) ? 1 : 0;
	//}
}

bool ModelViewer::OnCalculateSphericalHarmonics(UI::Widget*)
{
	m_uiMessageQueue.push_back(ECalculateSphericalHarmonics);
	return true;
}
#pragma optimize("",off)
void ModelViewer::CalculateSphericalHarmonics()
{
	Profiler::Timer timer;
	timer.Start();

	using SceneGraph::Node;
	using SceneGraph::Group;
	using SceneGraph::MatrixTransform;
	using SceneGraph::StaticGeometry;

	const int numBands = 3;
	const int numFunctions = numBands * numBands;
	const int sqrtNumSamples = 50;
	const int numSamples = sqrtNumSamples * sqrtNumSamples;
	std::vector<Sampler> samples;
	GenerateSamples(numSamples, sqrtNumSamples, numBands, numFunctions, samples);

	//This will find all StaticGeometry
	FindStaticGeometryVisitor StaticGeometryFinder;
	m_model->GetRoot()->Accept(StaticGeometryFinder);
	const std::vector<StaticGeometry*> &results = StaticGeometryFinder.GetResults();

	//This will find all StaticGeometry
	FindLODVisitor LODFinder;
	m_model->GetRoot()->Accept(LODFinder);
	const std::vector<LOD*> &LODresults = LODFinder.GetResults();
	for (auto lod : LODresults)
	{
		const Uint32 num = lod->GetNumChildren();
		for (Uint32 ichild = 0; ichild < num; ichild++) {
			auto child = lod->GetChildAt(ichild);
			Output("Type name: %s\n", child->GetTypeName());

			FindStaticGeometryVisitor StaticGeometryFinder;
			lod->Accept(StaticGeometryFinder);
			const std::vector<StaticGeometry*> &sgresults = StaticGeometryFinder.GetResults();
		}
	}

	//auto coll = m_model->CreateCollisionMesh();
	auto coll = m_model->GetCollisionMesh();

	//{
	//	m_model->
	//	//duplicate data again for geomtree...
	//	const size_t numVertices = m_vertices.size();
	//	const size_t numIndices = m_indices.size();
	//	const size_t numTris = numIndices / 3;
	//	std::vector<vector3f> vertices(numVertices);
	//	Uint32 *indices = new Uint32[numIndices];
	//	Uint32 *triFlags = new Uint32[numTris];

	//	m_totalTris += numTris;

	//	for (size_t i = 0; i < numVertices; i++)
	//		vertices[i] = m_vertices[i];

	//	for (size_t i = 0; i < numIndices; i++)
	//		indices[i] = m_indices[i];

	//	for (size_t i = 0; i < numTris; i++)
	//		triFlags[i] = m_flags[i];

	//	//create geomtree
	//	//takes ownership of data
	//	GeomTree *gt = new GeomTree(
	//		numVertices, numTris,
	//		vertices,
	//		indices, triFlags);
	//	m_collMesh->SetGeomTree(gt);
	//	m_collMesh->SetNumTriangles(m_totalTris);
	//}

	//Otherwise, regenerate the coefficients
	//Loop through the vertices and project the transfer function into SH space
	const size_t numObjects = results.size();
	for (int i = 0; i < numObjects; ++i)
	{
		const auto sg = results[i];
		const size_t numMeshes = sg->GetNumMeshes();
		for (size_t iMesh = 0; iMesh < numMeshes; iMesh++)
		{
			auto mesh = sg->GetMeshAt(iMesh);
			const Uint32 numVertices = mesh.vertexBuffer->GetSize();

			// get positions and normals from the mesh
			const auto& vbDesc = mesh.vertexBuffer->GetDesc();
			const Uint32 posOffset = vbDesc.GetOffset(Graphics::ATTRIB_POSITION);
			const Uint32 nrmOffset = vbDesc.GetOffset(Graphics::ATTRIB_NORMAL);
			const Uint32 uv0Offset = vbDesc.GetOffset(Graphics::ATTRIB_UV0);
			const Uint32 stride = vbDesc.stride;

			std::vector<vector3f> vertices(numVertices);
			std::vector<vector3f> normals(numVertices);
			std::vector<double> unshadowedCoeffs; unshadowedCoeffs.resize(numFunctions);
			std::vector<double> shadowedCoeffs; shadowedCoeffs.resize(numFunctions);
			Uint8 *vtxPtr = mesh.vertexBuffer->Map<Uint8>(Graphics::BUFFER_MAP_READ);
			for (Uint32 i = 0; i < vbDesc.numVertices; i++) {
				vertices.push_back(*reinterpret_cast<vector3f*>(vtxPtr + i * stride + posOffset));
				normals.push_back(*reinterpret_cast<vector3f*>(vtxPtr + i * stride + nrmOffset));
			}
			mesh.vertexBuffer->Unmap();

			for (int j = 0; j < numVertices; ++j)
			{
				const vector3f &position = vertices[j];
				const vector3f &normal = normals[j];

				//Clear the coefficients for this vertex
				for (int k = 0; k < numFunctions; ++k)
				{
					unshadowedCoeffs[k] = 0.0;
					shadowedCoeffs[k] = 0.0;
				}

				//Loop through samples
				for (int k = 0; k < numSamples; ++k)
				{
					//Calculate cosine term for this sample
					const double dot = samples[k].direction.Dot(normal);

					//Clamp to [0, 1]
					if (dot > 0.0)
					{
						//See if the ray is blocked by any object
						isect_t isect;
						isect.dist = float(1000.0f);
						isect.triIdx = -1;
						coll->GetGeomTree()->TraceRay(position + 2 * std::numeric_limits<float>::epsilon() * normal, samples[k].direction, &isect);
						bool rayBlocked = (isect.triIdx != -1);

						//Add the contribution of this sample to the coefficients
						for (int l = 0; l < numFunctions; ++l)
						{
							double contribution = dot * samples[k].shValues[l];

							unshadowedCoeffs[l] += contribution;

							if (!rayBlocked)
								shadowedCoeffs[l] += contribution;
						}
					}
				}

				//Rescale the coefficients
				for (int k = 0; k < numFunctions; ++k)
				{
					unshadowedCoeffs[k] *= 4 * M_PI / numSamples;
					shadowedCoeffs[k] *= 4 * M_PI / numSamples;
				}
			}
		}
	}
	timer.Stop();
	Output("\n\nCalculateSphericalHarmonics took: %lf milliseconds\n", timer.millicycles());
}

void ModelViewer::UpdateShield()
{
	if (m_shieldIsHit) {
		m_shieldHitPan += 0.05f;
	}
	if (m_shieldHitPan > 0.34f) {
		m_shieldHitPan = -1.48f;
		m_shieldIsHit = false;
	}
}

bool ModelViewer::OnHitIt(UI::Widget*)
{
	HitImpl();
	return true;
}

void ModelViewer::HitImpl()
{
	if(m_model) {
		assert(m_shields.get());
		// pick a point on the shield to serve as the point of impact.
		SceneGraph::StaticGeometry* sg = m_shields->GetFirstShieldMesh();
		if(sg) {
			SceneGraph::StaticGeometry::Mesh &mesh = sg->GetMeshAt(0);

			// Please don't do this in game, no speed guarantee
			const Uint32 posOffs = mesh.vertexBuffer->GetDesc().GetOffset(Graphics::ATTRIB_POSITION);
			const Uint32 stride  = mesh.vertexBuffer->GetDesc().stride;
			const Uint32 vtxIdx = m_rng.Int32() % mesh.vertexBuffer->GetSize();

			const Uint8 *vtxPtr = mesh.vertexBuffer->Map<Uint8>(Graphics::BUFFER_MAP_READ);
			const vector3f pos = *reinterpret_cast<const vector3f*>(vtxPtr + vtxIdx * stride + posOffs);
			mesh.vertexBuffer->Unmap();
			m_shields->AddHit(vector3d(pos));
		}
	}
	m_shieldHitPan = -1.48f;
	m_shieldIsHit = true;
}

void ModelViewer::AddLog(const std::string &line)
{
	m_log->AppendText(line+"\n");
	m_logScroller->SetScrollPosition(1.0f);
	Output("%s\n", line.c_str());
}

void ModelViewer::ChangeCameraPreset(SDL_Keycode key, SDL_Keymod mod)
{
	if (!m_model) return;

	// Like Blender, but a bit different because we like that
	// 1 - front (+ctrl back)
	// 7 - top (+ctrl bottom)
	// 3 - left (+ctrl right)
	// 2,4,6,8 incrementally rotate

	const bool invert = mod & KMOD_CTRL;

	switch (key)
	{
	case SDLK_KP_7: case SDLK_u:
		m_rotX = invert ? -90.f : 90.f;
		m_rotY = 0.f;
		AddLog(invert ? "Bottom view" : "Top view");
		break;
	case SDLK_KP_3: case SDLK_PERIOD:
		m_rotX = 0.f;
		m_rotY = invert ? -90.f : 90.f;
		AddLog(invert ? "Right view" : "Left view");
		break;
	case SDLK_KP_1: case SDLK_m:
		m_rotX = 0.f;
		m_rotY = invert ? 0.f : 180.f;
		AddLog(invert ? "Rear view" : "Front view");
		break;
	case SDLK_KP_4: case SDLK_j:
		m_rotY += 15.f;
		break;
	case SDLK_KP_6: case SDLK_l:
		m_rotY -= 15.f;
		break;
	case SDLK_KP_2: case SDLK_COMMA:
		m_rotX += 15.f;
		break;
	case SDLK_KP_8: case SDLK_i:
		m_rotX -= 15.f;
		break;
	default:
		break;
		//no others yet
	}
}

void ModelViewer::ToggleViewControlMode()
{
	m_options.mouselookEnabled = !m_options.mouselookEnabled;
	m_renderer->SetGrab(m_options.mouselookEnabled);

	if (m_options.mouselookEnabled) {
		m_viewRot = matrix3x3f::RotateY(DEG2RAD(m_rotY)) * matrix3x3f::RotateX(DEG2RAD(Clamp(m_rotX, -90.0f, 90.0f)));
		m_viewPos = zoom_distance(m_baseDistance, m_zoom) * m_viewRot.VectorZ();
	} else {
		// XXX re-initialise the turntable style view position from the current mouselook view
		ResetCamera();
	}
}

void ModelViewer::ClearLog()
{
	m_log->SetText("");
}

void ModelViewer::ClearModel()
{
	delete m_model; m_model = 0;
	m_gunModel.reset();
	m_scaleModel.reset();

	m_options.mouselookEnabled = false;
	m_renderer->SetGrab(false);
	m_viewPos = vector3f(0.0f, 0.0f, 10.0f);
	ResetCamera();
}

void ModelViewer::CreateTestResources()
{
	//load gun model for attachment test
	//landingpad model for scale test
	SceneGraph::Loader loader(m_renderer);
	try {
		SceneGraph::Model *m = loader.LoadModel("test_gun");
		m_gunModel.reset(m);

		m = loader.LoadModel("scale");
		m_scaleModel.reset(m);
	} catch (SceneGraph::LoadingError &) {
		AddLog("Could not load test_gun or scale model");
	}
}

void ModelViewer::DrawBackground()
{
	m_renderer->SetOrthographicProjection(0.f, 1.f, 0.f, 1.f, -1.f, 1.f);
	m_renderer->SetTransform(matrix4x4f::Identity());

	if(!m_bgBuffer.Valid())
	{
		const Color top = Color::BLACK;
		const Color bottom = Color(77, 77, 77);
		Graphics::VertexArray bgArr(Graphics::ATTRIB_POSITION | Graphics::ATTRIB_DIFFUSE, 6);
		// triangle 1
		bgArr.Add(vector3f(0.f, 0.f, 0.f), bottom);
		bgArr.Add(vector3f(1.f, 0.f, 0.f), bottom);
		bgArr.Add(vector3f(1.f, 1.f, 0.f), top);
		// triangle 2
		bgArr.Add(vector3f(0.f, 0.f, 0.f), bottom);
		bgArr.Add(vector3f(1.f, 1.f, 0.f), top);
		bgArr.Add(vector3f(0.f, 1.f, 0.f), top);

		Graphics::VertexBufferDesc vbd;
		vbd.attrib[0].semantic	= Graphics::ATTRIB_POSITION;
		vbd.attrib[0].format	= Graphics::ATTRIB_FORMAT_FLOAT3;
		vbd.attrib[1].semantic	= Graphics::ATTRIB_DIFFUSE;
		vbd.attrib[1].format	= Graphics::ATTRIB_FORMAT_UBYTE4;
		vbd.numVertices = 6;
		vbd.usage = Graphics::BUFFER_USAGE_STATIC;

		// VertexBuffer
		m_bgBuffer.Reset( m_renderer->CreateVertexBuffer(vbd) );
		m_bgBuffer->Populate(bgArr);
	}

	m_renderer->DrawBuffer(m_bgBuffer.Get(), m_bgState, Graphics::vtxColorMaterial, Graphics::TRIANGLES);
}

//Draw grid and axes
void ModelViewer::DrawGrid(const matrix4x4f &trans, float radius)
{
	assert(m_options.showGrid);

	const float dist = zoom_distance(m_baseDistance, m_zoom);

	const float max = std::min(powf(10, ceilf(log10f(dist))), ceilf(radius/m_options.gridInterval)*m_options.gridInterval);

	static std::vector<vector3f> points;
	points.clear();

	for (float x = -max; x <= max; x += m_options.gridInterval) {
		points.push_back(vector3f(x,0,-max));
		points.push_back(vector3f(x,0,max));
		points.push_back(vector3f(0,x,-max));
		points.push_back(vector3f(0,x,max));

		points.push_back(vector3f(x,-max,0));
		points.push_back(vector3f(x,max,0));
		points.push_back(vector3f(0,-max,x));
		points.push_back(vector3f(0,max,x));

		points.push_back(vector3f(-max,x,0));
		points.push_back(vector3f(max,x,0));
		points.push_back(vector3f(-max,0,x));
		points.push_back(vector3f(max,0,x));
	}

	m_renderer->SetTransform(trans);
	m_gridLines.SetData(points.size(), &points[0], Color(128, 128, 128));
	m_gridLines.Draw(m_renderer, m_bgState);

	// industry-standard red/green/blue XYZ axis indiactor
	m_renderer->SetTransform(trans * matrix4x4f::ScaleMatrix(radius));
	Graphics::Drawables::GetAxes3DDrawable(m_renderer)->Draw(m_renderer);
}

void ModelViewer::DrawModel(const matrix4x4f &mv)
{
	assert(m_model);

	m_model->UpdateAnimations();

	m_model->SetDebugFlags(
		(m_options.showAabb            ? SceneGraph::Model::DEBUG_BBOX      : 0x0) |
		(m_options.showCollMesh        ? SceneGraph::Model::DEBUG_COLLMESH  : 0x0) |
		(m_options.showTags            ? SceneGraph::Model::DEBUG_TAGS      : 0x0) |
		(m_options.showDockingLocators ? SceneGraph::Model::DEBUG_DOCKING   : 0x0) |
		(m_options.wireframe           ? SceneGraph::Model::DEBUG_WIREFRAME : 0x0)
	);

	m_model->Render(mv);
	m_navLights->Render(m_renderer);
}

void ModelViewer::MainLoop()
{
	double lastTime = SDL_GetTicks() * 0.001;
	while (!m_done)
	{
		const double ticks = SDL_GetTicks() * 0.001;
		m_frameTime = (ticks - lastTime);
		lastTime = ticks;

		// logic update
		PollEvents();

		m_renderer->ClearScreen();
		UpdateLights();
		UpdateCamera();
		UpdateShield();

		// render the gradient backdrop
		DrawBackground();

		//update animations, draw model etc.
		if (m_model) {
			m_navLights->Update(m_frameTime);
			m_shields->SetEnabled(m_options.showShields || m_shieldIsHit);

			//Calculate the impact's radius dependant on time
			const float dif1 = 0.34 - (-1.48f);
			const float dif2 = m_shieldHitPan - (-1.48f);
			//Range from start (0.0) to end (1.0)
			const float dif = dif2 / (dif1 * 1.0f);

			m_shields->Update(m_options.showShields ? 1.0f : (1.0f - dif), 1.0f);

			// setup rendering
			if (!m_options.orthoView) {
                m_renderer->SetPerspectiveProjection(85, Graphics::GetScreenWidth()/float(Graphics::GetScreenHeight()), 0.1f, 100000.f);
			} else {
                /* TODO: Zoom in ortho mode seems don't work as in perspective mode,
                 / I change "screen dimensions" to avoid the problem.
                 / However the zoom needs more care
                */
                if (m_zoom<=0.0) m_zoom = 0.01;
                float screenW = Graphics::GetScreenWidth()*m_zoom/10;
                float screenH = Graphics::GetScreenHeight()*m_zoom/10;
                matrix4x4f orthoMat = matrix4x4f::OrthoFrustum(-screenW,screenW, -screenH, screenH, 0.1f, 100000.0f);
                orthoMat.ClearToRotOnly();
                m_renderer->SetProjection(orthoMat);
			}

			m_renderer->SetTransform(matrix4x4f::Identity());

			// calc camera info
			matrix4x4f mv;
			float zd=0;
			if (m_options.mouselookEnabled) {
				mv = m_viewRot.Transpose() * matrix4x4f::Translation(-m_viewPos);
			} else {
				m_rotX = Clamp(m_rotX, -90.0f, 90.0f);
				matrix4x4f rot = matrix4x4f::Identity();
				rot.RotateX(DEG2RAD(-m_rotX));
				rot.RotateY(DEG2RAD(-m_rotY));
				if (m_options.orthoView) zd = -m_baseDistance;
				else zd = -zoom_distance(m_baseDistance, m_zoom);
				mv = matrix4x4f::Translation(0.0f, 0.0f, zd) * rot;
			}

			// draw the model itself
			DrawModel(mv);

			// helper rendering
			if (m_options.showLandingPad) {
				if (!m_scaleModel) CreateTestResources();
				m_scaleModel->Render(mv * matrix4x4f::Translation(0.f, m_landingMinOffset, 0.f));
			}

			if (m_options.showGrid) {
				DrawGrid(mv, m_model->GetDrawClipRadius());
			}
		}

		m_ui->Update();
		if (m_options.showUI && !m_screenshotQueued) {
			m_ui->Draw();
		}
		if (m_screenshotQueued) {
			m_screenshotQueued = false;
			Screenshot();
		}

		// end scene
		m_renderer->SwapBuffers();

		// if we've requested a different model then switch too it
		if(!m_requestedModelName.empty()) {
			SetModel(m_requestedModelName);
			m_requestedModelName.clear();
			ResetCamera();
		}
	}
}

void ModelViewer::OnAnimChanged(unsigned int, const std::string &name)
{
	m_currentAnimation = 0;
	// Find the animation matching the name (could also store the anims in a map
	// when the animationSelector is filled)
	if (!name.empty()) {
		const std::vector<SceneGraph::Animation*> &anims = m_model->GetAnimations();
		for (std::vector<SceneGraph::Animation*>::const_iterator anim = anims.begin(); anim != anims.end(); ++anim) {
			if ((*anim)->GetName() == name)
				m_currentAnimation = (*anim);
		}
	}
	if (m_currentAnimation)
		animSlider->SetValue(m_currentAnimation->GetProgress());
	else
		animSlider->SetValue(0.0);
}

void ModelViewer::OnAnimSliderChanged(float value)
{
	if (m_currentAnimation)
		m_currentAnimation->SetProgress(value);
	animValue->SetText(stringf("%0{f.2}", value));
}

void ModelViewer::OnDecalChanged(unsigned int index, const std::string &texname)
{
	if (!m_model) return;

	m_decalTexture = Graphics::TextureBuilder::Decal(stringf("textures/decals/%0.dds", texname)).GetOrCreateTexture(m_renderer, "decal");

	m_model->SetDecalTexture(m_decalTexture, 0);
	m_model->SetDecalTexture(m_decalTexture, 1);
	m_model->SetDecalTexture(m_decalTexture, 2);
	m_model->SetDecalTexture(m_decalTexture, 3);
}

void ModelViewer::OnLightPresetChanged(unsigned int index, const std::string&)
{
	m_options.lightPreset = std::min<unsigned int>(index, 3);
}

void ModelViewer::OnModelColorsChanged(float)
{
	if (!m_model || m_settingColourSliders) return;

	//don't care about the float. Fetch values from all sliders.
	std::vector<Color> colors;
	colors.push_back(get_slider_color(colorSliders[0], colorSliders[1], colorSliders[2]));
	colors.push_back(get_slider_color(colorSliders[3], colorSliders[4], colorSliders[5]));
	colors.push_back(get_slider_color(colorSliders[6], colorSliders[7], colorSliders[8]));
	m_model->SetColors(colors);
}

void ModelViewer::OnPatternChanged(unsigned int index, const std::string &value)
{
	if (!m_model) return;
	assert(index < m_model->GetPatterns().size());
	m_model->SetPattern(index);
}

void ModelViewer::OnThrustChanged(float)
{
	vector3f linthrust;
	vector3f angthrust;

	linthrust.x = get_thrust(thrustSliders[0]);
	linthrust.y = get_thrust(thrustSliders[1]);
	linthrust.z = get_thrust(thrustSliders[2]);

	// angthrusts are negated in ship.cpp for some reason
	angthrust.x = -get_thrust(thrustSliders[3]);
	angthrust.y = -get_thrust(thrustSliders[4]);
	angthrust.z = -get_thrust(thrustSliders[5]);

	m_model->SetThrust(linthrust, angthrust);
}

void ModelViewer::PollEvents()
{
	/*
	 * Special butans
	 *
	 * Space: reset camera
	 * Keyboard: rotate view
	 * plus/minus: zoom view
	 * Shift: zoom faster
	 * printscr - screenshots
	 * tab - toggle ui (always invisible on screenshots)
	 * g - grid
	 * o - switch orthographic<->perspective
	 *
	 */
	m_mouseMotion[0] = m_mouseMotion[1] = 0;
	m_mouseWheelUp = m_mouseWheelDown = false;

	SDL_Event event;
	while (SDL_PollEvent(&event)) {
		//ui gets all events
		if (m_options.showUI && m_ui->DispatchSDLEvent(event))
			continue;

		switch (event.type)
		{
		case SDL_QUIT:
			m_done = true;
			break;
		case SDL_MOUSEMOTION:
			m_mouseMotion[0] += event.motion.xrel;
			m_mouseMotion[1] += event.motion.yrel;
			break;
		case SDL_MOUSEBUTTONDOWN:
			m_mouseButton[event.button.button] = true;
			break;
		case SDL_MOUSEBUTTONUP:
			m_mouseButton[event.button.button] = false;
			break;
		case SDL_MOUSEWHEEL:
			if (event.wheel.y > 0) m_mouseWheelUp = true;
			if (event.wheel.y < 0) m_mouseWheelDown = true;
			break;
		case SDL_KEYDOWN:
			switch (event.key.keysym.sym)
			{
			case SDLK_ESCAPE:
				if (m_model) {
					ClearModel();
					onModelChanged.emit();
					PopulateFilePicker();
				} else {
					m_done = true;
				}
				break;
			case SDLK_SPACE:
				ResetCamera();
				ResetThrusters();
				break;
			case SDLK_TAB:
				m_options.showUI = !m_options.showUI;
				break;
			case SDLK_t:
				m_options.showTags = !m_options.showTags;
				break;
			case SDLK_PRINTSCREEN:
				m_screenshotQueued = true;
				break;
			case SDLK_g:
				OnToggleGrid(0);
				break;
            case SDLK_o:
                m_options.orthoView = !m_options.orthoView;
                break;
			case SDLK_z:
				m_options.wireframe = !m_options.wireframe;
				break;
			case SDLK_f:
				ToggleViewControlMode();
				break;
			case SDLK_F6:
				SaveModelToBinary();
				break;
			case SDLK_F11:
				if (event.key.keysym.mod & KMOD_SHIFT)
					m_renderer->ReloadShaders();
				break;
			case SDLK_KP_1: case SDLK_m:
			case SDLK_KP_2: case SDLK_COMMA:
			case SDLK_KP_3: case SDLK_PERIOD:
			case SDLK_KP_4: case SDLK_j:
			case SDLK_KP_6: case SDLK_l:
			case SDLK_KP_7: case SDLK_u:
			case SDLK_KP_8: case SDLK_i:
				ChangeCameraPreset(event.key.keysym.sym, SDL_Keymod(event.key.keysym.mod));
				break;
			case SDLK_p: //landing pad test
				m_options.showLandingPad = !m_options.showLandingPad;
				AddLog(stringf("Scale/landing pad test %0", m_options.showLandingPad ? "on" : "off"));
				break;
			case SDLK_r: //random colors, eastereggish
				for(unsigned int i=0; i<3*3; i++) {
					if (colorSliders[i])
						colorSliders[i]->SetValue(m_rng.Double());
				}
				break;
			default:
				break; //shuts up -Wswitch
			} //keysym switch
			m_keyStates[event.key.keysym.sym] = true;
			break;
		case SDL_KEYUP:
			m_keyStates[event.key.keysym.sym] = false;
			break;
		default:
			break;
		}
	}

	// handle all UI messages we've queued
	for (auto message : m_uiMessageQueue)
	{
		switch (message)
		{
		case ECalculateSphericalHarmonics:
			CalculateSphericalHarmonics();
			break;
		default:
			break;
		}
	}
	m_uiMessageQueue.clear();
}

static void collect_models(std::vector<std::string> &list)
{
	const std::string basepath("models");
	FileSystem::FileSource &fileSource = FileSystem::gameDataFiles;
	for (FileSystem::FileEnumerator files(fileSource, basepath, FileSystem::FileEnumerator::Recurse); !files.Finished(); files.Next())
	{
		const FileSystem::FileInfo &info = files.Current();
		const std::string &fpath = info.GetPath();

		//check it's the expected type
		if (info.IsFile()) {
			if (ends_with_ci(fpath, ".model"))
				list.push_back(info.GetName().substr(0, info.GetName().size()-6));
			else if (ends_with_ci(fpath, ".sgm"))
				list.push_back(info.GetName());
		}
	}
}

void ModelViewer::PopulateFilePicker()
{
	m_fileList->Clear();

	std::vector<std::string> models;
	collect_models(models);

	for (const auto& it : models)
		m_fileList->AddOption(it);
}

void ModelViewer::ResetCamera()
{
	m_baseDistance = m_model ? m_model->GetDrawClipRadius() * 1.5f : 100.f;
	m_rotX = m_rotY = 0.f;
	m_zoom = 0.f;
}

void ModelViewer::ResetThrusters()
{
	if (thrustSliders[0] == 0) return;

	for (unsigned int i=0; i<6; i++) {
		thrustSliders[i]->SetValue(0.5f);
	}
}

void ModelViewer::Screenshot()
{
	char buf[256];
	const time_t t = time(0);
	const struct tm *_tm = localtime(&t);
	strftime(buf, sizeof(buf), "modelviewer-%Y%m%d-%H%M%S.png", _tm);
	Graphics::ScreendumpState sd;
	m_renderer->Screendump(sd);
	write_screenshot(sd, buf);
	AddLog(stringf("Screenshot %0 saved", buf));
}

void ModelViewer::SaveModelToBinary()
{
	if (!m_model)
		return AddLog("No current model to binarize");

	//load the current model in a pristine state (no navlights, shields...)
	//and then save it into binary

	std::unique_ptr<SceneGraph::Model> model;
	try {
		SceneGraph::Loader ld(m_renderer);
		model.reset(ld.LoadModel(m_modelName));
	} catch (...) {
		//minimal error handling, this is not expected to happen since we got this far.
		AddLog("Could not load model");
		return;
	}

	try {
		SceneGraph::BinaryConverter bc(m_renderer);
		bc.Save(m_modelName, model.get());
		AddLog("Saved binary model file");
	} catch (const CouldNotOpenFileException&) {
		AddLog("Could not open file or directory for writing");
	} catch (const CouldNotWriteToFileException&) {
		AddLog("Error while writing to file");
	}
}

void ModelViewer::SetModel(const std::string &filename)
{
	AddLog(stringf("Loading model %0...", filename));

	//this is necessary to reload textures
	m_renderer->RemoveAllCachedTextures();

	ClearModel();

	try {
		if (ends_with_ci(filename, ".sgm")) {
			//binary loader expects extension-less name. Might want to change this.
			m_modelName = filename.substr(0, filename.size()-4);
			SceneGraph::BinaryConverter bc(m_renderer);
			m_model = bc.Load(m_modelName);
		} else {
			m_modelName = filename;
			SceneGraph::Loader loader(m_renderer, true);
			m_model = loader.LoadModel(filename);

			//dump warnings
			for (std::vector<std::string>::const_iterator it = loader.GetLogMessages().begin();
				it != loader.GetLogMessages().end(); ++it)
			{
				AddLog(*it);
				Output("%s\n", (*it).c_str());
			}
		}

		Shields::ReparentShieldNodes(m_model);

		//set decal textures, max 4 supported.
		//Identical texture at the moment
		OnDecalChanged(0, "pioneer");
		Output("\n\n");

		SceneGraph::DumpVisitor d(m_model);
		m_model->GetRoot()->Accept(d);
		AddLog(d.GetModelStatistics());

		// If we've got the tag_landing set then use it for an offset otherwise grab the AABB
		const SceneGraph::MatrixTransform *mt = m_model->FindTagByName("tag_landing");
		if (mt)
			m_landingMinOffset = mt->GetTransform().GetTranslate().y;
		else if (m_model->GetCollisionMesh())
			m_landingMinOffset = m_model->GetCollisionMesh()->GetAabb().min.y;
		else
			m_landingMinOffset = 0.0f;

		//note: stations won't demonstrate full docking light logic in MV
		m_navLights.reset(new NavLights(m_model));
		m_navLights->SetEnabled(true);

		m_shields.reset(new Shields(m_model));
	} catch (SceneGraph::LoadingError &err) {
		// report the error and show model picker.
		m_model = 0;
		AddLog(stringf("Could not load model %0: %1", filename, err.what()));
	}

	onModelChanged.emit();
}

void ModelViewer::SetupFilePicker()
{
	UI::Context *c = m_ui.Get();

	m_fileList = c->List();
	UI::Button *quitButton = c->Button();
	UI::Button *loadButton = c->Button();
	quitButton->SetInnerWidget(c->Label("Quit"));
	loadButton->SetInnerWidget(c->Label("Load"));

	PopulateFilePicker();

	UI::Widget *fp =
	c->Grid(UI::CellSpec(1,3,1), UI::CellSpec(1,3,1))
		->SetCell(1,1,
			c->VBox(10)
				->PackEnd(c->Label("Select a model"))
				->PackEnd(c->Expand(UI::Expand::BOTH)->SetInnerWidget(c->Scroller()->SetInnerWidget(m_fileList)))
				->PackEnd(c->Grid(2,1)->SetRow(0, UI::WidgetSet(
					c->Align(UI::Align::LEFT)->SetInnerWidget(loadButton),
					c->Align(UI::Align::RIGHT)->SetInnerWidget(quitButton)
				)))
		);

	m_logScroller->Layout(); //issues without this
	c->GetTopLayer()->SetInnerWidget(c->Grid(2,1)
		->SetRow(0, UI::WidgetSet(fp, m_logScroller.Get()))
	);

	c->Layout();
	m_logScroller->SetScrollPosition(1.f);

	loadButton->onClick.connect(sigc::bind(sigc::mem_fun(*this, &ModelViewer::OnPickModel), m_fileList));
	quitButton->onClick.connect(sigc::mem_fun(*this, &ModelViewer::OnQuit));
}

void ModelViewer::SetupUI()
{
	UI::Context *c = m_ui.Get();
	c->SetFont(UI::Widget::FONT_XSMALL);

	for (unsigned int i=0; i<9; i++)
		colorSliders[i] = 0;
	for (unsigned int i=0; i<6; i++)
		thrustSliders[i] = 0;

	animSlider = 0;
	animValue = 0;

	if (!m_model)
		return SetupFilePicker();

	const int spacing = 5;

	UI::SmallButton *reloadButton = nullptr;
	UI::SmallButton *toggleGridButton = nullptr;
	UI::SmallButton *hitItButton = nullptr;
	UI::SmallButton *randomColours = nullptr;
	UI::SmallButton *calculateSphericalHarmonics = nullptr;
	UI::CheckBox *collMeshCheck = nullptr;
	UI::CheckBox *showShieldsCheck = nullptr;
	UI::CheckBox *gunsCheck = nullptr;

	UI::VBox* outerBox = c->VBox();

	UI::VBox* mainBox = c->VBox(5);
	UI::VBox* bottomBox = c->VBox(5);

	UI::HBox* sliderBox = c->HBox();
	bottomBox->PackEnd(sliderBox);

	outerBox->PackEnd(UI::WidgetSet(
		c->Expand()->SetInnerWidget(c->Grid(UI::CellSpec(0.30f,0.8f,0.35f),1)
			->SetColumn(0, mainBox)
			->SetColumn(2, m_logScroller.Get())
		),
		bottomBox
	));

	c->GetTopLayer()->SetInnerWidget(c->Margin(spacing)->SetInnerWidget(outerBox));

	//model name + reload button: visible even if loading failed
	mainBox->PackEnd(nameLabel = c->Label(m_modelName));
	nameLabel->SetFont(UI::Widget::FONT_NORMAL);
	add_pair(c, mainBox, reloadButton = c->SmallButton(), "Reload model");
	reloadButton->onClick.connect(sigc::bind(sigc::mem_fun(*this, &ModelViewer::OnReloadModel), reloadButton));

	if (m_model == 0) {
		c->Layout();
		return;
	}

	add_pair(c, mainBox, toggleGridButton = c->SmallButton(), "Grid mode");
	add_pair(c, mainBox, collMeshCheck = c->CheckBox(), "Collision mesh");
	// not everything has a shield
	if( m_shields.get() && m_shields->GetFirstShieldMesh() ) {
		add_pair(c, mainBox, showShieldsCheck = c->CheckBox(), "Show Shields");
		add_pair(c, mainBox, hitItButton = c->SmallButton(), "Hit it!");
		hitItButton->onClick.connect(sigc::bind(sigc::mem_fun(*this, &ModelViewer::OnHitIt), hitItButton));
	}

	
	add_pair(c, mainBox, calculateSphericalHarmonics = c->SmallButton(), "Calculate SH");
	calculateSphericalHarmonics->onClick.connect(sigc::bind(sigc::mem_fun(*this, &ModelViewer::OnCalculateSphericalHarmonics), calculateSphericalHarmonics));

	add_pair(c, mainBox, randomColours = c->SmallButton(), "Random Colours");
	randomColours->onClick.connect(sigc::bind(sigc::mem_fun(*this, &ModelViewer::OnRandomColor), randomColours));


	//pattern selector
	if (m_model->SupportsPatterns()) {
		mainBox->PackEnd(c->Label("Pattern:"));
		mainBox->PackEnd(patternSelector = c->DropDown()->AddOption("Default"));

		sliderBox->PackEnd(
			c->Grid(3,4)
				->SetColumn(0, UI::WidgetSet(
					c->Label("Color 1"),
					c->HBox(spacing)->PackEnd(c->Label("R"))->PackEnd(colorSliders[0] = c->HSlider()),
					c->HBox(spacing)->PackEnd(c->Label("G"))->PackEnd(colorSliders[1] = c->HSlider()),
					c->HBox(spacing)->PackEnd(c->Label("B"))->PackEnd(colorSliders[2] = c->HSlider())
				))
				->SetColumn(1, UI::WidgetSet(
					c->Label("Color 2"),
					c->HBox(spacing)->PackEnd(c->Label("R"))->PackEnd(colorSliders[3] = c->HSlider()),
					c->HBox(spacing)->PackEnd(c->Label("G"))->PackEnd(colorSliders[4] = c->HSlider()),
					c->HBox(spacing)->PackEnd(c->Label("B"))->PackEnd(colorSliders[5] = c->HSlider())
				))
				->SetColumn(2, UI::WidgetSet(
					c->Label("Color 3"),
					c->HBox(spacing)->PackEnd(c->Label("R"))->PackEnd(colorSliders[6] = c->HSlider()),
					c->HBox(spacing)->PackEnd(c->Label("G"))->PackEnd(colorSliders[7] = c->HSlider()),
					c->HBox(spacing)->PackEnd(c->Label("B"))->PackEnd(colorSliders[8] = c->HSlider())
				))
		);

		//connect slider signals, set initial values (RGB)
		const float values[] = {
			1.f, 0.f, 0.f,
			0.f, 1.f, 0.f,
			0.f, 0.f, 1.f
		};
		for(unsigned int i=0; i<3*3; i++) {
			colorSliders[i]->SetValue(values[i]);
			colorSliders[i]->onValueChanged.connect(sigc::mem_fun(*this, &ModelViewer::OnModelColorsChanged));
		}
		//// slidems end

		patternSelector->onOptionSelected.connect(sigc::mem_fun(*this, &ModelViewer::OnPatternChanged));

		UpdatePatternList();
	}

	//decal selector
	//models support up to 4 but 1 is enough here
	if (m_model->SupportsDecals()) {
		mainBox->PackEnd(c->Label("Decal:"));
		mainBox->PackEnd(decalSelector = c->DropDown());

		decalSelector->onOptionSelected.connect(sigc::mem_fun(*this, &ModelViewer::OnDecalChanged));

		std::vector<std::string> decals;
		collect_decals(decals);

		for (std::vector<std::string>::const_iterator it = decals.begin(); it != decals.end(); ++it) {
			decalSelector->AddOption(*it);
		}
		if (decals.size() > 0)
			decalSelector->SetSelectedOption("pioneer");
	}

	//light dropdown
	UI::DropDown *lightSelector;
	mainBox->PackEnd(c->Label("Lights:"));
	mainBox->PackEnd(
		lightSelector = c->DropDown()
			->AddOption("1  Front white")
			->AddOption("2  Two-point")
			->AddOption("3  Backlight")
			//->AddOption("4  Nuts")
	);
	lightSelector->SetSelectedOption("1  Front white");
	m_options.lightPreset = 0;

	add_pair(c, mainBox, gunsCheck = c->CheckBox(), "Attach guns");

	//Animation controls
	if (!m_model->GetAnimations().empty()) {
		//UI::Button *playBtn;
		//UI::Button *revBtn;
		//UI::Button *stopBtn;
		UI::Box *animBox;
		mainBox->PackEnd(animBox = c->VBox(spacing));
		animBox->PackEnd(m_ui->Label("Animation:"));
		animBox->PackEnd(animSelector = m_ui->DropDown()->AddOption("None"));
		//add_pair(m_ui, animBox, playBtn = m_ui->Button(), "Play/Pause");
		//add_pair(m_ui, animBox, revBtn = m_ui->Button(), "Play reverse");
		//add_pair(m_ui, animBox, stopBtn = m_ui->Button(), "Stop");

		bottomBox->PackStart(c->HBox(10)->PackEnd(UI::WidgetSet(c->Label("Animation:"), animSlider = c->HSlider(), animValue = c->Label("0.0"))));
		animValue->SetFont(UI::Widget::FONT_NORMAL);

		//playBtn->onClick.connect(sigc::bind(sigc::mem_fun(*this, &ModelViewer::OnAnimPlay), playBtn, false));
		//revBtn->onClick.connect(sigc::bind(sigc::mem_fun(*this, &ModelViewer::OnAnimPlay), revBtn, true));
		//stopBtn->onClick.connect(sigc::bind(sigc::mem_fun(*this, &ModelViewer::OnAnimStop), stopBtn));
		animSlider->onValueChanged.connect(sigc::mem_fun(*this, &ModelViewer::OnAnimSliderChanged));
		animSelector->onOptionSelected.connect(sigc::mem_fun(*this, &ModelViewer::OnAnimChanged));

		//update anims from model
		UpdateAnimList();
	}

	//// Thrust sliders
	bool supportsThrusters = false;
	{
		SceneGraph::FindNodeVisitor fivi(SceneGraph::FindNodeVisitor::MATCH_NAME_STARTSWITH, "thruster_");
		m_model->GetRoot()->Accept(fivi);
		supportsThrusters = !fivi.GetResults().empty();
	}
	if (supportsThrusters) {
		sliderBox->PackStart(
			c->Grid(2,4)
				->SetColumn(0, UI::WidgetSet(
					// Column 1, Linear thrust sliders
					c->Label("Linear"),
					c->HBox(spacing)->PackEnd(c->Label("X"))->PackEnd(thrustSliders[0] = c->HSlider()),
					c->HBox(spacing)->PackEnd(c->Label("Y"))->PackEnd(thrustSliders[1] = c->HSlider()),
					c->HBox(spacing)->PackEnd(c->Label("Z"))->PackEnd(thrustSliders[2] = c->HSlider())
				))
				->SetColumn(1, UI::WidgetSet(
					//Column 2, Angular thrust sliders
					c->Label("Angular"),
					c->HBox(spacing)->PackEnd(c->Label("Pitch"))->PackEnd(thrustSliders[3] = c->HSlider()),
					c->HBox(spacing)->PackEnd(c->Label("Yaw"))->PackEnd(thrustSliders[4] = c->HSlider()),
					c->HBox(spacing)->PackEnd(c->Label("Roll"))->PackEnd(thrustSliders[5] = c->HSlider())
				))
		);
		for(unsigned int i=0; i<2*3; i++) {
			thrustSliders[i]->SetValue(0.5f);
			thrustSliders[i]->onValueChanged.connect(sigc::mem_fun(*this, &ModelViewer::OnThrustChanged));
		}
		////thruster sliders end
	}

	c->Layout();

	//event handlers
	collMeshCheck->onClick.connect(sigc::bind(sigc::mem_fun(*this, &ModelViewer::OnToggleCollMesh), collMeshCheck));
	if( m_shields.get() && showShieldsCheck ) { showShieldsCheck->onClick.connect(sigc::bind(sigc::mem_fun(*this, &ModelViewer::OnToggleShowShields), showShieldsCheck)); }
	gunsCheck->onClick.connect(sigc::bind(sigc::mem_fun(*this, &ModelViewer::OnToggleGuns), gunsCheck));
	lightSelector->onOptionSelected.connect(sigc::mem_fun(*this, &ModelViewer::OnLightPresetChanged));
	toggleGridButton->onClick.connect(sigc::bind(sigc::mem_fun(*this, &ModelViewer::OnToggleGrid), toggleGridButton));
}

void ModelViewer::UpdateAnimList()
{
	animSelector->Clear();
	if (m_model) {
		const std::vector<SceneGraph::Animation*> &anims = m_model->GetAnimations();
		for(unsigned int i=0; i<anims.size(); i++) {
			animSelector->AddOption(anims[i]->GetName());
		}
		if (anims.size())
			animSelector->SetSelectedOption(anims[0]->GetName());
	}
	animSelector->Layout();
	OnAnimChanged(0, animSelector->GetSelectedOption());
}

void ModelViewer::UpdateCamera()
{
	static const float BASE_ZOOM_RATE = 1.0f / 12.0f;
	float zoomRate = (BASE_ZOOM_RATE * 8.0f) * m_frameTime;
	float rotateRate = 25.f * m_frameTime;
	float moveRate = 10.0f * m_frameTime;

	if (m_keyStates[SDLK_LSHIFT]) {
		zoomRate *= 8.0f;
		moveRate *= 4.0f;
		rotateRate *= 4.0f;
	}
	else if (m_keyStates[SDLK_RSHIFT]) {
		zoomRate *= 3.0f;
		moveRate *= 2.0f;
		rotateRate *= 2.0f;
	}

	if (m_options.mouselookEnabled) {
		const float degrees_per_pixel = 0.2f;
		if (!m_mouseButton[SDL_BUTTON_RIGHT]) {
			// yaw and pitch
			const float rot_y = degrees_per_pixel*m_mouseMotion[0];
			const float rot_x = degrees_per_pixel*m_mouseMotion[1];
			const matrix3x3f rot =
				matrix3x3f::RotateX(DEG2RAD(rot_x)) *
				matrix3x3f::RotateY(DEG2RAD(rot_y));

			m_viewRot = m_viewRot * rot;
		} else {
			// roll
			m_viewRot = m_viewRot * matrix3x3f::RotateZ(DEG2RAD(degrees_per_pixel * m_mouseMotion[0]));
		}

		vector3f motion(0.0f);
		if (m_keyStates[SDLK_w]) motion.z -= moveRate;
		if (m_keyStates[SDLK_s]) motion.z += moveRate;
		if (m_keyStates[SDLK_a]) motion.x -= moveRate;
		if (m_keyStates[SDLK_d]) motion.x += moveRate;
		if (m_keyStates[SDLK_q]) motion.y -= moveRate;
		if (m_keyStates[SDLK_e]) motion.y += moveRate;

		m_viewPos += m_viewRot * motion;
	} else {
		//zoom
		if (m_keyStates[SDLK_EQUALS] || m_keyStates[SDLK_KP_PLUS]) m_zoom -= zoomRate;
		if (m_keyStates[SDLK_MINUS] || m_keyStates[SDLK_KP_MINUS]) m_zoom += zoomRate;

		//zoom with mouse wheel
		if (m_mouseWheelUp) m_zoom -= BASE_ZOOM_RATE;
		if (m_mouseWheelDown) m_zoom += BASE_ZOOM_RATE;

		m_zoom = Clamp(m_zoom, -10.0f, 10.0f); // distance range: [baseDistance * 1/1024, baseDistance * 1024]

		//rotate
		if (m_keyStates[SDLK_UP]) m_rotX += rotateRate;
		if (m_keyStates[SDLK_DOWN]) m_rotX -= rotateRate;
		if (m_keyStates[SDLK_LEFT]) m_rotY += rotateRate;
		if (m_keyStates[SDLK_RIGHT]) m_rotY -= rotateRate;

		//mouse rotate when right button held
		if (m_mouseButton[SDL_BUTTON_RIGHT]) {
			m_rotY += 0.2f*m_mouseMotion[0];
			m_rotX += 0.2f*m_mouseMotion[1];
		}
	}
}

void ModelViewer::UpdateLights()
{
	using Graphics::Light;
	std::vector<Light> lights;

	switch(m_options.lightPreset) {
	case 0:
		//Front white
		lights.push_back(Light(Light::LIGHT_DIRECTIONAL, az_el_to_dir(90,0), Color::WHITE, Color::WHITE));
		lights.push_back(Light(Light::LIGHT_DIRECTIONAL, az_el_to_dir(0,-90), Color(13, 13, 26), Color::WHITE));
		break;
	case 1:
		//Two-point
		lights.push_back(Light(Light::LIGHT_DIRECTIONAL, az_el_to_dir(120,0), Color(230, 204, 204), Color::WHITE));
		lights.push_back(Light(Light::LIGHT_DIRECTIONAL, az_el_to_dir(-30,-90), Color(178, 128, 0), Color::WHITE));
		break;
	case 2:
		//Backlight
		lights.push_back(Light(Light::LIGHT_DIRECTIONAL, az_el_to_dir(-75,20), Color::WHITE, Color::WHITE));
		lights.push_back(Light(Light::LIGHT_DIRECTIONAL, az_el_to_dir(0,-90), Color(13, 13, 26), Color::WHITE));
		break;
	case 3:
		//4 lights
		lights.push_back(Light(Light::LIGHT_DIRECTIONAL, az_el_to_dir(0, 90), Color::YELLOW, Color::WHITE));
		lights.push_back(Light(Light::LIGHT_DIRECTIONAL, az_el_to_dir(0, -90), Color::GREEN, Color::WHITE));
		lights.push_back(Light(Light::LIGHT_DIRECTIONAL, az_el_to_dir(0, 45), Color::BLUE, Color::WHITE));
		lights.push_back(Light(Light::LIGHT_DIRECTIONAL, az_el_to_dir(0, -45), Color::WHITE, Color::WHITE));
		break;
	};

	m_renderer->SetLights(int(lights.size()), &lights[0]);
}

void ModelViewer::UpdatePatternList()
{
	patternSelector->Clear();

	if (m_model) {
		const SceneGraph::PatternContainer &pats = m_model->GetPatterns();
		for(unsigned int i=0; i<pats.size(); i++) {
			patternSelector->AddOption(pats[i].name);
		}
		if (pats.size() > 0)
			patternSelector->SetSelectedOption(pats[0].name);
	}

	m_ui->Layout();
}
