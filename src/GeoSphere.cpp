// Copyright © 2008-2025 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "GeoSphere.h"

#include "GameConfig.h"
#include "GeoPatch.h"
#include "GeoPatchContext.h"
#include "GeoPatchJobs.h"
#include "Pi.h"
#include "RefCounted.h"
#include "galaxy/AtmosphereParameters.h"
#include "galaxy/StarSystem.h"
#include "graphics/Drawables.h"
#include "graphics/Frustum.h"
#include "graphics/Graphics.h"
#include "graphics/Material.h"
#include "graphics/RenderState.h"
#include "graphics/Renderer.h"
#include "graphics/Texture.h"
#include "graphics/TextureBuilder.h"
#include "graphics/Types.h"
#include "graphics/VertexArray.h"
#include "perlin.h"
#include "utils.h"
#include "vcacheopt/vcacheopt.h"
#include <algorithm>
#include <deque>

RefCountedPtr<GeoPatchContext> GeoSphere::s_patchContext;

// must be odd numbers
static const int detail_edgeLen[5] = {
	//7, 15, 25, 35, 55 -- old non power-of-2+1 values
	// some detail settings duplicated intentionally
	// in real terms provides only 3 settings
	// however this value is still used for gas giants
	// with 5 distinct settings elsewhere
	9, 17, 17, 33, 33
};

static const double gs_targetPatchTriLength(100.0);
static std::vector<GeoSphere *> s_allGeospheres;
static Uint32 s_debugFlags = GeoSphere::DebugFlags::DEBUG_NONE;

void GeoSphere::InitGeoSphere()
{
	s_patchContext.Reset(new GeoPatchContext(detail_edgeLen[Pi::detail.planets > 4 ? 4 : Pi::detail.planets]));
}

void GeoSphere::UninitGeoSphere()
{
	assert(s_patchContext.Unique());
	s_patchContext.Reset();
}

static void print_info(const SystemBody *sbody, const Terrain *terrain)
{
	Log::Verbose(
		"Geosphere Init for {}:\n" \
		"\theight fractal: {}\n" \
		"\tcolour fractal: {}\n" \
		"\tseed: {}\n",
		sbody->GetName(), terrain->GetHeightFractalName(), terrain->GetColorFractalName(), sbody->GetSeed());
}

// static
void GeoSphere::UpdateAllGeoSpheres()
{
	PROFILE_SCOPED()
	for (std::vector<GeoSphere *>::iterator i = s_allGeospheres.begin(); i != s_allGeospheres.end(); ++i) {
		(*i)->Update();
	}
}

// static
void GeoSphere::OnChangeGeoSphereDetailLevel()
{
	s_patchContext.Reset(new GeoPatchContext(detail_edgeLen[Pi::detail.planets > 4 ? 4 : Pi::detail.planets]));

	// reinit the geosphere terrain data
	for (std::vector<GeoSphere *>::iterator i = s_allGeospheres.begin(); i != s_allGeospheres.end(); ++i) {
		// clearout anything we don't need
		(*i)->Reset();

		// reinit the terrain with the new settings
		(*i)->m_terrain.Reset(Terrain::InstanceTerrain((*i)->GetSystemBody()));
		print_info((*i)->GetSystemBody(), (*i)->m_terrain.Get());

		// Reload the atmosphere material (scattering option)
		if ((*i)->m_atmosphereMaterial.Valid()) {
			(*i)->CreateAtmosphereMaterial();
		}
	}
}

//static
bool GeoSphere::OnAddQuadSplitResult(const SystemPath &path, SQuadSplitResult *res)
{
	// Find the correct GeoSphere via it's system path, and give it the split result
	for (std::vector<GeoSphere *>::iterator i = s_allGeospheres.begin(), iEnd = s_allGeospheres.end(); i != iEnd; ++i) {
		if (path == (*i)->GetSystemBody()->GetPath()) {
			(*i)->AddQuadSplitResult(res);
			return true;
		}
	}
	// GeoSphere not found to return the data to, cancel and delete it instead
	if (res) {
		res->OnCancel();
		delete res;
	}
	return false;
}

//static
bool GeoSphere::OnAddSingleSplitResult(const SystemPath &path, SSingleSplitResult *res)
{
	// Find the correct GeoSphere via it's system path, and give it the split result
	for (std::vector<GeoSphere *>::iterator i = s_allGeospheres.begin(), iEnd = s_allGeospheres.end(); i != iEnd; ++i) {
		if (path == (*i)->GetSystemBody()->GetPath()) {
			(*i)->AddSingleSplitResult(res);
			return true;
		}
	}
	// GeoSphere not found to return the data to, cancel and delete it instead
	if (res) {
		res->OnCancel();
		delete res;
	}
	return false;
}

//static
void GeoSphere::SetDebugFlags(Uint32 flags)
{
	s_debugFlags = flags;
}

//static
Uint32 GeoSphere::GetDebugFlags()
{
	return s_debugFlags;
}

void GeoSphere::Reset()
{
	{
		std::deque<SSingleSplitResult *>::iterator iter = mSingleSplitResults.begin();
		while (iter != mSingleSplitResults.end()) {
			// finally pass SplitResults
			SSingleSplitResult *psr = (*iter);
			assert(psr);

			psr->OnCancel();

			// tidyup
			delete psr;

			// Next!
			++iter;
		}
		mSingleSplitResults.clear();
	}

	{
		std::deque<SQuadSplitResult *>::iterator iter = mQuadSplitResults.begin();
		while (iter != mQuadSplitResults.end()) {
			// finally pass SplitResults
			SQuadSplitResult *psr = (*iter);
			assert(psr);

			psr->OnCancel();

			// tidyup
			delete psr;

			// Next!
			++iter;
		}
		mQuadSplitResults.clear();
	}

	for (int p = 0; p < NUM_PATCHES; p++) {
		// delete patches
		if (m_patches[p]) {
			m_patches[p].reset();
		}
	}

	CalculateMaxPatchDepth();

	m_visiblePatches.reserve(1024);

	m_initStage = eBuildFirstPatches;
}

GeoSphere::GeoSphere(const SystemBody *body) :
	BaseSphere(body),
	m_hasTempCampos(false),
	m_tempCampos(0.0),
	m_tempFrustum(800, 600, 0.5, 1.0, 1000.0),
	m_initStage(eBuildFirstPatches),
	m_maxDepth(0)
{
	print_info(body, m_terrain.Get());

	s_allGeospheres.emplace_back(this);

	CalculateMaxPatchDepth();

	m_visiblePatches.reserve(1024);

	if (Pi::config->Int("SortGeoPatches") == 0) {
		SetDebugFlags(GetDebugFlags() & ~DebugFlags::DEBUG_SORTGEOPATCHES);
	} else {
		SetDebugFlags(GetDebugFlags() | DebugFlags::DEBUG_SORTGEOPATCHES);
	}

	//SetUpMaterials is not called until first Render since light count is zero :)

	InitCityRegions(body);
}

GeoSphere::~GeoSphere()
{
	// update thread should not be able to access us now, so we can safely continue to delete
	assert(std::count(s_allGeospheres.begin(), s_allGeospheres.end(), this) == 1);
	s_allGeospheres.erase(std::find(s_allGeospheres.begin(), s_allGeospheres.end(), this));
}

bool GeoSphere::AddQuadSplitResult(SQuadSplitResult *res)
{
	bool result = false;
	assert(res);
	assert(mQuadSplitResults.size() < MAX_SPLIT_OPERATIONS);
	if (mQuadSplitResults.size() < MAX_SPLIT_OPERATIONS) {
		mQuadSplitResults.emplace_back(res);
		result = true;
	}
	return result;
}

bool GeoSphere::AddSingleSplitResult(SSingleSplitResult *res)
{
	bool result = false;
	assert(res);
	assert(mSingleSplitResults.size() < MAX_SPLIT_OPERATIONS);
	if (mSingleSplitResults.size() < MAX_SPLIT_OPERATIONS) {
		mSingleSplitResults.emplace_back(res);
		result = true;
	}
	return result;
}

void GeoSphere::ProcessSplitResults()
{
	PROFILE_SCOPED()
	// now handle the single split results that define the base level of the quad tree
	{
		std::deque<SSingleSplitResult *>::iterator iter = mSingleSplitResults.begin();
		while (iter != mSingleSplitResults.end()) {
			// finally pass SplitResults
			SSingleSplitResult *psr = (*iter);
			assert(psr);

			const int32_t faceIdx = psr->face();
			if (m_patches[faceIdx]) {
				m_patches[faceIdx]->ReceiveHeightmap(psr);
			} else {
				psr->OnCancel();
			}

			// tidyup
			delete psr;

			// Next!
			++iter;
		}
		mSingleSplitResults.clear();
	}

	// now handle the quad split results
	{
		std::deque<SQuadSplitResult *>::iterator iter = mQuadSplitResults.begin();
		while (iter != mQuadSplitResults.end()) {
			// finally pass SplitResults
			SQuadSplitResult *psr = (*iter);
			assert(psr);

			const int32_t faceIdx = psr->face();
			if (m_patches[faceIdx]) {
				m_patches[faceIdx]->ReceiveHeightmaps(psr);
			} else {
				psr->OnCancel();
			}

			// tidyup
			delete psr;

			// Next!
			++iter;
		}
		mQuadSplitResults.clear();
	}
}

void GeoSphere::BuildFirstPatches()
{
	PROFILE_SCOPED()
	assert(!m_patches[0]);
	if (m_patches[0])
		return;

	CalculateMaxPatchDepth();

	// generate root face patches of the cube/sphere
	static const vector3d p1 = (vector3d(1, 1, 1)).Normalized();
	static const vector3d p2 = (vector3d(-1, 1, 1)).Normalized();
	static const vector3d p3 = (vector3d(-1, -1, 1)).Normalized();
	static const vector3d p4 = (vector3d(1, -1, 1)).Normalized();
	static const vector3d p5 = (vector3d(1, 1, -1)).Normalized();
	static const vector3d p6 = (vector3d(-1, 1, -1)).Normalized();
	static const vector3d p7 = (vector3d(-1, -1, -1)).Normalized();
	static const vector3d p8 = (vector3d(1, -1, -1)).Normalized();

	const uint64_t maxShiftDepth = GeoPatchID::MAX_SHIFT_DEPTH;

	m_patches[0].reset(new GeoPatch(s_patchContext, this, p1, p2, p3, p4, 0, (0ULL << maxShiftDepth)));
	m_patches[1].reset(new GeoPatch(s_patchContext, this, p4, p3, p7, p8, 0, (1ULL << maxShiftDepth)));
	m_patches[2].reset(new GeoPatch(s_patchContext, this, p1, p4, p8, p5, 0, (2ULL << maxShiftDepth)));
	m_patches[3].reset(new GeoPatch(s_patchContext, this, p2, p1, p5, p6, 0, (3ULL << maxShiftDepth)));
	m_patches[4].reset(new GeoPatch(s_patchContext, this, p3, p2, p6, p7, 0, (4ULL << maxShiftDepth)));
	m_patches[5].reset(new GeoPatch(s_patchContext, this, p8, p7, p6, p5, 0, (5ULL << maxShiftDepth)));

	for (int i = 0; i < NUM_PATCHES; i++) {
		m_patches[i]->RequestSinglePatch();
	}

	m_initStage = eRequestedFirstPatches;
}

void GeoSphere::CalculateMaxPatchDepth()
{
	const double circumference = 2.0 * M_PI * m_sbody->GetRadius();
	// calculate length of each edge segment (quad) times 4 due to that being the number around the sphere (1 per side, 4 sides for Root).
	double edgeMetres = circumference / double(s_patchContext->GetEdgeLen() * 8);
	// find out what depth we reach the desired resolution
	while (edgeMetres > gs_targetPatchTriLength && m_maxDepth < GEOPATCH_MAX_DEPTH) {
		edgeMetres *= 0.5;
		++m_maxDepth;
	}
}

void GeoSphere::Update()
{
	switch (m_initStage) {
	case eBuildFirstPatches:
		BuildFirstPatches();
		break;
	case eRequestedFirstPatches: {
		ProcessSplitResults();
		uint8_t numValidPatches = 0;
		for (int i = 0; i < NUM_PATCHES; i++) {
			if (m_patches[i]->HasHeightData()) {
				++numValidPatches;
			}
		}
		m_initStage = (NUM_PATCHES == numValidPatches) ? eReceivedFirstPatches : eRequestedFirstPatches;
	} break;
	case eReceivedFirstPatches: {
		for (int i = 0; i < NUM_PATCHES; i++) {
			m_patches[i]->SetNeedToUpdateVBOs();
		}
		m_initStage = eDefaultUpdateState;
	} break;
	case eDefaultUpdateState:
		if (m_hasTempCampos) {
			ProcessSplitResults();
			for (int i = 0; i < NUM_PATCHES; i++) {
				m_patches[i]->LODUpdate(m_tempCampos, m_tempFrustum);
			}
			ProcessQuadSplitRequests();
		}
		break;
	}
}

void GeoSphere::AddQuadSplitRequest(double dist, SQuadSplitRequest *pReq, GeoPatch *pPatch)
{
	mQuadSplitRequests.emplace_back(dist, pReq, pPatch);
}

void GeoSphere::ProcessQuadSplitRequests()
{
	std::sort(mQuadSplitRequests.begin(), mQuadSplitRequests.end(), [](TDistanceRequest &a, TDistanceRequest &b) { return a.mDistance < b.mDistance; });

	for (auto iter : mQuadSplitRequests) {
		SQuadSplitRequest *ssrd = iter.mpRequest;
		iter.mpRequester->ReceiveJobHandle(Pi::GetAsyncJobQueue()->Queue(new QuadPatchJob(ssrd)));
	}
	mQuadSplitRequests.clear();
}

void GeoSphere::Render(Graphics::Renderer *renderer, const matrix4x4d &modelView, vector3d campos, const float radius, const std::vector<Camera::Shadow> &shadows)
{
	PROFILE_SCOPED()
	// store this for later usage in the update method.
	m_tempCampos = campos;
	m_hasTempCampos = true;

	if (m_initStage < eDefaultUpdateState)
		return;

	matrix4x4d trans = modelView;
	trans.Translate(-campos.x, -campos.y, -campos.z);
	renderer->SetTransform(matrix4x4f(trans)); //need to set this for the following line to work
	matrix4x4d modv = matrix4x4d(renderer->GetTransform());
	matrix4x4d proj = matrix4x4d(renderer->GetProjection());
	Graphics::Frustum frustum(modv, proj);
	m_tempFrustum = frustum;

	// no frustum test of entire geosphere, since Space::Render does this
	// for each body using its GetBoundingRadius() value

	//First draw - create materials (they do not change afterwards)
	if (!m_surfaceMaterial)
		SetUpMaterials();

	//Update material parameters
	SetMaterialParameters(trans, radius, shadows, m_atmosphereParameters);

	if (m_atmosphereMaterial.Valid() && m_atmosphereParameters.atmosDensity > 0.0) {
		// make atmosphere sphere slightly bigger than required so
		// that the edges of the pixel shader atmosphere jizz doesn't
		// show ugly polygonal angles
		DrawAtmosphereSurface(renderer, trans, campos,
			m_atmosphereParameters.atmosRadius * 1.02,
			m_atmosphereMaterial);
	}

	Color ambient;
	Color &emission = m_surfaceMaterial->emissive;

	// save old global ambient
	const Color oldAmbient = renderer->GetAmbientColor();

	if ((GetSystemBody()->GetSuperType() == SystemBody::SUPERTYPE_STAR) || (GetSystemBody()->GetType() == SystemBody::TYPE_BROWN_DWARF)) {
		// stars should emit light and terrain should be visible from distance
		ambient.r = ambient.g = ambient.b = 51;
		ambient.a = 255;
		emission = StarSystem::starRealColors[GetSystemBody()->GetType()];
		emission.a = 255;
	} else {
		// give planet some ambient lighting if the viewer is close to it
		double camdist = 0.1 / campos.LengthSqr();
		// why the fuck is this returning 0.1 when we are sat on the planet??
		// JJ: Because campos is relative to a unit-radius planet - 1.0 at the surface
		// XXX oh well, it is the value we want anyway...
		ambient.r = ambient.g = ambient.b = camdist * 255;
		ambient.a = 255;
	}

	renderer->SetAmbientColor(ambient);

	renderer->SetTransform(matrix4x4f(modelView));

	if (s_debugFlags & GeoSphere::DebugFlags::DEBUG_WIREFRAME)
		renderer->SetWireFrameMode(true);

	if (s_debugFlags & GeoSphere::DebugFlags::DEBUG_SORTGEOPATCHES) {
		// Gather the patches that could be rendered
		for (int i = 0; i < NUM_PATCHES; i++) {
			m_patches[i]->GatherRenderablePatches(m_visiblePatches, renderer, campos, frustum);
		}

		// distance sort the patches
		std::sort(m_visiblePatches.begin(), m_visiblePatches.end(), [&, campos](const std::pair<double, GeoPatch *> &a, const std::pair<double, GeoPatch *> &b) {
			return (a.first) < (b.first);
		});

		// cull occluded patches somehow?
		// create frustum from corner points, something vertical, and the campos??? Cull anything within that frustum?

		// render the sorted patches
		for (std::pair<double, GeoPatch *> &pPatch : m_visiblePatches) {
			pPatch.second->RenderImmediate(renderer, campos, modelView);
		}

		// must clear this after each render otherwise it just accumulates every patch ever drawn!
		m_visiblePatches.clear();
	} else {
		for (int i = 0; i < NUM_PATCHES; i++) {
			m_patches[i]->Render(renderer, campos, modelView, frustum);
		}
	}

	if (s_debugFlags & GeoSphere::DebugFlags::DEBUG_WIREFRAME)
		renderer->SetWireFrameMode(false);

	renderer->SetAmbientColor(oldAmbient);

	renderer->GetStats().AddToStatCount(Graphics::Stats::STAT_PLANETS, 1);
}

void GeoSphere::SetUpMaterials()
{
	// XXX: this has to be synced with the vertex format used in GeoPatch
	auto vtxFormat = Graphics::VertexFormatDesc::FromAttribSet(Graphics::ATTRIB_POSITION | Graphics::ATTRIB_NORMAL | Graphics::ATTRIB_DIFFUSE | Graphics::ATTRIB_UV0);

	m_atmosphereParameters = GetSystemBody()->CalcAtmosphereParams();
	// normal star has a different setup path than geosphere terrain does
	if (GetSystemBody()->GetSuperType() == SystemBody::SUPERTYPE_STAR) {
		Graphics::MaterialDescriptor surfDesc;
		surfDesc.lighting = false;
		Graphics::RenderStateDesc rsd;
		m_surfaceMaterial.Reset(Pi::renderer->CreateMaterial("geosphere_star", surfDesc, rsd, vtxFormat));
	} else {
		// Request material for this star or planet, with or without
		// atmosphere. Separate material for surface and sky.
		Graphics::MaterialDescriptor surfDesc;
		const Uint32 effect_flags = m_terrain->GetSurfaceEffects();
		if (effect_flags & Terrain::EFFECT_LAVA)
			surfDesc.effect = Graphics::EFFECT_GEOSPHERE_TERRAIN_WITH_LAVA;
		else if (effect_flags & Terrain::EFFECT_WATER)
			surfDesc.effect = Graphics::EFFECT_GEOSPHERE_TERRAIN_WITH_WATER;
		else
			surfDesc.effect = Graphics::EFFECT_GEOSPHERE_TERRAIN;

		if ((GetSystemBody()->GetType() == SystemBody::TYPE_BROWN_DWARF) ||
			(GetSystemBody()->GetType() == SystemBody::TYPE_STAR_M)) {
			//dim star (emits and receives light)
			surfDesc.lighting = true;
			// emits light, so atmosphere should not affect visible light
			surfDesc.quality &= ~Graphics::HAS_ATMOSPHERE;
		} else {
			//planetoid with or without atmosphere
			surfDesc.lighting = true;
			if (m_atmosphereParameters.atmosDensity > 0.0) {
				surfDesc.quality |= Graphics::HAS_ATMOSPHERE;
			}
		}

		//solid blendmode
		Graphics::RenderStateDesc rsd;
		surfDesc.quality |= Graphics::HAS_ECLIPSES;
		m_surfaceMaterial.Reset(Pi::renderer->CreateMaterial("geosphere_terrain", surfDesc, rsd, vtxFormat));

		m_texHi.Reset(Graphics::TextureBuilder::Model("textures/high.dds").GetOrCreateTexture(Pi::renderer, "model"));
		m_texLo.Reset(Graphics::TextureBuilder::Model("textures/low.dds").GetOrCreateTexture(Pi::renderer, "model"));
		m_surfaceMaterial->SetTexture("texture0"_hash, m_texHi.Get());
		m_surfaceMaterial->SetTexture("texture1"_hash, m_texLo.Get());
	}

	CreateAtmosphereMaterial();
}

void GeoSphere::CreateAtmosphereMaterial()
{
	Graphics::RenderStateDesc rsd;
	// atmosphere is blended over the background
	rsd.blendMode = Graphics::BLEND_ALPHA_ONE;
	rsd.cullMode = Graphics::CULL_FRONT;
	rsd.depthWrite = false;

	Graphics::MaterialDescriptor skyDesc;
	skyDesc.effect = Graphics::EFFECT_GEOSPHERE_SKY;
	skyDesc.lighting = true;
	skyDesc.quality |= Graphics::HAS_ECLIPSES;

	Graphics::VertexFormatDesc atmosVtxFmt = m_atmos->GetVertexFormat();

	const int scattering = Pi::config->Int("RealisticScattering");
	switch (scattering) {
	case 1:
		m_atmosphereMaterial.Reset(Pi::renderer->CreateMaterial("rayleigh_fast", skyDesc, rsd, atmosVtxFmt));
		break;
	case 2:
		m_atmosphereMaterial.Reset(Pi::renderer->CreateMaterial("rayleigh_accurate", skyDesc, rsd, atmosVtxFmt));
		break;
	default:
		m_atmosphereMaterial.Reset(Pi::renderer->CreateMaterial("geosphere_sky", skyDesc, rsd, atmosVtxFmt));
		break;
	}
}

static constexpr double TARGET_CITY_RADIUS = 600.0;

// Set up region data for each of the system body's child surface starports
void GeoSphere::InitCityRegions(const SystemBody *sb)
{
	m_regions.clear();

	if (!sb->HasChildren()) {
		return;
	}

	const double planetRadius = sb->GetRadius();
	m_dynamicRangeHeight = 60.0 / planetRadius; //in radii

	const double delta = (0.75 * TARGET_CITY_RADIUS / planetRadius);
	constexpr double citySize_m = TARGET_CITY_RADIUS * 1.1;

	constexpr size_t DEFAULT_RESERVE = 8;
	m_regions.reserve(DEFAULT_RESERVE);

	// step through the planet's sbody's children and set up regions for surface starports
	for (std::vector<SystemBody *>::const_iterator i = sb->GetChildren().begin(); i != sb->GetChildren().end(); i++) {
		if ((*i)->GetType() == SystemBody::TYPE_STARPORT_SURFACE) {
			// calculate position of starport
			const vector3d pos = ((*i)->GetOrbit().GetPlane() * vector3d(0, 1, 0));

			// set up regions which contain the details for region implementation
			Region rt;
			rt.position = pos;

			// height in planet radii
			rt.height = GetHeight(pos);

			// Calculate average variation of four points about star port
			// points do not need to be on the planet surface
			double avgVariation = fabs(GetHeight(vector3d(pos.x + delta, pos.y, pos.z)) - rt.height);
			avgVariation += fabs(GetHeight(vector3d(pos.x - delta, pos.y, pos.z)) - rt.height);
			avgVariation += fabs(GetHeight(vector3d(pos.x, pos.y, pos.z + delta)) - rt.height);
			avgVariation += fabs(GetHeight(vector3d(pos.x, pos.y, pos.z - delta)) - rt.height);
			avgVariation *= (1 / 4.0);
			rt.heightVariation = (1.0 / planetRadius) + 0.625 * avgVariation;

			// angle between city center/boundary = 2pi*city size/(perimeter great circle = 2pi r)
			// city center pos and current point will be dotted, and compared against size
			const double size = fabs(cos(std::min(citySize_m, 0.2 * sb->GetRadius()) / (sb->GetRadius())));
			rt.outer = size;
			rt.inner = (1.0 - size) * 0.5 + size;

			m_regions.emplace_back(rt);
		}
	}

	// reduce wasted space, maybe
	m_regions.resize(m_regions.size());
}

void GeoSphere::ApplySimpleHeightRegions(double &h, const vector3d &p) const
{
	for (size_t i = 0; i < m_regions.size(); i++) {
		const Region &rt = m_regions[i];
		const double posDotp = rt.position.Dot(p);
		if (posDotp > rt.outer) {
			// target height
			const double th = rt.height;

			// maximum variation in height with respect to target height
			const double delta_h = fabs(h - th);
			const double neg = (h - th > 0.0) ? 1.0 : -1.0;

			// Make up an expression to compress delta_h:
			// Compress delta_h between 0 and 1
			//    1.1 use compression of the form c = (delta_h+a)/(a+(delta_h+a)) (eqn. 1)
			//    1.2 this gives c in the interval [0.5, 1] for delta_h [0, +inf] with c=0.5 at delta_h=0.
			//  2.0 Use compressed_h = dynamic range*(sign(h-th)*(c-0.5)) (eqn. 2) to get h between th-0.5*dynamic range, th+0.5*dynamic range

			// Choosing a value for a
			//    3.1 c [0.5, 0.8] occurs when delta_h [a to 3a] (3x difference) which is roughly the expressible range (above or below that the function changes slowly)
			//    3.2 Find an expression for the expected variation and divide by around 3

			// It may become necessary calculate expected variation based on intermediate quantities generated (e.g. distribution fractals)
			// or to store a per planet estimation of variation when fracdefs are calculated.
			const double variationEstimate = rt.heightVariation;
			const double a = variationEstimate * (1.0 / 3.0); // point 3.2

			const double c = (delta_h + a) / (2.0 * a + delta_h);					 // point 1.1
			const double compressed_h = m_dynamicRangeHeight * (neg * (c - 0.5)) + th; // point 2.0

			// blends from compressed height-terrain height as pos goes inner to outer
			h = MathUtil::Lerp(h, compressed_h, Clamp((posDotp - rt.outer) / (rt.inner - rt.outer), 0.0, 1.0));
			break;
		}
	}
}

const GeoSphere::Region* GeoSphere::FindNearestRegion(const vector3d &p, double &posDotPOut) const
{
	bool bFound = false;
	size_t validIndex = 0;
	double lastPosDotP = -1.0; // Dot-product values closer to 1 are closer to alignment
	for (size_t i = 0; i < m_regions.size(); i++) {
		const Region &rt = m_regions[i];
		const double posDotp = rt.position.Dot(p);

		// if-less
		bFound = (posDotp > rt.outer) ? true : bFound; 
		const bool bestDot = (posDotp > lastPosDotP);
		validIndex = (bFound && bestDot) ? i : validIndex;
		lastPosDotP = (bFound && bestDot) ? posDotp : lastPosDotP;
	}
	posDotPOut = lastPosDotP;
	return (bFound) ? &m_regions[validIndex] : nullptr;
}

std::vector<const GeoSphere::Region *> GeoSphere::FindPatchRegions(const vector3d &p, const vector3d &v0, const vector3d &v1, const vector3d &v2, const vector3d &v3) const
{
	std::vector<const GeoSphere::Region *> foundRegions;
	for (size_t i = 0; i < m_regions.size(); i++) {
		// check each region against the four bounding points

		// build a frustum from the four points maybe?
		// test if spheres (pos + outer) are inside frustum?

		// or can i use Dot product somehow?

		// or construct a plane, project region->pos onto it and test if it's inside?
	}
	return foundRegions;
}

void GeoSphere::ApplyHeightRegion(double &h, const double posDotP, const Region *region) const
{
	if (!region)
		return;
	
	const Region &rt = *region;
	
	// target height
	const double th = rt.height;

	// maximum variation in height with respect to target height
	const double delta_h = fabs(h - th);
	const double neg = (h - th > 0.0) ? 1.0 : -1.0;

	// Make up an expression to compress delta_h:
	// Compress delta_h between 0 and 1
	//    1.1 use compression of the form c = (delta_h+a)/(a+(delta_h+a)) (eqn. 1)
	//    1.2 this gives c in the interval [0.5, 1] for delta_h [0, +inf] with c=0.5 at delta_h=0.
	//  2.0 Use compressed_h = dynamic range*(sign(h-th)*(c-0.5)) (eqn. 2) to get h between th-0.5*dynamic range, th+0.5*dynamic range

	// Choosing a value for a
	//    3.1 c [0.5, 0.8] occurs when delta_h [a to 3a] (3x difference) which is roughly the expressible range (above or below that the function changes slowly)
	//    3.2 Find an expression for the expected variation and divide by around 3

	// It may become necessary calculate expected variation based on intermediate quantities generated (e.g. distribution fractals)
	// or to store a per planet estimation of variation when fracdefs are calculated.
	const double variationEstimate = rt.heightVariation;
	const double a = variationEstimate * (1.0 / 3.0); // point 3.2

	const double c = (delta_h + a) / (2.0 * a + delta_h);					   // point 1.1
	const double compressed_h = m_dynamicRangeHeight * (neg * (c - 0.5)) + th; // point 2.0

	// blends from compressed height-terrain height as pos goes inner to outer
	h = MathUtil::Lerp(h, compressed_h, Clamp((posDotP - rt.outer) / (rt.inner - rt.outer), 0.0, 1.0));
}
