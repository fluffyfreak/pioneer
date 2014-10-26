// Copyright © 2008-2014 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "libs.h"
#include "GeoSphere.h"
#include "GeoPatchContext.h"
#include "GeoPatch.h"
#include "GeoPatchJobs.h"
#include "perlin.h"
#include "Pi.h"
#include "RefCounted.h"
#include "graphics/Material.h"
#include "graphics/Renderer.h"
#include "graphics/Frustum.h"
#include "graphics/Graphics.h"
#include "graphics/VertexArray.h"
#include "vcacheopt/vcacheopt.h"
#include <deque>
#include <algorithm>

//------------------------------------------------------------

static const int MAX_SUBDIVS = 5;
static const float ICOSX = 0.525731112119133f;
static const float ICOSZ = 0.850650808352039f;

static const vector3f icosahedron_vertices[12] = {
	vector3f(-ICOSX, 0.0, ICOSZ), vector3f(ICOSX, 0.0, ICOSZ), vector3f(-ICOSX, 0.0, -ICOSZ), vector3f(ICOSX, 0.0, -ICOSZ),
	vector3f(0.0, ICOSZ, ICOSX), vector3f(0.0, ICOSZ, -ICOSX), vector3f(0.0, -ICOSZ, ICOSX), vector3f(0.0, -ICOSZ, -ICOSX),
	vector3f(ICOSZ, ICOSX, 0.0), vector3f(-ICOSZ, ICOSX, 0.0), vector3f(ICOSZ, -ICOSX, 0.0), vector3f(-ICOSZ, -ICOSX, 0.0)
};

static const int icosahedron_faces[20][3] = {
	{0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},
	{8,10,1}, {8,3,10},{5,3,8}, {5,2,3}, {2,7,3},
	{7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6},
	{6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11}
};

// XXX merge these with their fracdef versions
inline double octavenoise(int octaves, const double persistence, const double lacunarity, const vector3d &p) {
	//assert(persistence <= (1.0 / lacunarity));
	double n = 0;
	double amplitude = persistence;
	double frequency = 1.0;
	while (octaves--) {
		n += amplitude * noise(frequency*p);
		amplitude *= persistence;
		frequency *= lacunarity;
	}
	return n;
}

//#define PERTURB_VERTICES
Icosahedron::Icosahedron(Graphics::Renderer *renderer, RefCountedPtr<Graphics::Material> mat, Graphics::RenderState *state, int subdivs, float scale)
{
	using namespace Graphics;

	PROFILE_SCOPED()
	m_material = mat;
	m_renderState = state;

	subdivs = Clamp(subdivs, 0, MAX_SUBDIVS);
	scale = fabs(scale);
	matrix4x4f trans = matrix4x4f::Identity();
	trans.Scale(scale, scale, scale);

	//reserve some data
	VertexArray vts(ATTRIB_POSITION | ATTRIB_NORMAL | ATTRIB_UV0, 256);
	std::vector<Uint32> indices;

#ifdef PERTURB_VERTICES
	const vector3f v0 = trans * icosahedron_vertices[0];
	const vector3f v1 = trans * icosahedron_vertices[1];
	const float len = (v1 - v0).Length();

	vector3f ico_verts[12];
#endif // PERTURB_VERTICES

	//initial vertices
	int vi[12];
	for (int i=0; i<12; i++) {
		const vector3f &v = icosahedron_vertices[i];

#ifdef PERTURB_VERTICES
		// basis vectors
		vector3f up(0.0f, 1.0f, 0.0f); // dummy initial value of up
		const vector3f right(v.Cross(up));
		up = v.Cross(right);

		const float s = (len * octavenoise(8, 0.5, 2.0, vector3d(v))) * 0.3;

		const vector3f fv = (v + (right * s) + (up * s)).Normalized();
		ico_verts[i] = fv;
		vi[i] = AddVertex(vts, trans * fv, fv);
#endif // PERTURB_VERTICES
		vi[i] = AddVertex(vts, trans * v, v);
	}

	//subdivide
	for (int i=0; i<20; i++) {
#ifdef PERTURB_VERTICES
		Subdivide(vts, indices, trans, ico_verts[icosahedron_faces[i][0]],
				ico_verts[icosahedron_faces[i][1]],
				ico_verts[icosahedron_faces[i][2]],
#endif // PERTURB_VERTICES
		Subdivide(vts, indices, trans, icosahedron_vertices[icosahedron_faces[i][0]],
				icosahedron_vertices[icosahedron_faces[i][1]],
				icosahedron_vertices[icosahedron_faces[i][2]],
				vi[icosahedron_faces[i][0]],
				vi[icosahedron_faces[i][1]],
				vi[icosahedron_faces[i][2]],
				subdivs);
	}

	//Create vtx & index buffers and copy data
	VertexBufferDesc vbd;
	vbd.attrib[0].semantic = ATTRIB_POSITION;
	vbd.attrib[0].format   = ATTRIB_FORMAT_FLOAT3;
	vbd.attrib[1].semantic = ATTRIB_NORMAL;
	vbd.attrib[1].format   = ATTRIB_FORMAT_FLOAT3;
	vbd.attrib[2].semantic = ATTRIB_UV0;
	vbd.attrib[2].format   = ATTRIB_FORMAT_FLOAT2;
	vbd.numVertices = vts.GetNumVerts();
	vbd.usage = BUFFER_USAGE_STATIC;
	m_vertexBuffer.reset(renderer->CreateVertexBuffer(vbd));
	m_vertexBuffer->Populate(vts);

	m_indexBuffer.reset(renderer->CreateIndexBuffer(indices.size(), BUFFER_USAGE_STATIC));
	Uint16 *idxPtr = m_indexBuffer->Map(Graphics::BUFFER_MAP_WRITE);
	for (auto it : indices) {
		*idxPtr = it;
		idxPtr++;
	}
	m_indexBuffer->Unmap();
}

void Icosahedron::Draw(Graphics::Renderer *r)
{
	PROFILE_SCOPED()
	r->DrawBufferIndexed(m_vertexBuffer.get(), m_indexBuffer.get(), m_renderState, m_material.Get());
}

int Icosahedron::AddVertex(Graphics::VertexArray &vts, const vector3f &v, const vector3f &n)
{
	PROFILE_SCOPED()
	vts.position.push_back(v);
	vts.normal.push_back(n);
	//http://www.mvps.org/directx/articles/spheremap.htm
	vts.uv0.push_back(vector2f(asinf(n.x)/M_PI+0.5f, asinf(n.y)/M_PI+0.5f));
	return vts.GetNumVerts() - 1;
}

void Icosahedron::AddTriangle(std::vector<Uint32> &indices, int i1, int i2, int i3)
{
	PROFILE_SCOPED()
	indices.push_back(i1);
	indices.push_back(i2);
	indices.push_back(i3);
}

void Icosahedron::Subdivide(Graphics::VertexArray &vts, std::vector<Uint32> &indices,
		const matrix4x4f &trans, const vector3f &v1, const vector3f &v2, const vector3f &v3,
		const int i1, const int i2, const int i3, const int depth)
{
	PROFILE_SCOPED()
	if (depth == 0) {
		AddTriangle(indices, i1, i3, i2);
		return;
	}

	const vector3f v12 = (v1+v2).Normalized();
	const vector3f v23 = (v2+v3).Normalized();
	const vector3f v31 = (v3+v1).Normalized();
	const int i12 = AddVertex(vts, trans * v12, v12);
	const int i23 = AddVertex(vts, trans * v23, v23);
	const int i31 = AddVertex(vts, trans * v31, v31);
	Subdivide(vts, indices, trans, v1, v12, v31, i1, i12, i31, depth-1);
	Subdivide(vts, indices, trans, v2, v23, v12, i2, i23, i12, depth-1);
	Subdivide(vts, indices, trans, v3, v31, v23, i3, i31, i23, depth-1);
	Subdivide(vts, indices, trans, v12, v23, v31, i12, i23, i31, depth-1);
}
//------------------------------------------------------------

RefCountedPtr<GeoPatchContext> GeoSphere::s_patchContext;

// must be odd numbers
static const int detail_edgeLen[5] = {
	7, 15, 25, 35, 55
};

static const double gs_targetPatchTriLength(100.0);

#define PRINT_VECTOR(_v) Output("%f,%f,%f\n", (_v).x, (_v).y, (_v).z);

static const int geo_sphere_edge_friends[NUM_PATCHES][4] = {
	{ 3, 4, 1, 2 },
	{ 0, 4, 5, 2 },
	{ 0, 1, 5, 3 },
	{ 0, 2, 5, 4 },
	{ 0, 3, 5, 1 },
	{ 1, 4, 3, 2 }
};

static std::vector<GeoSphere*> s_allGeospheres;

void GeoSphere::Init()
{
	s_patchContext.Reset(new GeoPatchContext(detail_edgeLen[Pi::detail.planets > 4 ? 4 : Pi::detail.planets]));
	assert(s_patchContext->GetEdgeLen() <= detail_edgeLen[4]);
}

void GeoSphere::Uninit()
{
	assert (s_patchContext.Unique());
	s_patchContext.Reset();
}

static void print_info(const SystemBody *sbody, const Terrain *terrain)
{
	Output(
		"%s:\n"
		"    height fractal: %s\n"
		"    colour fractal: %s\n"
		"    seed: %u\n",
		sbody->GetName().c_str(), terrain->GetHeightFractalName(), terrain->GetColorFractalName(), sbody->GetSeed());
}

// static
void GeoSphere::UpdateAllGeoSpheres()
{
	PROFILE_SCOPED()
	for(std::vector<GeoSphere*>::iterator i = s_allGeospheres.begin(); i != s_allGeospheres.end(); ++i)
	{
		(*i)->Update();
	}
}

// static
void GeoSphere::OnChangeDetailLevel()
{
	s_patchContext.Reset(new GeoPatchContext(detail_edgeLen[Pi::detail.planets > 4 ? 4 : Pi::detail.planets]));
	assert(s_patchContext->GetEdgeLen() <= detail_edgeLen[4]);

	// reinit the geosphere terrain data
	for(std::vector<GeoSphere*>::iterator i = s_allGeospheres.begin(); i != s_allGeospheres.end(); ++i)
	{
		// clearout anything we don't need
		(*i)->Reset();

		// reinit the terrain with the new settings
		(*i)->m_terrain.Reset(Terrain::InstanceTerrain((*i)->GetSystemBody()));
		print_info((*i)->GetSystemBody(), (*i)->m_terrain.Get());
	}
}

//static
bool GeoSphere::OnAddQuadSplitResult(const SystemPath &path, SQuadSplitResult *res)
{
	// Find the correct GeoSphere via it's system path, and give it the split result
	for(std::vector<GeoSphere*>::iterator i=s_allGeospheres.begin(), iEnd=s_allGeospheres.end(); i!=iEnd; ++i) {
		if( path == (*i)->GetSystemBody()->GetPath() ) {
			(*i)->AddQuadSplitResult(res);
			return true;
		}
	}
	// GeoSphere not found to return the data to, cancel and delete it instead
	if( res ) {
		res->OnCancel();
		delete res;
	}
	return false;
}

//static
bool GeoSphere::OnAddSingleSplitResult(const SystemPath &path, SSingleSplitResult *res)
{
	// Find the correct GeoSphere via it's system path, and give it the split result
	for(std::vector<GeoSphere*>::iterator i=s_allGeospheres.begin(), iEnd=s_allGeospheres.end(); i!=iEnd; ++i) {
		if( path == (*i)->GetSystemBody()->GetPath() ) {
			(*i)->AddSingleSplitResult(res);
			return true;
		}
	}
	// GeoSphere not found to return the data to, cancel and delete it instead
	if( res ) {
		res->OnCancel();
		delete res;
	}
	return false;
}

void GeoSphere::Reset()
{
	{
		std::deque<SSingleSplitResult*>::iterator iter = mSingleSplitResults.begin();
		while(iter!=mSingleSplitResults.end())
		{
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
		std::deque<SQuadSplitResult*>::iterator iter = mQuadSplitResults.begin();
		while(iter!=mQuadSplitResults.end())
		{
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

	for (int p=0; p<NUM_PATCHES; p++) {
		// delete patches
		if (m_patches[p]) {
			m_patches[p].reset();
		}
	}

	CalculateMaxPatchDepth();

	m_initStage = eBuildFirstPatches;
}

#define GEOSPHERE_TYPE	(GetSystemBody()->type)

GeoSphere::GeoSphere(const SystemBody *body) : BaseSphere(body),
	m_hasTempCampos(false), m_tempCampos(0.0), m_initStage(eBuildFirstPatches), m_maxDepth(0)
{
	print_info(body, m_terrain.Get());

	s_allGeospheres.push_back(this);

	CalculateMaxPatchDepth();

	//SetUpMaterials is not called until first Render since light count is zero :)
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
	assert(mQuadSplitResults.size()<MAX_SPLIT_OPERATIONS);
	if(mQuadSplitResults.size()<MAX_SPLIT_OPERATIONS) {
		mQuadSplitResults.push_back(res);
		result = true;
	}
	return result;
}

bool GeoSphere::AddSingleSplitResult(SSingleSplitResult *res)
{
	bool result = false;
	assert(res);
	assert(mSingleSplitResults.size()<MAX_SPLIT_OPERATIONS);
	if(mSingleSplitResults.size()<MAX_SPLIT_OPERATIONS) {
		mSingleSplitResults.push_back(res);
		result = true;
	}
	return result;
}

void GeoSphere::ProcessSplitResults()
{
	// now handle the single split results that define the base level of the quad tree
	{
		std::deque<SSingleSplitResult*>::iterator iter = mSingleSplitResults.begin();
		while(iter!=mSingleSplitResults.end())
		{
			// finally pass SplitResults
			SSingleSplitResult *psr = (*iter);
			assert(psr);

			const int32_t faceIdx = psr->face();
			if( m_patches[faceIdx] ) {
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
		std::deque<SQuadSplitResult*>::iterator iter = mQuadSplitResults.begin();
		while(iter!=mQuadSplitResults.end())
		{
			// finally pass SplitResults
			SQuadSplitResult *psr = (*iter);
			assert(psr);

			const int32_t faceIdx = psr->face();
			if( m_patches[faceIdx] ) {
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
	assert(!m_patches[0]);
	if(m_patches[0])
		return;

	// generate root face patches of the cube/sphere
	static const vector3d p1 = (vector3d( 1, 1, 1)).Normalized();
	static const vector3d p2 = (vector3d(-1, 1, 1)).Normalized();
	static const vector3d p3 = (vector3d(-1,-1, 1)).Normalized();
	static const vector3d p4 = (vector3d( 1,-1, 1)).Normalized();
	static const vector3d p5 = (vector3d( 1, 1,-1)).Normalized();
	static const vector3d p6 = (vector3d(-1, 1,-1)).Normalized();
	static const vector3d p7 = (vector3d(-1,-1,-1)).Normalized();
	static const vector3d p8 = (vector3d( 1,-1,-1)).Normalized();

	const uint64_t maxShiftDepth = GeoPatchID::MAX_SHIFT_DEPTH;

	m_patches[0].reset(new GeoPatch(s_patchContext, this, p1, p2, p3, p4, 0, (0ULL << maxShiftDepth)));
	m_patches[1].reset(new GeoPatch(s_patchContext, this, p4, p3, p7, p8, 0, (1ULL << maxShiftDepth)));
	m_patches[2].reset(new GeoPatch(s_patchContext, this, p1, p4, p8, p5, 0, (2ULL << maxShiftDepth)));
	m_patches[3].reset(new GeoPatch(s_patchContext, this, p2, p1, p5, p6, 0, (3ULL << maxShiftDepth)));
	m_patches[4].reset(new GeoPatch(s_patchContext, this, p3, p2, p6, p7, 0, (4ULL << maxShiftDepth)));
	m_patches[5].reset(new GeoPatch(s_patchContext, this, p8, p7, p6, p5, 0, (5ULL << maxShiftDepth)));
	for (int i=0; i<NUM_PATCHES; i++) {
		for (int j=0; j<4; j++) {
			m_patches[i]->SetEdgeFriend(j, m_patches[geo_sphere_edge_friends[i][j]].get());
		}
	}

	for (int i=0; i<NUM_PATCHES; i++) {
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
	while (edgeMetres>gs_targetPatchTriLength && m_maxDepth<GEOPATCH_MAX_DEPTH) {
		edgeMetres *= 0.5;
		++m_maxDepth;
	}
}

void GeoSphere::Update()
{
	switch(m_initStage)
	{
	case eBuildFirstPatches:
		BuildFirstPatches();
		break;
	case eRequestedFirstPatches:
		{
			ProcessSplitResults();
			uint8_t numValidPatches = 0;
			for (int i=0; i<NUM_PATCHES; i++) {
				if(m_patches[i]->HasHeightData()) {
					++numValidPatches;
				}
			}
			m_initStage = (NUM_PATCHES==numValidPatches) ? eReceivedFirstPatches : eRequestedFirstPatches;
		} break;
	case eReceivedFirstPatches:
		{
			for (int i=0; i<NUM_PATCHES; i++) {
				m_patches[i]->UpdateVBOs();
			}
			m_initStage = eDefaultUpdateState;
		} break;
	case eDefaultUpdateState:
		if(m_hasTempCampos) {
			ProcessSplitResults();
			for (int i=0; i<NUM_PATCHES; i++) {
				m_patches[i]->LODUpdate(m_tempCampos);
			}
		}
		break;
	}
}
#pragma optimize("",off)
void GeoSphere::Render(Graphics::Renderer *renderer, const matrix4x4d &modelView, vector3d campos, const float radius, const float scale, const std::vector<Camera::Shadow> &shadows)
{
	// store this for later usage in the update method.
	m_tempCampos = campos;
	m_hasTempCampos = true;

	if(m_initStage < eDefaultUpdateState)
		return;

	matrix4x4d trans = modelView;
	trans.Translate(-campos.x, -campos.y, -campos.z);
	renderer->SetTransform(trans); //need to set this for the following line to work
	matrix4x4d modv;
	matrix4x4d proj;
	matrix4x4ftod(renderer->GetCurrentModelView(), modv);
	matrix4x4ftod(renderer->GetCurrentProjection(), proj);
	Graphics::Frustum frustum( modv, proj );

	// no frustum test of entire geosphere, since Space::Render does this
	// for each body using its GetBoundingRadius() value

	//First draw - create materials (they do not change afterwards)
	if (!m_surfaceMaterial)
		SetUpMaterials();

	{
		//Update material parameters
		//XXX no need to calculate AP every frame
		m_materialParameters.atmosphere = GetSystemBody()->CalcAtmosphereParams();
		m_materialParameters.atmosphere.center = trans * vector3d(0.0, 0.0, 0.0);
		m_materialParameters.atmosphere.planetRadius = radius;
		m_materialParameters.atmosphere.scale = scale;

		m_materialParameters.shadows = shadows;

		m_surfaceMaterial->specialParameter0 = &m_materialParameters;

		if (m_materialParameters.atmosphere.atmosDensity > 0.0) {
			m_atmosphereMaterial->specialParameter0 = &m_materialParameters;

			// make atmosphere sphere slightly bigger than required so
			// that the edges of the pixel shader atmosphere jizz doesn't
			// show ugly polygonal angles
			DrawAtmosphereSurface(renderer, trans, campos,
				m_materialParameters.atmosphere.atmosRadius*1.01,
				m_atmosRenderState, m_atmosphereMaterial.get());
		}
	}

	// display the terrain height control-mesh
	{
		const float rad = m_materialParameters.atmosphere.atmosRadius * 1.01f;

		static std::unique_ptr<Icosahedron> m_ball;

		{
			renderer->SetTransform(trans * matrix4x4d::ScaleMatrix(rad, rad, rad));
			Pi::renderer->SetWireFrameMode(true);
			if( !m_ball.get() ) {
				RefCountedPtr<Graphics::Material> mat(Pi::renderer->CreateMaterial(Graphics::MaterialDescriptor()));
				mat->diffuse = Color4f(0.7f, 0.7f, 0.7f, 0.5f);
				m_ball.reset( new Icosahedron(Pi::renderer, mat, Pi::renderer->CreateRenderState(Graphics::RenderStateDesc()), 5, 1.0) );
			}

			if( m_ball.get() ) {
				m_ball->Draw( Pi::renderer );
			}
			Pi::renderer->SetWireFrameMode(false);
		}
	}

	Color ambient;
	Color &emission = m_surfaceMaterial->emissive;

	// save old global ambient
	const Color oldAmbient = renderer->GetAmbientColor();

	if ((GetSystemBody()->GetSuperType() == SystemBody::SUPERTYPE_STAR) || (GetSystemBody()->GetType() == SystemBody::TYPE_BROWN_DWARF)) {
		// stars should emit light and terrain should be visible from distance
		ambient.r = ambient.g = ambient.b = 51;
		ambient.a = 255;
		emission.r = StarSystem::starRealColors[GetSystemBody()->GetType()][0];
		emission.g = StarSystem::starRealColors[GetSystemBody()->GetType()][1];
		emission.b = StarSystem::starRealColors[GetSystemBody()->GetType()][2];
		emission.a = 255;
	}

	else {
		// give planet some ambient lighting if the viewer is close to it
		double camdist = campos.Length();
		camdist = 0.1 / (camdist*camdist);
		// why the fuck is this returning 0.1 when we are sat on the planet??
		// JJ: Because campos is relative to a unit-radius planet - 1.0 at the surface
		// XXX oh well, it is the value we want anyway...
		ambient.r = ambient.g = ambient.b = camdist * 255;
		ambient.a = 255;
	}

	renderer->SetAmbientColor(ambient);

	renderer->SetTransform(modelView);

	for (int i=0; i<NUM_PATCHES; i++) {
		m_patches[i]->Render(renderer, campos, modelView, frustum);
	}

	renderer->SetAmbientColor(oldAmbient);
}

void GeoSphere::SetUpMaterials()
{
	//solid
	Graphics::RenderStateDesc rsd;
	m_surfRenderState = Pi::renderer->CreateRenderState(rsd);

	//blended
	rsd.blendMode = Graphics::BLEND_ALPHA_ONE;
	rsd.depthWrite = false;
	m_atmosRenderState = Pi::renderer->CreateRenderState(rsd);

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
		surfDesc.quality &= ~Graphics::HAS_ATMOSPHERE;
	}
	else if (GetSystemBody()->GetSuperType() == SystemBody::SUPERTYPE_STAR) {
		//normal star
		surfDesc.lighting = false;
		surfDesc.quality &= ~Graphics::HAS_ATMOSPHERE;
	} else {
		//planetoid with or without atmosphere
		const SystemBody::AtmosphereParameters ap(GetSystemBody()->CalcAtmosphereParams());
		surfDesc.lighting = true;
		if(ap.atmosDensity > 0.0) {
			surfDesc.quality |= Graphics::HAS_ATMOSPHERE;
		} else {
			surfDesc.quality &= ~Graphics::HAS_ATMOSPHERE;
		}
	}

	const bool bEnableEclipse = (Pi::config->Int("DisableEclipse") == 0);
	if (bEnableEclipse) {
		surfDesc.quality |= Graphics::HAS_ECLIPSES;
	}
	m_surfaceMaterial.reset(Pi::renderer->CreateMaterial(surfDesc));

	{
		Graphics::MaterialDescriptor skyDesc;
		skyDesc.effect = Graphics::EFFECT_GEOSPHERE_SKY;
		skyDesc.lighting = true;
		if (bEnableEclipse) {
			skyDesc.quality |= Graphics::HAS_ECLIPSES;
		}
		m_atmosphereMaterial.reset(Pi::renderer->CreateMaterial(skyDesc));
	}
}
