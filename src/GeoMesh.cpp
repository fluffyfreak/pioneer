// Copyright © 2008-2012 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include <cassert>
#include "utils.h"
#include "GeoMesh.h"

#include "GeoContext.h"
#include "GeoPatch.h"

#include "perlin.h"
#include "Pi.h"
#include "RefCounted.h"
#include "graphics/Material.h"
#include "graphics/Renderer.h"
#include "graphics/Frustum.h"
#include "graphics/Graphics.h"
#include "graphics/VertexArray.h"
#include "graphics/gl2/GeoSphereMaterial.h"

static std::vector<GeoMesh*> s_allGeospheres;
// must be odd numbers
static const int detail_edgeLen[5] = {
	29, 61, 125, 253, 509
};

GeoMesh::GeoMesh(const SystemBody *body) : mSystemBody(body)
{
	mTerrain = Terrain::InstanceTerrain(body);

	for (int i=0; i<NUM_PATCHES; i++) {
		mGeoPatches[i] = nullptr;
	}

	s_allGeospheres.push_back(this);
}

GeoMesh::~GeoMesh()
{
	// update thread should not be able to access us now, so we can safely continue to delete
	assert(std::count(s_allGeospheres.begin(), s_allGeospheres.end(), this) == 1);
	s_allGeospheres.erase(std::find(s_allGeospheres.begin(), s_allGeospheres.end(), this));

	for (int i=0; i<NUM_PATCHES; i++) {
		delete mGeoPatches[i];
		mGeoPatches[i] = nullptr;
	}

	delete mTerrain;
}

void GeoMesh::Update(const vector3f &campos)
{
	if(nullptr==mGeoPatches[0]) {
		BuildFirstPatches();
	} else if(mSplitRequestDescriptions.empty()) {
		ProcessSplitResults();
		for (int i=0; i<NUM_PATCHES; i++) {
			mGeoPatches[i]->LODUpdate(campos);
		}
	} else {
		ProcessSplitRequests();
	}
}

void GeoMesh::Render(const matrix4x4f &ViewMatrix, const matrix4x4f &ModelMatrix, const matrix4x4f &MVP, 
					   Graphics::Renderer *renderer, const vector3f& campos, const float radius, const float scale)
{
	// setup the basics for the patch shader,
	// individual patches will change settings to match their own parameters
	sPatchContext->UsePatchShader(ViewMatrix, ModelMatrix, MVP);

	Render(renderer, campos, radius, scale);
}

static const float g_ambient[4] = { 0, 0, 0, 1.0 };

static void DrawAtmosphereSurface(Graphics::Renderer *renderer,
	const vector3f &campos, const float rad, Graphics::Material *mat)
{
	const int LAT_SEGS = 20;
	const int LONG_SEGS = 20;
	vector3f yaxis = campos.Normalized();
	vector3f zaxis = vector3f(1.0,0.0,0.0).Cross(yaxis).Normalized();
	vector3f xaxis = yaxis.Cross(zaxis);
	const matrix4x4f m = matrix4x4f::MakeRotMatrix(xaxis, yaxis, zaxis).InverseOf();

	glPushMatrix();
	glScalef(rad, rad, rad);
	glMultMatrixf(&m[0]);

	// what is this? Well, angle to the horizon is:
	// acos(planetRadius/viewerDistFromSphereCentre)
	// and angle from this tangent on to atmosphere is:
	// acos(planetRadius/atmosphereRadius) ie acos(1.0/1.01244blah)
	double endAng = acos(1.0/campos.Length())+acos(1.0/rad);
	double latDiff = endAng / double(LAT_SEGS);

	double rot = 0.0;
	float sinCosTable[LONG_SEGS+1][2];
	for (int i=0; i<=LONG_SEGS; i++, rot += 2.0*M_PI/double(LONG_SEGS)) {
		sinCosTable[i][0] = float(sin(rot));
		sinCosTable[i][1] = float(cos(rot));
	}

	/* Tri-fan above viewer */
	Graphics::VertexArray va(Graphics::ATTRIB_POSITION);
	va.Add(vector3f(0.f, 1.f, 0.f));
	for (int i=0; i<=LONG_SEGS; i++) {
		va.Add(vector3f(
			sin(latDiff)*sinCosTable[i][0],
			cos(latDiff),
			-sin(latDiff)*sinCosTable[i][1]));
	}
	renderer->DrawTriangles(&va, mat, Graphics::TRIANGLE_FAN);

	/* and wound latitudinal strips */
	double lat = latDiff;
	for (int j=1; j<LAT_SEGS; j++, lat += latDiff) {
		Graphics::VertexArray v(Graphics::ATTRIB_POSITION);
		float cosLat = cos(lat);
		float sinLat = sin(lat);
		float cosLat2 = cos(lat+latDiff);
		float sinLat2 = sin(lat+latDiff);
		for (int i=0; i<=LONG_SEGS; i++) {
			v.Add(vector3f(sinLat*sinCosTable[i][0], cosLat, -sinLat*sinCosTable[i][1]));
			v.Add(vector3f(sinLat2*sinCosTable[i][0], cosLat2, -sinLat2*sinCosTable[i][1]));
		}
		renderer->DrawTriangles(&v, mat, Graphics::TRIANGLE_STRIP);
	}

	glPopMatrix();
}

void GeoMesh::Render(Graphics::Renderer *renderer, const vector3f& campos, const float radius, const float scale)
{
	glPushMatrix();
	glTranslated(-campos.x, -campos.y, -campos.z);
	Graphics::Frustum frustum = Graphics::Frustum::FromGLState();

	// no frustum test of entire geosphere, since Space::Render does this
	// for each body using its GetBoundingRadius() value

	//First draw - create materials (they do not change afterwards)
	if (!m_surfaceMaterial.Valid())
		SetUpMaterials();

	if (Graphics::AreShadersEnabled()) {
		matrix4x4d modelMatrix;
		glGetDoublev (GL_MODELVIEW_MATRIX, &modelMatrix[0]);

		//Update material parameters
		//XXX no need to calculate AP every frame
		m_atmosphereParameters = mSystemBody->CalcAtmosphereParams();
		m_atmosphereParameters.center = modelMatrix * vector3d(0.0, 0.0, 0.0);
		m_atmosphereParameters.planetRadius = radius;
		m_atmosphereParameters.scale = scale;

		m_surfaceMaterial->specialParameter0 = &m_atmosphereParameters;

		if (m_atmosphereParameters.atmosDensity > 0.0) {
			m_atmosphereMaterial->specialParameter0 = &m_atmosphereParameters;

			renderer->SetBlendMode(Graphics::BLEND_ALPHA_ONE);
			renderer->SetDepthWrite(false);
			// make atmosphere sphere slightly bigger than required so
			// that the edges of the pixel shader atmosphere jizz doesn't
			// show ugly polygonal angles
			DrawAtmosphereSurface(renderer, campos, m_atmosphereParameters.atmosRadius*1.01, m_atmosphereMaterial.Get());
			renderer->SetDepthWrite(true);
			renderer->SetBlendMode(Graphics::BLEND_SOLID);
		}
	}
	glPopMatrix();

	Color ambient;
	Color &emission = m_surfaceMaterial->emissive;

	// save old global ambient
	const Color oldAmbient = renderer->GetAmbientColor();

	if ((mSystemBody->GetSuperType() == SystemBody::SUPERTYPE_STAR) || (mSystemBody->type == SystemBody::TYPE_BROWN_DWARF)) {
		// stars should emit light and terrain should be visible from distance
		ambient.r = ambient.g = ambient.b = 0.2f;
		ambient.a = 1.0f;
		emission.r = StarSystem::starRealColors[mSystemBody->type][0];
		emission.g = StarSystem::starRealColors[mSystemBody->type][1];
		emission.b = StarSystem::starRealColors[mSystemBody->type][2];
		emission.a = 1.f;
	}

	else {
		// give planet some ambient lighting if the viewer is close to it
		double camdist = campos.Length();
		camdist = 0.1 / (camdist*camdist);
		// why the fuck is this returning 0.1 when we are sat on the planet??
		// JJ: Because campos is relative to a unit-radius planet - 1.0 at the surface
		// XXX oh well, it is the value we want anyway...
		ambient.r = ambient.g = ambient.b = float(camdist);
		ambient.a = 1.0f;
	}

	renderer->SetAmbientColor(ambient);
	// this is pretty much the only place where a non-renderer is allowed to call Apply()
	// to be removed when someone rewrites terrain
	m_surfaceMaterial->Apply();

	for (int i=0; i<6; i++) {
		mGeoPatches[i]->Render(campos, frustum);
	}

	m_surfaceMaterial->Unapply();

	renderer->SetAmbientColor(oldAmbient);
}

void GeoMesh::SetUpMaterials()
{
	// Request material for this star or planet, with or without
	// atmosphere. Separate material for surface and sky.
	Graphics::MaterialDescriptor surfDesc;
	surfDesc.effect = Graphics::EFFECT_GEOSPHERE_TERRAIN;
	if ((mSystemBody->type == SystemBody::TYPE_BROWN_DWARF) ||
		(mSystemBody->type == SystemBody::TYPE_STAR_M)) {
		//dim star (emits and receives light)
		surfDesc.lighting = true;
		surfDesc.atmosphere = false;
	}
	else if (mSystemBody->GetSuperType() == SystemBody::SUPERTYPE_STAR) {
		//normal star
		surfDesc.lighting = false;
		surfDesc.atmosphere = false;
	} else {
		//planetoid with or without atmosphere
		const SystemBody::AtmosphereParameters ap(mSystemBody->CalcAtmosphereParams());
		surfDesc.lighting = true;
		surfDesc.atmosphere = (ap.atmosDensity > 0.0);
	}
	m_surfaceMaterial.Reset(Pi::renderer->CreateMaterial(surfDesc));

	//Shader-less atmosphere is drawn in Planet
	if (Graphics::AreShadersEnabled()) {
		Graphics::MaterialDescriptor skyDesc;
		skyDesc.effect = Graphics::EFFECT_GEOSPHERE_SKY;
		skyDesc.lighting = true;
		m_atmosphereMaterial.Reset(Pi::renderer->CreateMaterial(skyDesc));
	}
}


bool GeoMesh::AddSplitRequest(SSplitRequestDescription *desc)
{
	if(mSplitRequestDescriptions.size()<MAX_SPLIT_REQUESTS) {
		mSplitRequestDescriptions.push_back(desc);
		return true;
	}
	return false;
}

void GeoMesh::ProcessSplitRequests()
{
	std::deque<SSplitRequestDescription*>::const_iterator iter = mSplitRequestDescriptions.begin();
	while (iter!=mSplitRequestDescriptions.end())
	{
		const SSplitRequestDescription* srd = (*iter);

		const vector3f v01	= (srd->v0+srd->v1).Normalized();
		const vector3f v12	= (srd->v1+srd->v2).Normalized();
		const vector3f v23	= (srd->v2+srd->v3).Normalized();
		const vector3f v30	= (srd->v3+srd->v0).Normalized();
		const vector3f cn	= (srd->centroid).Normalized();

		// 
		const vector3f vecs[4][4] = {
			{srd->v0,	v01,		cn,			v30},
			{v01,		srd->v1,	v12,		cn},
			{cn,		v12,		srd->v2,	v23},
			{v30,		cn,			v23,		srd->v3}
		};

		SSplitResult *sr = new SSplitResult(srd->patchID.GetPatchFaceIdx(), srd->depth);
		for (int i=0; i<4; i++)
		{
			// Now we need to create the texture which will contain the heightmap. 
			GLuint texID;
			glGenTextures(1, &texID);
			//checkGLError();
 
			// Bind the newly created texture : all future texture functions will modify this texture
			glBindTexture(GL_TEXTURE_2D, texID);
			//checkGLError();
 
			// Create the texture itself without any data
			glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, 
				sPatchContext->fboWidth(), sPatchContext->fboWidth(), 
				0, GL_LUMINANCE, GL_FLOAT, nullptr);
			//checkGLError();
		
			// Bad filtering needed
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
			//checkGLError();

			// render the heightmap to a framebuffer.
			sPatchContext->renderHeightmap(vecs[i][0], vecs[i][1], vecs[i][2], vecs[i][3], texID);

			sr->addResult(texID, vecs[i][0], vecs[i][1], vecs[i][2], vecs[i][3], srd->patchID.NextPatchID(srd->depth+1, i));
		}

		// store result
		mSplitResult.push_back( sr );

		// cleanup after ourselves
		delete srd;

		// next!
		++iter;
	}
	mSplitRequestDescriptions.clear();
}

void GeoMesh::ProcessSplitResults()
{
	std::deque<SSplitResult*>::const_iterator iter = mSplitResult.begin();
	while(iter!=mSplitResult.end())
	{
		// finally pass SplitResults
		const SSplitResult *psr = (*iter);

		const int32_t faceIdx = psr->face;
		mGeoPatches[faceIdx]->ReceiveHeightmaps(psr);

		// tidyup
		delete psr;

		// Next!
		++iter;
	}
	mSplitResult.clear();
}


static const int geo_sphere_edge_friends[6][4] = {
	{ 3, 4, 1, 2 },
	{ 0, 4, 5, 2 },
	{ 0, 1, 5, 3 },
	{ 0, 2, 5, 4 },
	{ 0, 3, 5, 1 },
	{ 1, 4, 3, 2 }
};

void GeoMesh::BuildFirstPatches()
{
	// generate root face patches of the cube/sphere
	static const vector3f p1 = (vector3f( 1, 1, 1)).Normalized();
	static const vector3f p2 = (vector3f(-1, 1, 1)).Normalized();
	static const vector3f p3 = (vector3f(-1,-1, 1)).Normalized();
	static const vector3f p4 = (vector3f( 1,-1, 1)).Normalized();
	static const vector3f p5 = (vector3f( 1, 1,-1)).Normalized();
	static const vector3f p6 = (vector3f(-1, 1,-1)).Normalized();
	static const vector3f p7 = (vector3f(-1,-1,-1)).Normalized();
	static const vector3f p8 = (vector3f( 1,-1,-1)).Normalized();

	const uint64_t maxShiftDepth = GeoPatchID::MAX_SHIFT_DEPTH;

	mGeoPatches[0] = new GeoPatch(*sPatchContext, this, p1, p2, p3, p4, 0, (0i64 << maxShiftDepth));
	mGeoPatches[1] = new GeoPatch(*sPatchContext, this, p4, p3, p7, p8, 0, (1i64 << maxShiftDepth));
	mGeoPatches[2] = new GeoPatch(*sPatchContext, this, p1, p4, p8, p5, 0, (2i64 << maxShiftDepth));
	mGeoPatches[3] = new GeoPatch(*sPatchContext, this, p2, p1, p5, p6, 0, (3i64 << maxShiftDepth));
	mGeoPatches[4] = new GeoPatch(*sPatchContext, this, p3, p2, p6, p7, 0, (4i64 << maxShiftDepth));
	mGeoPatches[5] = new GeoPatch(*sPatchContext, this, p8, p7, p6, p5, 0, (5i64 << maxShiftDepth));

	for (int i=0; i<NUM_PATCHES; i++) {
		for (int j=0; j<4; j++) {
			mGeoPatches[i]->OnEdgeFriendChanged(j, mGeoPatches[geo_sphere_edge_friends[i][j]]);
		}
	}
}

void GeoMesh::Init()
{
	sPatchContext.Reset(new GeoPatchContext(detail_edgeLen[Pi::detail.planets > 4 ? 4 : Pi::detail.planets]));
	assert(sPatchContext->edgeLen() <= GEOPATCH_MAX_EDGELEN);
}

void GeoMesh::Uninit()
{
	assert (sPatchContext.Unique());
	sPatchContext.Reset();
}

void GeoMesh::OnChangeDetailLevel()
{
	sPatchContext.Reset(new GeoPatchContext(detail_edgeLen[Pi::detail.planets > 4 ? 4 : Pi::detail.planets]));
	assert(sPatchContext->edgeLen() <= GEOPATCH_MAX_EDGELEN);

	// reinit the geosphere terrain data
	for(std::vector<GeoMesh*>::iterator i = s_allGeospheres.begin(); i != s_allGeospheres.end(); ++i) 
	{
		for (int p=0; p<6; p++) {
			// delete patches
			if ((*i)->mGeoPatches[p]) {
				delete (*i)->mGeoPatches[p];
				(*i)->mGeoPatches[p] = NULL;
			}
		}

		// reinit the terrain with the new settings
		delete (*i)->mTerrain;
		(*i)->mTerrain = Terrain::InstanceTerrain((*i)->mSystemBody);
	}
}