// Copyright © 2008-2014 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _GEOSPHERE_H
#define _GEOSPHERE_H

#include <SDL_stdinc.h>

#include "vector3.h"
#include "Random.h"
#include "Camera.h"
#include "galaxy/StarSystem.h"
#include "graphics/RenderState.h"
#include "graphics/Material.h"
#include "terrain/Terrain.h"
#include "GeoPatchID.h"
#include "BaseSphere.h"

#include <deque>

namespace Graphics { class Renderer; }
class SystemBody;
class GeoPatch;
class GeoPatchContext;
class SQuadSplitRequest;
class SQuadSplitResult;
class SSingleSplitResult;

//------------------------------------------------------------

// Three dimensional sphere (subdivided icosahedron) with normals
// and spherical texture coordinates.
class Icosahedron {
public:
	//subdivisions must be 0-4
	Icosahedron(Graphics::Renderer*, RefCountedPtr<Graphics::Material> material, Graphics::RenderState*, int subdivisions=0, float scale=1.f);
	virtual void Draw(Graphics::Renderer *r);

	RefCountedPtr<Graphics::Material> GetMaterial() const { return m_material; }

private:
	std::unique_ptr<Graphics::VertexBuffer> m_vertexBuffer;
	std::unique_ptr<Graphics::IndexBuffer> m_indexBuffer;
	RefCountedPtr<Graphics::Material> m_material;
	Graphics::RenderState *m_renderState;

	//std::unique_ptr<Surface> m_surface;
	//add a new vertex, return the index
	int AddVertex(Graphics::VertexArray&, const vector3f &v, const vector3f &n);
	//add three vertex indices to form a triangle
	void AddTriangle(std::vector<Uint32>&, int i1, int i2, int i3);
	void Subdivide(Graphics::VertexArray&, std::vector<Uint32>&,
		const matrix4x4f &trans, const vector3f &v1, const vector3f &v2, const vector3f &v3,
		const int i1, const int i2, const int i3, const int depth);
};
//------------------------------------------------------------

#define NUM_PATCHES 6

class GeoSphere : public BaseSphere {
public:
	GeoSphere(const SystemBody *body);
	virtual ~GeoSphere();

	virtual void Update();
	virtual void Render(Graphics::Renderer *renderer, const matrix4x4d &modelView, vector3d campos, const float radius, const float scale, const std::vector<Camera::Shadow> &shadows);

	virtual double GetHeight(const vector3d &p) const {
		const double h = m_terrain->GetHeight(p);
#ifdef DEBUG
		// XXX don't remove this. Fix your fractals instead
		// Fractals absolutely MUST return heights >= 0.0 (one planet radius)
		// otherwise atmosphere and other things break.
		if (h < 0.0) {
			Output("GetHeight({ %f, %f, %f }) returned %f\n", p.x, p.y, p.z, h);
			m_terrain->DebugDump();
			assert(h >= 0.0);
		}
#endif /* DEBUG */
		return h;
	}
	
	static void Init();
	static void Uninit();
	static void UpdateAllGeoSpheres();
	static void OnChangeDetailLevel();
	static bool OnAddQuadSplitResult(const SystemPath &path, SQuadSplitResult *res);
	static bool OnAddSingleSplitResult(const SystemPath &path, SSingleSplitResult *res);
	// in sbody radii
	virtual double GetMaxFeatureHeight() const { return m_terrain->GetMaxHeight(); }

	bool AddQuadSplitResult(SQuadSplitResult *res);
	bool AddSingleSplitResult(SSingleSplitResult *res);
	void ProcessSplitResults();

	virtual void Reset();

	inline Sint32 GetMaxDepth() const { return m_maxDepth; }

private:
	void BuildFirstPatches();
	void CalculateMaxPatchDepth();
	inline vector3d GetColor(const vector3d &p, double height, const vector3d &norm) const {
		return m_terrain->GetColor(p, height, norm);
	}

	std::unique_ptr<GeoPatch> m_patches[6];

	static const uint32_t MAX_SPLIT_OPERATIONS = 128;
	std::deque<SQuadSplitResult*> mQuadSplitResults;
	std::deque<SSingleSplitResult*> mSingleSplitResults;

	bool m_hasTempCampos;
	vector3d m_tempCampos;

	static RefCountedPtr<GeoPatchContext> s_patchContext;

	virtual void SetUpMaterials();

	enum EGSInitialisationStage {
		eBuildFirstPatches=0,
		eRequestedFirstPatches,
		eReceivedFirstPatches,
		eDefaultUpdateState
	};
	EGSInitialisationStage m_initStage;

	Sint32 m_maxDepth;
};

#endif /* _GEOSPHERE_H */
