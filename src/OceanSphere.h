#ifndef _OCEANSPHERE_H
#define _OCEANSPHERE_H

#include <SDL_stdinc.h>

#include "vector3.h"
#include "mtrand.h"
#include "galaxy/StarSystem.h"
#include "graphics/Material.h"
#include "terrain/Terrain.h"

namespace Graphics { class Renderer; }
class SystemBody;
class OceanPatch;
class OceanPatchContext;
class OceanSphere {
public:
	OceanSphere(const SystemBody *body);
	~OceanSphere();
	int UpdateLODs(const vector3d campos);
	void Render(Graphics::Renderer *r, const vector3d campos, const float radius, const float scale);
	void BuildFirstPatches();
	void DestroyPatches();
	inline double GetHeight(const vector3d &p) const {
		const double h = 0.00001;//m_terrain->GetHeight(p);
		s_vtxGenCount++;
#ifdef DEBUG
		// XXX don't remove this. Fix your fractals instead
		// Fractals absolutely MUST return heights >= 0.0 (one planet radius)
		// otherwise atmosphere and other things break.
		assert(h >= 0.0);
#endif /* DEBUG */
		return h;
	}
	void SetUpMaterials();
	friend class OceanPatch;
	static void Init();
	static void Uninit();
	static void OnChangeDetailLevel();
	// in sbody radii
	double GetMaxFeatureHeight() const { return m_terrain->GetMaxHeight(); }
	static int GetVtxGenCount() { return s_vtxGenCount; }
	static void ClearVtxGenCount() { s_vtxGenCount = 0; }

private:
	OceanPatch *m_patches[6];
	const SystemBody *m_sbody;

	/* all variables for GetHeight(), GetColor() */
	Terrain *m_terrain;

	///////////////////////////
	std::list<GLuint> m_vbosToDestroy;
	SDL_mutex *m_vbosToDestroyLock;
	void AddVBOToDestroy(GLuint vbo);
	void DestroyVBOs();
	//////////////////////////////

	inline vector3d GetColor(const vector3d &p, double height, const vector3d &norm) {
		return vector3d(1.0,0.0,0.0);//m_terrain->GetColor(p, height, norm);
	}

	static int s_vtxGenCount;

	static RefCountedPtr<OceanPatchContext> s_patchContext;

	
	ScopedPtr<Graphics::Material> m_oceanMaterial;
	//special parameters for shaders
	SystemBody::AtmosphereParameters m_atmosphereParameters;
};

#endif /* _OCEANSPHERE_H */
