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
		++s_vtxGenCount;
		return 0.0;
	}
	void SetUpMaterials();
	friend class OceanPatch;
	static void Init();
	static void Uninit();
	static void OnChangeDetailLevel();
	static int GetVtxGenCount() { return s_vtxGenCount; }
	static void ClearVtxGenCount() { s_vtxGenCount = 0; }

private:
	OceanPatch *m_patches[6];
	const SystemBody *m_sbody;

	/* all variables for GetHeight(), GetColor() */
	//Terrain *m_terrain;

	///////////////////////////
	std::list<GLuint> m_vbosToDestroy;
	SDL_mutex *m_vbosToDestroyLock;
	void AddVBOToDestroy(GLuint vbo);
	void DestroyVBOs();
	//////////////////////////////

	static int s_vtxGenCount;

	static RefCountedPtr<OceanPatchContext> s_patchContext;
	
	ScopedPtr<Graphics::Material> m_oceanMaterial;
	//special parameters for shaders
	SystemBody::AtmosphereParameters m_atmosphereParameters;
};

#endif /* _OCEANSPHERE_H */
