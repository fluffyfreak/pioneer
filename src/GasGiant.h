// Copyright © 2008-2014 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _GASGIANT_H
#define _GASGIANT_H

#include <SDL_stdinc.h>

#include "vector3.h"
#include "Random.h"
#include "Camera.h"
#include "graphics/Drawables.h"
#include "graphics/Material.h"
#include "terrain/Terrain.h"
#include "BaseSphere.h"
#include "JobQueue.h"
#include "GPUJobQueue.h"

#include <deque>

namespace Graphics { class Renderer; }
class SystemBody;
class GasGiant;
class GasPatch;
class GasPatchContext;
namespace { 
	class STextureFaceResult; 
	class STextureFaceGPUResult;
}

#define NUM_PATCHES 6

class GasGiant : public BaseSphere {
public:
	GasGiant(const SystemBody *body);
	virtual ~GasGiant();

	virtual void Update();
	virtual void Render(Graphics::Renderer *renderer, const matrix4x4d &modelView, vector3d campos, const float radius, const float scale, const std::vector<Camera::Shadow> &shadows);

	virtual double GetHeight(const vector3d &p) const { return 0.0; }

	// in sbody radii
	virtual double GetMaxFeatureHeight() const { return 0.0; }

	virtual void Reset() {};

	static bool OnAddTextureFaceResult(const SystemPath &path, STextureFaceResult *res);
	static bool OnAddTextureFaceResult(const SystemPath &path, STextureFaceGPUResult *res);
	static void Init();
	static void Uninit();
	static void UpdateAllGasGiants();

	static void CreateRenderTarget(const Uint16 width, const Uint16 height);
	static void SetRenderTargetCubemap(const Uint32, Graphics::Texture*);
	static void BeginRenderTarget();
	static void EndRenderTarget();

private:
	void BuildFirstPatches();
	void GenerateTexture();
	bool AddTextureFaceResult(STextureFaceResult *res);
	bool AddTextureFaceResult(STextureFaceGPUResult *res);

	static RefCountedPtr<GasPatchContext> s_patchContext;

	static Graphics::RenderTarget *s_renderTarget;
	static Graphics::RenderState *s_quadRenderState;

	std::unique_ptr<GasPatch> m_patches[NUM_PATCHES];

	bool m_hasTempCampos;
	vector3d m_tempCampos;

	virtual void SetUpMaterials();
	RefCountedPtr<Graphics::Texture> m_surfaceTextureSmall;
	RefCountedPtr<Graphics::Texture> m_surfaceTexture;
	RefCountedPtr<Graphics::Texture> m_builtTexture;
	
	std::unique_ptr<Color[]> m_jobColorBuffers[NUM_PATCHES];
	JobHandle m_job[NUM_PATCHES];
	bool m_hasJobRequest[NUM_PATCHES];

	GPUJobHandle m_gpuJob[NUM_PATCHES];
	bool m_hasGpuJobRequest[NUM_PATCHES];

	float m_timeDelay;
};

#endif /* _GASGIANT_H */
