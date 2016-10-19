#ifndef _GEORING_H
#define _GEORING_H

#include "vector3.h"
#include "Camera.h"
#include "terrain/Terrain.h"

extern int GEOPATCH_EDGELEN;
#define ATMOSPHERE_RADIUS 1.015

namespace Graphics { 
	class Renderer; 
	class RenderState;
	class Material;
}
class SystemBody;
class GeoPlate;
class GeoPlateHull;
class GeoPlateWall;
class GeoPatchContext;

class SQuadPlateRequest;
class SQuadPlateResult;
class SSinglePlateResult;

class GeoRing {
public:
	GeoRing(const SystemBody *body);
	~GeoRing();

	void Update();
	void Render(Graphics::Renderer *renderer, const matrix4x4d &modelView, vector3d campos, const float radius, const std::vector<Camera::Shadow> &shadows);
	//double GetDistFromSurface(const vector3d p);
	double GetRingWidth() const { return mRingWidth; }
	friend class GeoPlate;
	friend class GeoPlateHull;
	friend class GeoPlateWall;
#if WITH_OBJECTVIEWER
	friend class ObjectViewerView;
#endif
	static void Init();
	static void UpdateAllGeoRings();
	static void OnChangeDetailLevel();
	static bool OnAddQuadPlateResult(const SystemPath &path, SQuadPlateResult *res);
	static bool OnAddSinglePlateResult(const SystemPath &path, SSinglePlateResult *res);
	
	bool AddQuadPlateResult(SQuadPlateResult *res);
	void AddQuadPlateRequest(double, SQuadPlateRequest*, GeoPlate*);
	bool AddSinglePlateResult(SSinglePlateResult *res);
	void ProcessSplitResults();
	void ProcessQuadPlateRequests();

	void GetAtmosphereFlavor(Color *outColor, double *outDensity) const {
		m_sbody->GetAtmosphereFlavor(outColor, outDensity);
	}
	// in sbody radii
	double GetMaxFeatureHeight() const { return m_terrain->GetMaxHeight(); }

	inline double GetHeight(vector3d p) {
		/*return 0.0;*/
		const double h = m_terrain->GetHeight(p);
		// XXX don't remove this. Fix your fractals instead
		// Fractals absolutely MUST return heights >= 0.0 (one planet radius)
		// otherwise atmosphere and other things break.
		assert(h >= 0.0);
		return h;
	}

	inline Sint32 GetMaxDepth() const { return m_maxDepth; }
	inline const SystemBody *GetSystemBody() const { return m_sbody; }
	Terrain* GetTerrain() const { return m_terrain.Get(); }

	struct MaterialParameters {
		SystemBody::AtmosphereParameters atmosphere;
		std::vector<Camera::Shadow> shadows;
		Sint32 patchDepth;
		Sint32 maxPatchDepth;
	};

	Graphics::RenderState* GetSurfRenderState() const { return m_surfRenderState; }
	RefCountedPtr<Graphics::Material> GetSurfaceMaterial() const { return m_surfaceMaterial; }
	MaterialParameters& GetMaterialParameters() { return m_materialParameters; }

private:
	// only called from fishy thread
	void _UpdateLODs();
	void BuildFirstPatches();
	void SetUpMaterials();
	void CalculateMaxPatchDepth();
	typedef std::vector<GeoPlate*>::iterator PlateIter;
	std::vector<GeoPlate*>		m_plates;
	std::vector<GeoPlateHull*>	m_hull;
	std::vector<GeoPlateWall*>	m_wallInner;
	std::vector<GeoPlateWall*>	m_wallOuter;
	float m_diffColor[4], m_ambColor[4];
	const SystemBody *m_sbody;
	RefCountedPtr<Terrain> m_terrain;
	double mRingWidth;

	Sint32 m_maxDepth;

	RefCountedPtr<Graphics::Texture> m_texHi;
	RefCountedPtr<Graphics::Texture> m_texLo;

	enum EGSInitialisationStage {
		eBuildFirstPatches=0,
		eRequestedFirstPatches,
		eReceivedFirstPatches,
		eDefaultUpdateState
	};
	EGSInitialisationStage m_initStage;

	Graphics::RenderState *m_surfRenderState;
	RefCountedPtr<Graphics::Material> m_surfaceMaterial;
	//special parameters for shaders
	MaterialParameters m_materialParameters;
	
	bool m_hasTempCampos;
	vector3d m_tempCampos;
	Graphics::Frustum m_tempFrustum;
	//////////////////////////////

	struct TDistanceRequest {
		TDistanceRequest(double dist, SQuadPlateRequest *pRequest, GeoPlate *pRequester) :
			mDistance(dist), mpRequest(pRequest), mpRequester(pRequester) {}
		double mDistance;
		SQuadPlateRequest *mpRequest;
		GeoPlate *mpRequester;
	};
	std::deque<TDistanceRequest> mQuadPlateRequests;

	static RefCountedPtr<GeoPatchContext> s_patchContext;
	static const uint32_t MAX_SPLIT_OPERATIONS = 256;
	std::deque<SQuadPlateResult*> mQuadPlateResults;
	std::deque<SSinglePlateResult*> mSinglePlateResults;

	inline vector3d GetColor(const vector3d &p, double height, const vector3d &norm) {
		//return vector3d(0.5, 0.5, 0.5);
		return m_terrain->GetColor(p, height, norm);
	}
};

#endif /* _GEORING_H */
