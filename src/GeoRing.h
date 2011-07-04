#ifndef _GEORING_H
#define _GEORING_H

#include "vector3.h"
#include "mtrand.h"
#include "GeoRingStyle.h"

extern int GEOPATCH_EDGELEN;
#define ATMOSPHERE_RADIUS 1.015

class SBody;
class GeoPlate;
class GeoPlateHull;
class GeoRing {
public:
	GeoRing(const SBody *body);
	~GeoRing();
	void Render(vector3d campos, const float radius, const float scale);
	inline double GetHeight(vector3d p) {
		const double h = m_style.GetHeight(p);
#ifdef DEBUG
		// XXX don't remove this. Fix your fractals instead
		// Fractals absolutely MUST return heights >= 0.0 (one planet radius)
		// otherwise atmosphere and other things break.
		assert(h >= 0.0);
#endif /* DEBUG */
		return h;
	}
	// only called from fishy thread
	void _UpdateLODs();
	friend class GeoPlate;
	friend class GeoPlateHull;
#if OBJECTVIEWER
	friend class OrbitalViewerView;
#endif /* DEBUG */
	static void Init();
	static void OnChangeDetailLevel();
	void GetAtmosphereFlavor(Color *outColor, double *outDensity) const {
		m_style.GetAtmosphereFlavor(outColor, outDensity);
	}
	// in sbody radii
	double GetMaxFeatureHeight() const { return m_style.GetMaxHeight(); }
private:
	void BuildFirstPatches(const int numSegments = 16);
	std::vector<GeoPlate*>		m_plates;
	std::vector<GeoPlateHull*>	m_hull;
	float m_diffColor[4], m_ambColor[4];
	const SBody *m_sbody;

	/* all variables for GetHeight(), GetColor() */
	GeoRingStyle m_style;

	///////////////////////////
	// threading rubbbbbish
	// update thread can't do it since only 1 thread can molest opengl
	static int UpdateLODThread(void *data) __attribute((noreturn));
	std::list<GLuint> m_vbosToDestroy;
	SDL_mutex *m_vbosToDestroyLock;
	void AddVBOToDestroy(GLuint vbo);
	void DestroyVBOs();
	
	vector3d m_tempCampos;
	int m_runUpdateThread;
	//////////////////////////////

	inline vector3d GetColor(const vector3d &p, double height, const vector3d &norm) {
		return m_style.GetColor(p, height, norm);
	}
};

#endif /* _GEORING_H */
