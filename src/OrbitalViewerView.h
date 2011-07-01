#ifndef _ORBITALVIEWERVIEW_H
#define _ORBITALVIEWERVIEW_H

#include "libs.h"
#include "gui/Gui.h"
#include "View.h"

#if ORBITALVIEWER

class Body;

class OrbitalViewerView: public View {
public:
	OrbitalViewerView();
	virtual void Update();
	virtual void Draw3D();
	virtual void OnSwitchTo();
private:
	float viewingDist;
	Gui::Label *m_infoLabel;
	const Body* lastTarget;
	matrix4x4d m_camRot;

	Gui::TextEntry *m_sbodyMass;
	Gui::TextEntry *m_sbodyRadius;
	Gui::TextEntry *m_sbodySeed;
	Gui::TextEntry *m_sbodyVolatileGas;
	Gui::TextEntry *m_sbodyVolatileLiquid;
	Gui::TextEntry *m_sbodyVolatileIces;
	Gui::TextEntry *m_sbodyLife;
	Gui::TextEntry *m_sbodyVolcanicity;
	Gui::TextEntry *m_sbodyMetallicity;
	void OnChangeGeoRingStyle();
};

#endif

#endif /* _ORBITALVIEWERVIEW_H */
