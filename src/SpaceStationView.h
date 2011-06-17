#ifndef _SPACESTATIONVIEW_H
#define _SPACESTATIONVIEW_H

#include "libs.h"
#include "gui/Gui.h"
#include "View.h"
#include "GenericChatForm.h"

class StationViewShipView;
class StationShipPaintView;

class SpaceStationView: public View {
public:
	SpaceStationView();
	virtual ~SpaceStationView() {}
	virtual void Update();
	virtual void Draw3D();
	virtual void OnSwitchTo();
	virtual void JumpToForm(GenericChatForm *f);
	// friends
	friend class StationViewShipView;
	friend class StationShipPaintView;
private:
	// hack so StationViewShipView can draw its 3d shit
	sigc::signal<void> onDraw3D;
	GenericChatForm *m_baseSubView;
	GenericChatForm *m_jumpToForm;
};

#endif /* _SPACESTATIONVIEW_H */
