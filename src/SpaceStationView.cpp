#include "SpaceStation.h"
#include "SpaceStationView.h"
#include "Pi.h"
#include "Player.h"
#include "WorldView.h"
#include "ShipFlavour.h"
#include "ShipCpanel.h"
#include "CommodityTradeWidget.h"
#include "GenericChatForm.h"
#include "LuaChatForm.h"
#include "PoliceChatForm.h"
#include "LmrModel.h"
#include "utils.h"
#include "PaintJob.h"

////////////////////////////////////////////////////////////////////

class StationCommoditiesView: public GenericChatForm {
public:
	StationCommoditiesView();
	virtual void ShowAll();
private:
	void OnClickBuy(int commodity_type) {
		if (m_station->SellTo(Pi::player, Equip::Type(commodity_type), true)) {
			Pi::cpan->MsgLog()->Message("", stringf(512, "You have bought 1t of %s.", EquipType::types[commodity_type].name));
		}
		m_commodityTradeWidget->UpdateStock(commodity_type);
		UpdateBaseDisplay();
	}
	void OnClickSell(int commodity_type) {
		if (m_station->BuyFrom(Pi::player, Equip::Type(commodity_type), true)) {
			Pi::cpan->MsgLog()->Message("", stringf(512, "You have sold 1t of %s.", EquipType::types[commodity_type].name));
		}
		m_commodityTradeWidget->UpdateStock(commodity_type);
		UpdateBaseDisplay();
	}
	SpaceStation *m_station;
	CommodityTradeWidget *m_commodityTradeWidget;
};

StationCommoditiesView::StationCommoditiesView(): GenericChatForm()
{
	SetTransparency(false);
}

void StationCommoditiesView::ShowAll()
{
	ReInit();

	m_station = Pi::player->GetDockedWith();
	assert(m_station);
	SetTransparency(false);
	{
		char buf[256];
		snprintf(buf, sizeof(buf), "Welcome to %s commodities market", m_station->GetLabel().c_str());
		SetTitle(buf);
	}

	Gui::Button *backButton = new Gui::SolidButton();
	backButton->onClick.connect(sigc::mem_fun(this, &StationCommoditiesView::Close));
	Add(backButton,680,500);
	Add(new Gui::Label("Go back"), 700, 500);

	m_commodityTradeWidget = new CommodityTradeWidget(m_station);
	m_commodityTradeWidget->SetSizeRequest(470.0, 450.0);
	m_commodityTradeWidget->onClickBuy.connect(sigc::mem_fun(this, &StationCommoditiesView::OnClickBuy));
	m_commodityTradeWidget->onClickSell.connect(sigc::mem_fun(this, &StationCommoditiesView::OnClickSell));
	Add(m_commodityTradeWidget, 320, 40);

	AddBaseDisplay();
	AddVideoWidget();

	Gui::Fixed::ShowAll();
}

////////////////////////////////////////////////////////////////////

#define REMOVAL_VALUE_PERCENT 90

class StationLaserPickMount: public GenericChatForm {
public:
	StationLaserPickMount(Equip::Type t, bool doFit) {
		m_equipType = t;
		m_doFit = doFit;
	}
	virtual void ShowAll();
private:
	void SelectMount(int mount) {
		SpaceStation *station = Pi::player->GetDockedWith();
		Equip::Slot slot = EquipType::types[m_equipType].slot;
		// duplicates some code in StationShipUpgradesView
		if (m_doFit) {
			// fit
			Pi::player->m_equipment.Set(slot, mount, m_equipType);
			Pi::player->UpdateMass();
			Pi::player->SetMoney(Pi::player->GetMoney() - station->GetPrice(m_equipType));
			Pi::cpan->MsgLog()->Message("", "Fitting "+std::string(EquipType::types[m_equipType].name));
			Gui::Container *p = GetParent();
			Close();
			p->ShowAll();
		} else {
			// remove 
			const Sint64 value = station->GetPrice(m_equipType) * REMOVAL_VALUE_PERCENT / 100;
			Pi::player->m_equipment.Set(slot, mount, Equip::NONE);
			Pi::player->UpdateMass();
			Pi::player->SetMoney(Pi::player->GetMoney() + value);
			station->AddEquipmentStock(m_equipType, 1);
			Pi::cpan->MsgLog()->Message("", "Removing "+std::string(EquipType::types[m_equipType].name));
			Gui::Container *p = GetParent();
			Close();
			p->ShowAll();
		}
	}
	Equip::Type m_equipType;
	bool m_doFit;
};

void StationLaserPickMount::ShowAll()
{
	ReInit();
	SetTransparency(false);
	
	if (m_doFit) Add(new Gui::Label("Fit laser to which gun mount?"), 320, 200);
	else Add(new Gui::Label("Remove laser from which gun mount?"), 320, 200);

	Equip::Slot slot = EquipType::types[m_equipType].slot;

	int xpos = 400;
	for (int i=0; i<ShipType::GUNMOUNT_MAX; i++) {
		if (m_doFit && (Pi::player->m_equipment.Get(slot, i) != Equip::NONE)) continue;
		if ((!m_doFit) && (Pi::player->m_equipment.Get(slot, i) != m_equipType)) continue;
		Gui::Button *b = new Gui::SolidButton();
		b->onClick.connect(sigc::bind(sigc::mem_fun(this, &StationLaserPickMount::SelectMount), i));
		Add(b, float(xpos), 250);
		Add(new Gui::Label(ShipType::gunmountNames[i]), float(xpos), 270);

		xpos += 50;
	}
	Gui::Fixed::ShowAll();
}

////////////////////////////////////////////////////////////////////

class StationShipUpgradesView: public GenericChatForm {
public:
	StationShipUpgradesView();
	virtual void ShowAll();
private:
	void RemoveItem(Equip::Type t) {
		SpaceStation *station = Pi::player->GetDockedWith();
		Equip::Slot s = EquipType::types[t].slot;
		Sint64 value = station->GetPrice(t) * REMOVAL_VALUE_PERCENT / 100;
		int num = Pi::player->m_equipment.Count(s, t);
		
		if (num) {
			if ((s == Equip::SLOT_LASER) && (num > 1)) {
				/* you have a choice of mount points for lasers */
				OpenChildChatForm(new StationLaserPickMount(t, false));
			} else {
				Pi::player->m_equipment.Remove(t, 1);
				Pi::player->UpdateMass();
				Pi::player->SetMoney(Pi::player->GetMoney() + value);
				station->AddEquipmentStock(t, 1);
				Pi::cpan->MsgLog()->Message("", "Removing "+std::string(EquipType::types[t].name));
				UpdateBaseDisplay();
				ShowAll();
			}
		}
	}
	void FitItem(Equip::Type t) {
		SpaceStation *station = Pi::player->GetDockedWith();
		Equip::Slot s = EquipType::types[t].slot;
		const shipstats_t *stats = Pi::player->CalcStats();
		int freespace = Pi::player->m_equipment.FreeSpace(s);
		
		if (Pi::player->GetMoney() < station->GetPrice(t)) {
			Pi::cpan->MsgLog()->Message("", "You do not have enough money");
		} else if (stats->free_capacity < EquipType::types[t].mass) {
			Pi::cpan->MsgLog()->Message("", "There is no space on your ship");
		} else if (freespace) {
			if ((freespace > 1) && (s == Equip::SLOT_LASER)) {
				/* you have a choice of mount points for lasers */
				OpenChildChatForm(new StationLaserPickMount(t, true));
			} else {
				Pi::player->m_equipment.Add(t);
				Pi::player->UpdateMass();
				Pi::player->SetMoney(Pi::player->GetMoney() - station->GetPrice(t));
				Pi::cpan->MsgLog()->Message("", "Fitting "+std::string(EquipType::types[t].name));
				UpdateBaseDisplay();
				ShowAll();
			}
		} else {
			Pi::cpan->MsgLog()->Message("", "There is no space on your ship");
		}
	}
};

StationShipUpgradesView::StationShipUpgradesView(): GenericChatForm()
{
	SetTransparency(false);
}

void StationShipUpgradesView::ShowAll()
{
	ReInit();

	SpaceStation *station = Pi::player->GetDockedWith();
	assert(station);
	SetTransparency(false);
	SetTitle(stringf(256, "%s Shipyard", station->GetLabel().c_str()).c_str());
	
	Gui::Button *backButton = new Gui::SolidButton();
	backButton->onClick.connect(sigc::mem_fun(this, &StationShipUpgradesView::Close));
	Add(backButton,680,470);
	Add(new Gui::Label("Go back"), 700, 470);

	Gui::Fixed *fbox = new Gui::Fixed(470, 400);
	Add(fbox, 320, 40);

	Gui::VScrollBar *scroll = new Gui::VScrollBar();
	Gui::VScrollPortal *portal = new Gui::VScrollPortal(450,400);
	scroll->SetAdjustment(&portal->vscrollAdjust);
	//int GetStock(Equip::Type t) const { return m_equipmentStock[t]; }

	int NUM_ITEMS = 0;
	const float YSEP = floor(Gui::Screen::GetFontHeight() * 1.5f);
	for (int i=Equip::FIRST_SHIPEQUIP; i<=Equip::LAST_SHIPEQUIP; i++) {
		if (station->GetStock(static_cast<Equip::Type>(i))) NUM_ITEMS++;
	}

	Gui::Fixed *innerbox = new Gui::Fixed(450, NUM_ITEMS*YSEP);
	for (int i=Equip::FIRST_SHIPEQUIP, num=0; i<=Equip::LAST_SHIPEQUIP; i++) {
		Equip::Type type = static_cast<Equip::Type>(i);
		int stock = station->GetStock(type);
		if (!stock) continue;
		Gui::Label *l = new Gui::Label(EquipType::types[i].name);
		if (EquipType::types[i].description) {
			l->SetToolTip(EquipType::types[i].description);
		}
		innerbox->Add(l,0,num*YSEP);
		
		innerbox->Add(new Gui::Label(format_money(station->GetPrice(type))), 200, num*YSEP);

		innerbox->Add(new Gui::Label(format_money(REMOVAL_VALUE_PERCENT * station->GetPrice(type) / 100)),
				275, num*YSEP);
		
		innerbox->Add(new Gui::Label(stringf(64, "%dt", EquipType::types[i].mass)), 360, num*YSEP);
		
		Gui::Button *b;
		b = new Gui::SolidButton();
		b->onClick.connect(sigc::bind(sigc::mem_fun(this, &StationShipUpgradesView::FitItem), type));
		innerbox->Add(b, 400, num*YSEP);
		// only have remove button if we have this item installed
		if (Pi::player->m_equipment.Count(EquipType::types[i].slot, type )) {
			b = new Gui::SolidButton();
			innerbox->Add(b, 420, num*YSEP);
			b->onClick.connect(sigc::bind(sigc::mem_fun(this, &StationShipUpgradesView::RemoveItem), type));
		}
		num++;
	}
	innerbox->ShowAll();

	const float *col = Gui::Theme::Colors::tableHeading;
	fbox->Add((new Gui::Label("Item"))->Color(col), 0, 0);
	fbox->Add((new Gui::Label("$ to fit"))->Color(col), 200, 0);
	fbox->Add((new Gui::Label("$ for removal"))->Color(col), 275, 0);
	fbox->Add((new Gui::Label("Wt"))->Color(col), 360, 0);
	fbox->Add((new Gui::Label("Fit"))->Color(col), 400, 0);
	fbox->Add((new Gui::Label("Remove"))->Color(col), 420, 0);
	fbox->Add(portal, 0, YSEP);
	fbox->Add(scroll, 455, YSEP);
	portal->Add(innerbox);
	portal->ShowAll();
	fbox->ShowAll();
	AddBaseDisplay();
	AddVideoWidget();

	Gui::Fixed::ShowAll();
}

////////////////////////////////////////////////////////////////////

class StationViewShipView: public GenericChatForm {
public:
	StationViewShipView(int flavour_idx): GenericChatForm() {
		SpaceStation *station = Pi::player->GetDockedWith();
		m_flavourIdx = flavour_idx;
		m_flavour = station->GetShipsOnSale()[flavour_idx];
		m_lmrModel = LmrLookupModelByName(ShipType::types[m_flavour.type].lmrModelName.c_str());
		m_ondraw3dcon = Pi::spaceStationView->onDraw3D.connect(
				sigc::mem_fun(this, &StationViewShipView::Draw3D));
	}
	virtual ~StationViewShipView() {
		m_ondraw3dcon.disconnect();
	}
	virtual void ShowAll();
	void Draw3D();
private:
	void BuyShip() {
		// nasty. signal station->onShipsForSaleChanged will fire
		// and 'this' will be erased when StationShipyardView is
		// updated with new contents, so m_flavour must be kept on stack
		// alloc to avoid trouble
		ShipFlavour f = m_flavour;
		ShipFlavour _old = *Pi::player->GetFlavour();
		int cost = f.price - Pi::player->GetFlavour()->price;
		if (Pi::player->GetMoney() < cost) {
			Pi::cpan->MsgLog()->Message("", "You do not have enough money");
		} else {
			Pi::player->SetMoney(Pi::player->GetMoney() - cost);
			Pi::player->ResetFlavour(&f);
			Pi::player->m_equipment.Set(Equip::SLOT_ENGINE, 0, ShipType::types[f.type].hyperdrive);
			Pi::player->UpdateMass();
			
			SpaceStation *station = Pi::player->GetDockedWith();
			station->ReplaceShipOnSale(m_flavourIdx, &_old);
		}
	}
	ShipFlavour m_flavour;
	LmrModel *m_lmrModel;
	int m_flavourIdx;
	sigc::connection m_ondraw3dcon;
};

void StationViewShipView::Draw3D()
{
	/* XXX duplicated code in InfoView.cpp */
	LmrObjParams params = {
		{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ "HELLO" },
		{ 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 0.0f },

		{	// pColor[3]
		{ { 1.0f, 0.0f, 1.0f, 1.0f }, { 0, 0, 0 }, { 0, 0, 0 }, 0 },
		{ { 0.8f, 0.6f, 0.5f, 1.0f }, { 0, 0, 0 }, { 0, 0, 0 }, 0 },
		{ { 0.5f, 0.5f, 0.5f, 1.0f }, { 0, 0, 0 }, { 0, 0, 0 }, 0 } },
	};

	m_flavour.ApplyTo(&params);

	float guiscale[2];
	Gui::Screen::GetCoords2Pixels(guiscale);
	static float rot1, rot2;
	if (Pi::MouseButtonState(3)) {
		int m[2];
		Pi::GetMouseMotion(m);
		rot1 += -0.002*m[1];
		rot2 += -0.002*m[0];
	}
	else
	{
		rot1 += .5f*Pi::GetFrameTime();
		rot2 += Pi::GetFrameTime();
	}
	glClearColor(0.25,.37,.63,0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	
	const float bx = 5;
	const float by = 40;
	Gui::Screen::EnterOrtho();
	glColor3f(0,0,0);
	glBegin(GL_QUADS); {
		glVertex2f(bx,by);
		glVertex2f(bx,by+400);
		glVertex2f(bx+400,by+400);
		glVertex2f(bx+400,by);
	} glEnd();
	Gui::Screen::LeaveOrtho();
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glFrustum(-.5, .5, -.5, .5, 1.0f, 10000.0f);
	glDepthRange (0.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	Render::State::SetZnearZfar(1.0f, 10000.0f);

	float lightCol[] = { .5,.5,.5,0 };
	float lightDir[] = { 1,1,0,0 };

	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glLightfv(GL_LIGHT0, GL_POSITION, lightDir);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightCol);
	glLightfv(GL_LIGHT0, GL_AMBIENT, lightCol);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightCol);
	glEnable(GL_LIGHT0);
//	sbreSetDirLight (lightCol, lightDir);
	glViewport(
		GLint(bx/guiscale[0]),
		GLint((Gui::Screen::GetHeight() - by - 400)/guiscale[1]),
		GLsizei(400/guiscale[0]),
		GLsizei(400/guiscale[1]));
	
	matrix4x4f rot = matrix4x4f::RotateXMatrix(rot1);
	rot.RotateY(rot2);
	rot[14] = -1.5f * m_lmrModel->GetDrawClipRadius();

	m_lmrModel->Render(rot, &params);
	Render::State::UseProgram(0);
	Render::UnbindAllBuffers();
	glPopAttrib();
}

void StationViewShipView::ShowAll()
{
	const ShipType &t = ShipType::types[m_flavour.type];
	ReInit();

	SetTransparency(true);
	SpaceStation *station = Pi::player->GetDockedWith();
	SetTitle(stringf(256, "%s Shipyard", station->GetLabel().c_str()).c_str());
	
	Gui::Button *b = new Gui::SolidButton();
	b->onClick.connect(sigc::mem_fun(this, &StationViewShipView::Close));
	Add(b,680,470);
	Add(new Gui::Label("Go back"), 700, 470);
	
	b = new Gui::SolidButton();
	b->onClick.connect(sigc::mem_fun(this, &StationViewShipView::BuyShip));
	Add(b,480,470);
	Add(new Gui::Label("Buy this ship"), 500, 470);


	const float YSEP = floor(Gui::Screen::GetFontHeight() * 1.25f);
	float y = 40;
	Add(new Gui::Label("Ship type"), 420, y);
	Add(new Gui::Label(t.name), 600, y);
	y+=YSEP;
	Add(new Gui::Label("Price"), 420, y);
	Add(new Gui::Label(format_money(m_flavour.price)), 600, y);
	y+=YSEP;
	Add(new Gui::Label("Registration id"), 420, y);
	Add(new Gui::Label(m_flavour.regid), 600, y);
	y+=YSEP;
	y+=YSEP;
	Add(new Gui::Label("Weight empty"), 420, y);
	Add(new Gui::Label(stringf(64, "%d t", t.hullMass)), 600, y);
	y+=YSEP;
	Add(new Gui::Label("Weight fully loaded"), 420, y);
	Add(new Gui::Label(stringf(64, "%d t", t.hullMass + t.capacity)), 600, y);
	y+=YSEP;
	Add(new Gui::Label("Capacity"), 420, y);
	Add(new Gui::Label(stringf(64, "%d t", t.capacity)), 600, y);
	y+=YSEP;
	y+=YSEP;
	// forward accel
	float accel = t.linThrust[ShipType::THRUSTER_FORWARD] / (-9.81f*1000.0f*(t.hullMass));
	Add(new Gui::Label("Forward accel (empty)"), 420, y);
	Add(new Gui::Label(stringf(64, "%.1f G", accel)), 600, y);
	y+=YSEP;
	accel = t.linThrust[ShipType::THRUSTER_FORWARD] / (-9.81f*1000.0f*(t.hullMass + t.capacity));
	Add(new Gui::Label("Forward accel (laden)"), 420, y);
	Add(new Gui::Label(stringf(64, "%.1f G", accel)), 600, y);
	y+=YSEP;
	// rev accel
	accel = t.linThrust[ShipType::THRUSTER_REVERSE] / (9.81f*1000.0f*(t.hullMass));
	Add(new Gui::Label("Reverse accel (empty)"), 420, y);
	Add(new Gui::Label(stringf(64, "%.1f G", accel)), 600, y);
	y+=YSEP;
	accel = t.linThrust[ShipType::THRUSTER_REVERSE] / (9.81f*1000.0f*(t.hullMass + t.capacity));
	Add(new Gui::Label("Reverse accel (laden)"), 420, y);
	Add(new Gui::Label(stringf(64, "%.1f G", accel)), 600, y);
	y+=YSEP;
	y+=YSEP;
	Add(new Gui::Label("Hyperdrive fitted:"), 420, y);
	Add(new Gui::Label(EquipType::types[t.hyperdrive].name), 600, y);
	y+=YSEP;
	y+=YSEP;
	Add(new Gui::Label("Hyperspace range (fully laden):"), 420, y);
	y+=YSEP*1.5;

	{
		int row_size = 5, pos = 0;
		for (int x = 420, drivetype = Equip::DRIVE_CLASS1; drivetype <= Equip::DRIVE_CLASS9; x += 52, drivetype++) {
			if (t.capacity < EquipType::types[drivetype].mass)
				break;

			int hyperclass = EquipType::types[drivetype].pval;
			// for the sake of hyperspace range, we count ships mass as 60% of original.
			float range = Pi::CalcHyperspaceRange(hyperclass, t.hullMass + t.capacity);


			Add(new Gui::Label(stringf(128, "Class %d", hyperclass)), float(x), y);
			if (t.capacity < EquipType::types[drivetype].mass) {
				Add(new Gui::Label("---"), float(x), float(y+YSEP));
			} else {
				Add(new Gui::Label(stringf(128, "%.2f ly", range)), float(x), float(y+YSEP));
			}

			if (++pos == row_size) {
				pos = 0;
				x = 420-52;
				y += YSEP*2.5;
			}
		}
	}

	//AddBaseDisplay();
	//AddVideoWidget();

	Gui::Fixed::ShowAll();
}

////////////////////////////////////////////////////////////////////
#define RENDER_PICK_LINES
#ifdef RENDER_PICK_LINES
struct Ts2e {
	vector3f s_start;
	vector3f s_end;
};
static std::vector<Ts2e> s_pickLines;
#endif

////////////////////////////////////////////////////////////////////
class StationShipPaintView: public GenericChatForm {
public:
	StationShipPaintView(): GenericChatForm(), m_space(0), m_geom(0), m_cmesh(0), m_lmrModel(0), rot1(-0.6f), rot2(2.5f) {
		m_lmrModel = Pi::player->GetLmrModel();
		/* XXX duplicated code in InfoView.cpp */
		LmrObjParams params = {
			{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			{ /*"HELLO"*/ },
			{ 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 0.0f },

			{	// pColor[3]
			{ { 1.0f, 0.0f, 1.0f, 1.0f }, { 0, 0, 0 }, { 0, 0, 0 }, 0 },
			{ { 0.8f, 0.6f, 0.5f, 1.0f }, { 0, 0, 0 }, { 0, 0, 0 }, 0 },
			{ { 0.5f, 0.5f, 0.5f, 1.0f }, { 0, 0, 0 }, { 0, 0, 0 }, 0 } },
		};
		m_lmrParams = params;
		m_space = new CollisionSpace();
		Pi::player->GetFlavour()->ApplyTo(&m_lmrParams);
		m_ondraw3dcon = Pi::spaceStationView->onDraw3D.connect(sigc::mem_fun(this, &StationShipPaintView::Draw3D));
		m_colourIdx = 0;	// 1st colour

		RebuildCollMesh();
	}
	virtual ~StationShipPaintView() {
		m_ondraw3dcon.disconnect();

		// free our pointers
		delete m_cmesh;
		delete m_geom;
		delete m_space;
		
		// NOT owned by this, pointer is just a reference
		m_lmrModel=  NULL;
	
#ifdef RENDER_PICK_LINES
		s_pickLines.clear();
#endif
	}
	virtual void ShowAll();
	void Draw3D();
	virtual bool OnMouseDown(Gui::MouseButtonEvent *e);
private:
	void ChooseColour(int colourIdx) {
		m_colourIdx = colourIdx;
		if( m_colourIdx>=0 && m_colourIdx <3 ) {
			m_lmrParams.pMat[m_colourIdx].diffuse[0] = m_rAdjust.GetValue();
			m_lmrParams.pMat[m_colourIdx].diffuse[1] = m_gAdjust.GetValue();
			m_lmrParams.pMat[m_colourIdx].diffuse[2] = m_bAdjust.GetValue();
			m_lmrParams.pMat[m_colourIdx].diffuse[3] = 1.0f;
		}
	}
	void RebuildCollMesh();
	CollisionSpace* m_space;
	Geom* m_geom;
	LmrCollMesh* m_cmesh;
	LmrModel* m_lmrModel;
	LmrObjParams m_lmrParams;
	sigc::connection m_ondraw3dcon;
	uint8_t m_colourIdx;

	// Colour scrolling value things
	Gui::Adjustment m_rAdjust;
	Gui::Adjustment m_gAdjust;
	Gui::Adjustment m_bAdjust;

	// local copy of players PaintJob, all changes are made to this then copied
	// over the players _IF_ they like what they've done.
	CPaintJob m_paintJob;

	float rot1;
	float rot2;
};

void StationShipPaintView::RebuildCollMesh() 
{
	if( !m_space ) return;
	if( m_geom ) {
		m_space->RemoveGeom(m_geom);
		delete m_geom;
		delete m_cmesh;
	}

	m_cmesh = new LmrCollMesh(m_lmrModel, &m_lmrParams);
	m_geom = new Geom(m_cmesh->geomTree);
	m_space->AddGeom(m_geom);
	m_space->FlagRebuildObjectTrees();
	m_space->RebuildObjectTrees();
}

static void render_coll_mesh(const LmrCollMesh *m)
{
	glDisable(GL_LIGHTING);
	glColor3f(1,0,1);
	glDepthRange(0.0+(2.0/(1<<16)),1.0);
	glBegin(GL_TRIANGLES);
	for (int i=0; i<m->ni; i+=3) {
		glVertex3fv(&m->pVertex[3*m->pIndex[i]]);
		glVertex3fv(&m->pVertex[3*m->pIndex[i+1]]);
		glVertex3fv(&m->pVertex[3*m->pIndex[i+2]]);
	}
	glEnd();
	glColor3f(1,1,1);
	glDepthRange(0,1.0f-(2.0/(1<<16)));
	for (int i=0; i<m->ni; i+=3) {
		glBegin(GL_LINE_LOOP);
		glVertex3fv(&m->pVertex[3*m->pIndex[i]]);
		glVertex3fv(&m->pVertex[3*m->pIndex[i+1]]);
		glVertex3fv(&m->pVertex[3*m->pIndex[i+2]]);
		glEnd();
	}
	glDepthRange(0,1);
	glEnable(GL_LIGHTING);
}

static vector3f GetOGLPos(const int x, const int y)
{
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	GLfloat winX, winY, winZ;
	GLdouble posX, posY, posZ;

	glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
	glGetDoublev( GL_PROJECTION_MATRIX, projection );
	glGetIntegerv( GL_VIEWPORT, viewport );

	winX = (float)x;
	winY = (float)y;
	winZ = 1.0f;
	//glReadPixels( x, y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );

	gluUnProject( winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);

	return vector3f(posX, posY, posZ);
}

static void calcRay(const int x, const int y, vector3f &p1, vector3f &p2)
{
	float guiscale[2];
	Gui::Screen::GetCoords2Pixels(guiscale);
	matrix4x4f invMatrix;
	matrix4x4f viewMatrix;
	matrix4x4f projMatrix;
	const float	fNEAR	= 1.0f;
	const float	fFAR	= 10000.0f;
	//const float	fFOV = 0.8f;
	//glFrustum(-.5, .5, -.5, .5, 1.0f, 10000.0f);
	const float tminusb	= (.5)-(-.5);
	//const float	fFOV	= atan((tminusb*0.5) / fNEAR) * 2.0;
	const float	fFOV	= atan((tminusb*0.5) / fNEAR);
	const float	fWIDTH	= 400.0f;
	const float	fHEIGHT	= 400.0f;
	const float	fWIDTH_DIV_2 = (fWIDTH*0.5f);
	const float	fHEIGHT_DIV_2 = (fHEIGHT*0.5f);
	const float	fASPECT = 1.0f;

	const float dx=tanf(fFOV*0.5f)*(x/fWIDTH_DIV_2-1.0f)/fASPECT;
	const float dy=tanf(fFOV*0.5f)*(1.0f-y/fHEIGHT_DIV_2);

	glGetFloatv( GL_MODELVIEW_MATRIX, &viewMatrix[0] );
	glGetFloatv( GL_PROJECTION_MATRIX, &projMatrix[0] );
	invMatrix = (viewMatrix * projMatrix).InverseOf();

	p1=vector3f(dx*fNEAR, dy*fNEAR, fNEAR);
	p2=vector3f(dx*fFAR, dy*fFAR, fFAR);

	p1 = p1 * invMatrix;
	p2 = p2 * invMatrix;
}

static vector3f calcDir(const int x, const int y)
{
	vector3f p1,p2;
	calcRay(x,y,p1,p2);
	return (p2-p1).Normalized();
}

#ifdef RENDER_PICK_LINES
static void render_pick_line()
{
	glDisable(GL_LIGHTING);
	glColor3f(1,1,1);
	glBegin(GL_LINES);
	for( std::vector<Ts2e>::const_iterator iter = s_pickLines.begin() ; iter != s_pickLines.end() ; ++iter ) {
		glVertex3fv(&iter->s_start[0]);
		glVertex3fv(&iter->s_end[0]);
	}
	glEnd();
	glEnable(GL_LIGHTING);
}
#endif

// override the OnMouseDown virutal method so I can capture the mouse position
bool StationShipPaintView::OnMouseDown(Gui::MouseButtonEvent *e)
{
	if( e->button != 1 ) return false;
	const float bx = 5.0f;
	const float by = 40.0f;
	// work out the screen dimensions
	const float left = bx;
	const float top = by;
	const float bottom = by+400.0f;
	const float right = bx+400.0f;

	const float sx = e->screenX;
	const float sy = e->screenY;

	if(sx >= left && sx <= right) {
		if(sy >= top && sy <= bottom) {

			vector3f res1;
			glPushAttrib(GL_ALL_ATTRIB_BITS); {
				glPushMatrix(); {
					float guiscale[2];
					Gui::Screen::GetCoords2Pixels(guiscale);
					const float bx = 5;
					const float by = 40;
	
					glMatrixMode(GL_PROJECTION);
					glLoadIdentity();
					glFrustum(-.5, .5, -.5, .5, 1.0f, 10000.0f);
					glDepthRange (0.0, 1.0);

					glPushMatrix(); {
						glMatrixMode(GL_MODELVIEW);
						glLoadIdentity();
	
						glViewport(
							GLint(bx/guiscale[0]),
							GLint((Gui::Screen::GetHeight() - by - 400)/guiscale[1]),
							GLsizei(400/guiscale[0]),
							GLsizei(400/guiscale[1]));

						//vector3f p1, p2;
						const int sl = sx-left;
						const int st = sy-top;
						res1 = calcDir( sl, st );

					} glPopMatrix();
				} glPopMatrix();
			} glPopAttrib();
			

			// Origin of the ray is just the camera position
			vector3f camPos( 0.0, 0.0, 2.0f * m_lmrModel->GetDrawClipRadius() );
			//vector3f camPos( 0.0, 0.0, 0.0 );

			// inverted rotation the ship goes through
			matrix4x4f rot = matrix4x4f::RotateXMatrix(rot1);
			rot.RotateY(rot2);
			//rot[14] = -2.0f * m_lmrModel->GetDrawClipRadius();
			const matrix4x4f invrot = rot.InverseOf();
			res1 = invrot * res1;
			camPos = invrot * camPos;

#ifdef RENDER_PICK_LINES
			Ts2e t;
			t.s_start = camPos;
			t.s_end = camPos + (res1 * 10000.0f);
			s_pickLines.push_back(t);
#endif

			// trace laser beam through frame to see who it hits
			CollisionContact c;
			m_geom->Enable();
			m_space->TraceRay(camPos, res1, 10000.0, &c);
			m_geom->Disable();
			if (c.triIdx != -1) {
				vector3f tangent = c.normal;
				m_paintJob.AddDecal(c.pos, c.normal, tangent, 1.0f, 1.0f, 1.0f);
			}
		}
	}

	return Gui::Container::OnMouseDown(e);
}

void StationShipPaintView::Draw3D()
{
	float guiscale[2];
	Gui::Screen::GetCoords2Pixels(guiscale);
	
	if (Pi::MouseButtonState(3)) {
		int m[2];
		Pi::GetMouseMotion(m);
		rot1 += -0.002*m[1];
		rot2 += -0.002*m[0];
	}
	glClearColor(0.25,.37,.63,0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	
	const float bx = 5;
	const float by = 40;
	Gui::Screen::EnterOrtho();
	glColor3f(0,0,0);
	glBegin(GL_QUADS); {
		glVertex2f(bx,by);
		glVertex2f(bx,by+400);
		glVertex2f(bx+400,by+400);
		glVertex2f(bx+400,by);
	} glEnd();
	Gui::Screen::LeaveOrtho();
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glFrustum(-.5, .5, -.5, .5, 1.0f, 10000.0f);
	glDepthRange (0.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	Render::State::SetZnearZfar(1.0f, 10000.0f);

	float lightCol[] = { .5,.5,.5,0 };
	float lightDir[] = { 1,1,0,0 };

	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glLightfv(GL_LIGHT0, GL_POSITION, lightDir);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightCol);
	glLightfv(GL_LIGHT0, GL_AMBIENT, lightCol);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightCol);
	glEnable(GL_LIGHT0);
//	sbreSetDirLight (lightCol, lightDir);
	glViewport(
		GLint(bx/guiscale[0]),
		GLint((Gui::Screen::GetHeight() - by - 400)/guiscale[1]),
		GLsizei(400/guiscale[0]),
		GLsizei(400/guiscale[1]));
	
	matrix4x4f rot = matrix4x4f::RotateXMatrix(rot1);
	rot.RotateY(rot2);
	rot[14] = -2.0f * m_lmrModel->GetDrawClipRadius();

#ifdef RENDER_PICK_LINES
	m_lmrModel->Render(rot, &m_lmrParams);
	glPushMatrix();
		glMultMatrixf(&rot[0]);
		/*render_coll_mesh(m_cmesh);*/
		render_pick_line();
	glPopMatrix();
#else
	m_lmrModel->Render(rot, &m_lmrParams);
#endif

	Render::State::UseProgram(0);
	Render::UnbindAllBuffers();
	glPopAttrib();
}

void StationShipPaintView::ShowAll()
{
	ReInit();

	SetTransparency(true);
	SpaceStation *station = Pi::player->GetDockedWith();
	SetTitle(stringf(256, "%s Shipyard", station->GetLabel().c_str()).c_str());
	
	Gui::Button *b = new Gui::SolidButton();
	b->onClick.connect(sigc::mem_fun(this, &StationShipPaintView::Close));
	Add(b,680,470);
	Add(new Gui::Label("Go back"), 700, 470);
	
	// ----------------------------
	Gui::Button *i = new Gui::SolidButton();
	i->onClick.connect(sigc::bind(sigc::mem_fun(this, &StationShipPaintView::ChooseColour), 0));
	Add(i,480,70);
	Add(new Gui::Label("1st"), 500, 70);

	i = new Gui::SolidButton();
	i->onClick.connect(sigc::bind(sigc::mem_fun(this, &StationShipPaintView::ChooseColour), 1));
	Add(i,540,70);
	Add(new Gui::Label("2nd"), 560, 70);

	i = new Gui::SolidButton();
	i->onClick.connect(sigc::bind(sigc::mem_fun(this, &StationShipPaintView::ChooseColour), 2));
	Add(i,600,70);
	Add(new Gui::Label("3rd"), 620, 70);

	// ----------------------------
	Gui::HScrollBar *scroll = new Gui::HScrollBar();
	scroll->SetAdjustment(&m_rAdjust);
	Add(scroll, 455, 160);

	scroll = new Gui::HScrollBar();
	scroll->SetAdjustment(&m_gAdjust);
	Add(scroll, 455, 200);

	scroll = new Gui::HScrollBar();
	scroll->SetAdjustment(&m_bAdjust);
	Add(scroll, 455, 240);

	Gui::Fixed::ShowAll();
}

////////////////////////////////////////////////////////////////////

class StationShipRepairsView: public GenericChatForm {
public:
	StationShipRepairsView();
	virtual ~StationShipRepairsView() {
	}
	virtual void ShowAll();
private:
	int GetCostOfFixingHull(float percent) {
		return int(Pi::player->GetFlavour()->price * 0.001 * percent);
	}

	void RepairHull(float percent) {
		int cost = GetCostOfFixingHull(percent);
		if (Pi::player->GetMoney() < cost) {
			Pi::cpan->MsgLog()->Message("", "You do not have enough money");
		} else {
			Pi::player->SetMoney(Pi::player->GetMoney() - cost);
			Pi::player->SetPercentHull(Pi::player->GetPercentHull() + percent);
			ShowAll();
		}
	}
//	sigc::connection m_onShipsForSaleChangedConnection;
};

StationShipRepairsView::StationShipRepairsView(): GenericChatForm()
{
}

void StationShipRepairsView::ShowAll()
{
	ReInit();

	SpaceStation *station = Pi::player->GetDockedWith();
	assert(station);
	SetTransparency(false);
	
	SetTitle(stringf(256, "%s Shipyard", station->GetLabel().c_str()).c_str());
	
	Gui::Button *b = new Gui::SolidButton();
	b->onClick.connect(sigc::mem_fun(this, &StationShipRepairsView::Close));
	Add(b,680,470);
	Add(new Gui::Label("Go back"), 700, 470);

	Gui::Fixed *fbox = new Gui::Fixed(470, 400);
	Add(fbox, 320, 40);

	const float YSEP = floor(Gui::Screen::GetFontHeight() * 1.5f);
	float ypos = YSEP;

	float hullPercent = Pi::player->GetPercentHull();
	if (hullPercent >= 100.0f) {
		fbox->Add(new Gui::Label("Your ship is in perfect working condition."), 0, ypos);
	} else {
		int costAll = GetCostOfFixingHull(100.0f - hullPercent);
		int cost1 = GetCostOfFixingHull(1.0f);
		if (cost1 < costAll) {
			fbox->Add(new Gui::Label("Repair 1.0% of hull damage"), 0, ypos);
			fbox->Add(new Gui::Label(format_money(cost1)), 350, ypos);
			
			Gui::SolidButton *sb = new Gui::SolidButton();
			sb->onClick.connect(sigc::bind(sigc::mem_fun(this, &StationShipRepairsView::RepairHull), 1.0f));
			fbox->Add(sb, 430, ypos);
			ypos += YSEP;
		}
		fbox->Add(new Gui::Label(stringf(128, "Repair all hull damage (%.1f%%)", 100.0f-hullPercent)), 0, ypos);
		fbox->Add(new Gui::Label(format_money(costAll)), 350, ypos);
		
		Gui::SolidButton *sb = new Gui::SolidButton();
		sb->onClick.connect(sigc::bind(sigc::mem_fun(this, &StationShipRepairsView::RepairHull), 100.0f-hullPercent));
		fbox->Add(sb, 430, ypos);
	}

	const float *col = Gui::Theme::Colors::tableHeading;
	fbox->Add((new Gui::Label("Item"))->Color(col), 0, 0);
	fbox->Add((new Gui::Label("Price"))->Color(col), 350, 0);
	fbox->Add((new Gui::Label("Repair"))->Color(col), 430, 0);
	fbox->ShowAll();
	AddBaseDisplay();
	AddVideoWidget();

	Gui::Fixed::ShowAll();
}


////////////////////////////////////////////////////////////////////

class StationBuyShipsView: public GenericChatForm {
public:
	StationBuyShipsView();
	virtual ~StationBuyShipsView() {
		m_onShipsForSaleChangedConnection.disconnect();
	}
	virtual void ShowAll();
private:
	void ViewShip(int idx) {
		OpenChildChatForm(new StationViewShipView(idx));
	}
	sigc::connection m_onShipsForSaleChangedConnection;
};

StationBuyShipsView::StationBuyShipsView(): GenericChatForm()
{
	SpaceStation *station = Pi::player->GetDockedWith();
	m_onShipsForSaleChangedConnection = station->onShipsForSaleChanged.connect(
			sigc::mem_fun(this, &StationBuyShipsView::ShowAll));
	SetTransparency(false);
}

void StationBuyShipsView::ShowAll()
{
	ReInit();

	SpaceStation *station = Pi::player->GetDockedWith();
	assert(station);
	SetTransparency(false);
	
	SetTitle(stringf(256, "%s Shipyard", station->GetLabel().c_str()).c_str());
	
	Gui::Button *b = new Gui::SolidButton();
	b->onClick.connect(sigc::mem_fun(this, &StationBuyShipsView::Close));
	Add(b,680,470);
	Add(new Gui::Label("Go back"), 700, 470);

	Gui::Fixed *fbox = new Gui::Fixed(470, 400);
	Add(fbox, 320, 40);

	Gui::VScrollBar *scroll = new Gui::VScrollBar();
	Gui::VScrollPortal *portal = new Gui::VScrollPortal(450,400);
	scroll->SetAdjustment(&portal->vscrollAdjust);

	std::vector<ShipFlavour> &ships = station->GetShipsOnSale();
	int NUM_ITEMS = ships.size();
	const float YSEP = floor(Gui::Screen::GetFontHeight() * 1.5f);

	int num = 0;
	Gui::Fixed *innerbox = new Gui::Fixed(450, NUM_ITEMS*YSEP);
	for (std::vector<ShipFlavour>::iterator i = ships.begin(); i!=ships.end(); ++i) {
		Gui::Label *l = new Gui::Label(ShipType::types[(*i).type].name);
		innerbox->Add(l,0,num*YSEP);
		innerbox->Add(new Gui::Label(format_money((*i).price)), 200, num*YSEP);
		innerbox->Add(new Gui::Label(format_money((*i).price - Pi::player->GetFlavour()->price) ), 275, num*YSEP);
		innerbox->Add(new Gui::Label(stringf(16, "%dt", ShipType::types[(*i).type].capacity)), 370, num*YSEP);
		
		Gui::SolidButton *sb = new Gui::SolidButton();
		sb->onClick.connect(sigc::bind(sigc::mem_fun(this, &StationBuyShipsView::ViewShip), num));
		innerbox->Add(sb, 430, num*YSEP);
		num++;
	}
	innerbox->ShowAll();

	const float *col = Gui::Theme::Colors::tableHeading;
	fbox->Add((new Gui::Label("Ship"))->Color(col), 0, 0);
	fbox->Add((new Gui::Label("Price"))->Color(col), 200, 0);
	fbox->Add((new Gui::Label("Part exchange"))->Color(col), 275, 0);
	fbox->Add((new Gui::Label("Capacity"))->Color(col), 370, 0);
	fbox->Add((new Gui::Label("View"))->Color(col), 430, 0);
	fbox->Add(portal, 0, YSEP);
	fbox->Add(scroll, 455, YSEP);
	portal->Add(innerbox);
	portal->ShowAll();
	fbox->ShowAll();
	AddBaseDisplay();
	AddVideoWidget();

	Gui::Fixed::ShowAll();
}


////////////////////////////////////////////////////////////////////

class StationShipyardView: public GenericChatForm {
public:
	StationShipyardView();
	virtual void ShowAll();
private:
	void GotoUpgradesView() {
		OpenChildChatForm(new StationShipUpgradesView());
	}
	void GotoBuyShipsView() {
		OpenChildChatForm(new StationBuyShipsView());
	}
	void GotoShipRepairsView() {
		OpenChildChatForm(new StationShipRepairsView());
	}
	void GotoShipPainterView() {
		OpenChildChatForm(new StationShipPaintView());
	}
};

StationShipyardView::StationShipyardView(): GenericChatForm()
{
	SetTransparency(false);
}

void StationShipyardView::ShowAll()
{
	ReInit();

	SpaceStation *station = Pi::player->GetDockedWith();
	assert(station);
	SetTransparency(false);
	SetTitle(stringf(256, "%s Shipyard", station->GetLabel().c_str()).c_str());
	
	Gui::Button *backButton = new Gui::SolidButton();
	backButton->onClick.connect(sigc::mem_fun(this, &StationShipyardView::Close));
	Add(backButton,680,470);
	Add(new Gui::Label("Go back"), 700, 470);

	Gui::SolidButton *b = new Gui::SolidButton();
	b->SetShortcut(SDLK_1, KMOD_NONE);
	b->onClick.connect(sigc::mem_fun(this, &StationShipyardView::GotoUpgradesView));
	Add(b, 340, 240);
	Gui::Label *l = new Gui::Label("Ship equipment");
	Add(l, 365, 240);
	
	b = new Gui::SolidButton();
	b->SetShortcut(SDLK_2, KMOD_NONE);
	b->onClick.connect(sigc::mem_fun(this, &StationShipyardView::GotoShipRepairsView));
	Add(b, 340, 300);
	l = new Gui::Label("Repairs and servicing");
	Add(l, 365, 300);
	
	b = new Gui::SolidButton();
	b->SetShortcut(SDLK_3, KMOD_NONE);
	b->onClick.connect(sigc::mem_fun(this, &StationShipyardView::GotoBuyShipsView));
	Add(b, 340, 360);
	l = new Gui::Label("New and reconditioned ships");
	Add(l, 365, 360);

	b = new Gui::SolidButton();
	b->SetShortcut(SDLK_3, KMOD_NONE);
	b->onClick.connect(sigc::mem_fun(this, &StationShipyardView::GotoShipPainterView));
	Add(b, 340, 420);
	l = new Gui::Label("Paint-ship Pro!");	// part of me is sorry for this name... a small and totally ignored part :D
	Add(l, 365, 420);

	AddBaseDisplay();
	AddVideoWidget();

	Gui::Fixed::ShowAll();
}

////////////////////////////////////////////////////////////////////

class StationBBView: public GenericChatForm {
public:
	StationBBView();
	virtual ~StationBBView() {
		m_onBBChangedConnection.disconnect();
		m_onBBAdvertDeleted.disconnect();
	}
	virtual void ShowAll();
private:
	void OpenAdvertDialog(int ref);
	void OnCloseAdvertDialog(GenericChatForm *advertDlg);
	void OnAdvertDeleted(BBAdvert *);
	void OnBBChanged();
	sigc::connection m_onBBChangedConnection, m_onBBAdvertDeleted;
	BBAdvertChatForm *m_advertChatForm;
};

StationBBView::StationBBView(): GenericChatForm()
{
	m_advertChatForm = 0;
	SpaceStation *station = Pi::player->GetDockedWith();
	m_onBBChangedConnection = station->onBulletinBoardChanged.connect(
			sigc::mem_fun(this, &StationBBView::OnBBChanged));
	m_onBBAdvertDeleted = station->onBulletinBoardAdvertDeleted.connect(
			sigc::mem_fun(this, &StationBBView::OnAdvertDeleted));
	SetTransparency(false);
}

void StationBBView::OnBBChanged()
{
	if (!m_advertChatForm) {
		ShowAll();
	}
}

void StationBBView::OnAdvertDeleted(BBAdvert *ad)
{
	if (m_advertChatForm) {
	       if (ad == m_advertChatForm->GetAdvert()) m_advertChatForm->OnAdvertDeleted();
	} else {
		ShowAll();
	}
}

void StationBBView::OnCloseAdvertDialog(GenericChatForm *advertDlg)
{
	m_advertChatForm = 0;
}

void StationBBView::OpenAdvertDialog(int ref)
{
	if (m_advertChatForm) m_advertChatForm->Close(); // shouldn't happen...

	SpaceStation *station = Pi::player->GetDockedWith();
	const BBAdvert *ad = station->GetBBAdvert(ref);
	
	m_advertChatForm = ad->builder(station, ad);

	m_advertChatForm->onClose.connect(sigc::mem_fun(this, &StationBBView::OnCloseAdvertDialog));
	m_advertChatForm->AddBaseDisplay();
	m_advertChatForm->AddVideoWidget();

    m_advertChatForm->CallDialogHandler(0);

	OpenChildChatForm(m_advertChatForm);
}

void StationBBView::ShowAll()
{
	ReInit();

	SpaceStation *station = Pi::player->GetDockedWith();
	assert(station);
	SetTransparency(false);
	
	SetTitle((station->GetLabel() + " Bulletin Board").c_str());
	
	Gui::Button *b = new Gui::SolidButton();
	b->onClick.connect(sigc::mem_fun(this, &StationBBView::Close));
	Add(b,680,470);
	Add(new Gui::Label("Go back"), 700, 470);

	Gui::Fixed *fbox = new Gui::Fixed(470, 400);
	Add(fbox, 320, 40);

	Gui::VScrollBar *scroll = new Gui::VScrollBar();
	Gui::VScrollPortal *portal = new Gui::VScrollPortal(450,400);
	scroll->SetAdjustment(&portal->vscrollAdjust);

	const std::list<const BBAdvert*> bbadverts = station->GetBBAdverts();

	int NUM_ITEMS = bbadverts.size();
	const float YSEP = floor(Gui::Screen::GetFontHeight() * 5);

	int num = 0;
	Gui::Fixed *innerbox = new Gui::Fixed(450, NUM_ITEMS*YSEP);
	for (std::list<const BBAdvert*>::const_iterator i = bbadverts.begin(); i != bbadverts.end(); i++) {
		Gui::SolidButton *sb = new Gui::SolidButton();
		sb->onClick.connect(sigc::bind(sigc::mem_fun(this, &StationBBView::OpenAdvertDialog), (*i)->ref));
		innerbox->Add(sb, 10, num*YSEP);
		
		Gui::Label *l = new Gui::Label((*i)->description);
		innerbox->Add(l,40,num*YSEP);
		
		num++;
	}
	innerbox->ShowAll();

	fbox->Add(portal, 0, 10);
	fbox->Add(scroll, 455, 10);
	portal->Add(innerbox);
	portal->ShowAll();
	fbox->ShowAll();
	AddBaseDisplay();
	AddVideoWidget();

	Gui::Fixed::ShowAll();
}

/////////////////////////////////////////////////////////////////////

class StationRootView: public GenericChatForm {
public:
	StationRootView() {}
	virtual ~StationRootView() {}
	virtual void ShowAll();
private:
	void GotoShipyard();
	void GotoCommodities();
	void GotoBB();
	void GotoPolis();
	void OnClickRequestLaunch();
};

void StationRootView::ShowAll()
{

	SetTransparency(false);
	ReInit();
	SpaceStation *station = Pi::player->GetDockedWith();
	assert(station);
	{
		char buf[256];
		snprintf(buf, sizeof(buf), "Welcome to %s", station->GetLabel().c_str());
		SetTitle(buf);
	}
	Gui::Label *l = new Gui::Label("Hello friend! Thankyou for docking with this space station!\n"
	"We regret to inform you that due to a spacetime fissure you have "
	"ended up in the terrible mirror universe of boringness, and that there "
	"is nothing to do. We hope to have things back to normal within a "
	"few weeks, and in the mean time would like to offer our apologies for "
	"any loss of earnings, immersion or lunch.  "
	"Currently the usual space station services are not available, but we "
	"can offer you this promotional message from one of the station's sponsors:\n"
	"                       DIET STEAKETTE: IT'S BAD");

	Gui::Fixed *fbox = new Gui::Fixed(450, 400);
	fbox->Add(l, 0, 0);
	Add(fbox, 320, 40);
	fbox->ShowAll();

	Gui::SolidButton *b = new Gui::SolidButton();
	b->SetShortcut(SDLK_1, KMOD_NONE);
	b->onClick.connect(sigc::mem_fun(this, &StationRootView::OnClickRequestLaunch));
	Add(b, 340, 220);
	l = new Gui::Label("Request Launch");
	Add(l, 365, 220);

	b = new Gui::SolidButton();
	b->SetShortcut(SDLK_2, KMOD_NONE);
	b->onClick.connect(sigc::mem_fun(this, &StationRootView::GotoShipyard));
	Add(b, 340, 280);
	l = new Gui::Label("Shipyard");
	Add(l, 365, 280);

	b = new Gui::SolidButton();
	b->SetShortcut(SDLK_3, KMOD_NONE);
	b->onClick.connect(sigc::mem_fun(this, &StationRootView::GotoCommodities));
	Add(b, 340, 340);
	l = new Gui::Label("Commodity market");
	Add(l, 365, 340);

	b = new Gui::SolidButton();
	b->SetShortcut(SDLK_4, KMOD_NONE);
	b->onClick.connect(sigc::mem_fun(this, &StationRootView::GotoBB));
	Add(b, 340, 400);
	l = new Gui::Label("Bulletin board");
	Add(l, 365, 400);

	b = new Gui::SolidButton();
	b->SetShortcut(SDLK_5, KMOD_NONE);
	b->onClick.connect(sigc::mem_fun(this, &StationRootView::GotoPolis));
	Add(b, 340, 460);
	l = new Gui::Label("Contact local police");
	Add(l, 365, 460);

	AddBaseDisplay();
	AddVideoWidget();

	Gui::Fixed::ShowAll();
}

void StationRootView::OnClickRequestLaunch()
{
	Pi::worldView->OnClickBlastoff();
	Pi::SetView(Pi::worldView);
}

void StationRootView::GotoCommodities()
{
	OpenChildChatForm(new StationCommoditiesView());
}

void StationRootView::GotoShipyard()
{
	OpenChildChatForm(new StationShipyardView());
}

void StationRootView::GotoBB()
{
	OpenChildChatForm(new StationBBView());
}

void StationRootView::GotoPolis()
{
	OpenChildChatForm(new PoliceChatForm());
}

/////////////////////////////////////////////////////////////////////

SpaceStationView::SpaceStationView(): View(), m_jumpToForm(0)
{
	Gui::Label *l = new Gui::Label("Comms Link");
	l->Color(1,.7,0);
	m_rightRegion2->Add(l, 10, 0);
}

void SpaceStationView::JumpToForm(GenericChatForm *form)
{
	m_jumpToForm = form;
}

void SpaceStationView::OnSwitchTo()
{
	DeleteAllChildren();
	SetTransparency(true);
	m_baseSubView = new StationRootView();
	Add(m_baseSubView, 0, 0);
	m_baseSubView->ShowAll();
}

void SpaceStationView::Draw3D()
{
	onDraw3D.emit();
}

void SpaceStationView::Update()
{
	if (m_jumpToForm) {
		OnSwitchTo();
		m_baseSubView->OpenChildChatForm(m_jumpToForm);
		m_jumpToForm = 0;
	}
}
