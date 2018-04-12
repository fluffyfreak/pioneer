// Copyright © 2008-2018 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "HyperspaceCloud.h"
#include "libs.h"
#include "Game.h"
#include "Lang.h"
#include "perlin.h"
#include "Pi.h"
#include "Player.h"
#include "Ship.h"
#include "Space.h"
#include "graphics/Graphics.h"
#include "graphics/Material.h"
#include "graphics/Renderer.h"
#include "graphics/VertexArray.h"
#include "graphics/RenderState.h"
#include "json/JsonUtils.h"
#include "GameSaveError.h"

using namespace Graphics;

/** How long does a hyperspace cloud last for? 2 Days? */
#define HYPERCLOUD_DURATION (60.0*60.0*24.0*2.0)

HyperspaceCloud::HyperspaceCloud(Ship *s, double dueDate, bool isArrival) : m_isBeingKilled(false)
{
	m_flags = Body::FLAG_CAN_MOVE_FRAME |
		  Body::FLAG_LABEL_HIDDEN;
	m_ship = s;
	SetPhysRadius(0.0);
	SetClipRadius(1200.0);
	m_vel = (s ? s->GetVelocity() : vector3d(0.0));
	m_birthdate = Pi::game->GetTime();
	m_due = dueDate;
	SetIsArrival(isArrival);
	InitGraphics();
}

HyperspaceCloud::HyperspaceCloud() : m_isBeingKilled(false)
{
	m_ship = 0;
	SetPhysRadius(0.0);
	SetClipRadius(1200.0);
	InitGraphics();
}

void HyperspaceCloud::InitGraphics()
{
	Graphics::MaterialDescriptor desc;
	desc.vertexColors = true;
	material.Reset(Pi::renderer->CreateMaterial(desc));

	Graphics::RenderStateDesc rsd;
	rsd.blendMode  = BLEND_ALPHA_ONE;
	rsd.depthWrite = false;
	renderState = Pi::renderer->CreateRenderState(rsd);

	m_graphic.reset(new TriangleGraphic(Pi::renderer, renderState, new Graphics::VertexArray(ATTRIB_POSITION | ATTRIB_DIFFUSE), material, Graphics::TRIANGLE_FAN));
}

HyperspaceCloud::~HyperspaceCloud()
{
	if (m_ship) delete m_ship;
}

void HyperspaceCloud::SetIsArrival(bool isArrival)
{
	m_isArrival = isArrival;
	SetLabel(isArrival ? Lang::HYPERSPACE_ARRIVAL_CLOUD : Lang::HYPERSPACE_DEPARTURE_CLOUD);
}

void HyperspaceCloud::SaveToJson(Json::Value &jsonObj, Space *space)
{
	Body::SaveToJson(jsonObj, space);

	Json::Value hyperspaceCloudObj(Json::objectValue); // Create JSON object to contain hyperspace cloud data.

	VectorToJson(hyperspaceCloudObj, m_vel, "vel");
	hyperspaceCloudObj["birth_date"] = DoubleToStr(m_birthdate);
	hyperspaceCloudObj["due"] = DoubleToStr(m_due);
	hyperspaceCloudObj["is_arrival"] = m_isArrival;
	if (m_ship)
	{
		Json::Value shipObj(Json::objectValue); // Create JSON object to contain ship data.
		m_ship->ToJson(shipObj, space);
		hyperspaceCloudObj["ship"] = shipObj; // Add ship object to hyperpace cloud object.
	}

	jsonObj["hyperspace_cloud"] = hyperspaceCloudObj; // Add hyperspace cloud object to supplied object.
}

void HyperspaceCloud::LoadFromJson(const Json::Value &jsonObj, Space *space)
{
	Body::LoadFromJson(jsonObj, space);

	if (!jsonObj.isMember("hyperspace_cloud")) throw SavedGameCorruptException();
	Json::Value hyperspaceCloudObj = jsonObj["hyperspace_cloud"];

	if (!hyperspaceCloudObj.isMember("vel")) throw SavedGameCorruptException();
	if (!hyperspaceCloudObj.isMember("birth_date")) throw SavedGameCorruptException();
	if (!hyperspaceCloudObj.isMember("due")) throw SavedGameCorruptException();
	if (!hyperspaceCloudObj.isMember("is_arrival")) throw SavedGameCorruptException();

	JsonToVector(&m_vel, hyperspaceCloudObj, "vel");
	m_birthdate = StrToDouble(hyperspaceCloudObj["birth_date"].asString());
	m_due = StrToDouble(hyperspaceCloudObj["due"].asString());
	m_isArrival = hyperspaceCloudObj["is_arrival"].asBool();

	if (hyperspaceCloudObj.isMember("ship"))
	{
		Json::Value shipObj = hyperspaceCloudObj["ship"];
		m_ship = static_cast<Ship*>(Body::FromJson(shipObj, space));
	}
}

void HyperspaceCloud::PostLoadFixup(Space *space)
{
	Body::PostLoadFixup(space);
	if (m_ship) m_ship->PostLoadFixup(space);
}

void HyperspaceCloud::TimeStepUpdate(const float timeStep)
{
	if( m_isBeingKilled )
		return;

	SetPosition(GetPosition() + m_vel * timeStep);

	if (m_isArrival && m_ship && (m_due < Pi::game->GetTime())) {
		// spawn ship
		// XXX some overlap with Space::DoHyperspaceTo(). should probably all
		// be moved into EvictShip()
		m_ship->SetPosition(GetPosition());
		m_ship->SetVelocity(m_vel);
		m_ship->SetOrient(matrix3x3d::Identity());
		m_ship->SetFrame(GetFrame());
		Pi::game->GetSpace()->AddBody(m_ship);

		if (Pi::player->GetNavTarget() == this && !Pi::player->GetCombatTarget())
			Pi::player->SetCombatTarget(m_ship, Pi::player->GetSetSpeedTarget() == this);

		m_ship->EnterSystem();

		m_ship = 0;
	}

	// cloud expiration
	if( m_birthdate + HYPERCLOUD_DURATION <= Pi::game->GetTime() )
	{
		Pi::game->RemoveHyperspaceCloud(this);
		Pi::game->GetSpace()->KillBody(this);
		m_isBeingKilled = true;
	}
}

Ship *HyperspaceCloud::EvictShip()
{
	Ship *s = m_ship;
	m_ship = 0;
	return s;
}

static void make_circle_thing(VertexArray *va, float radius, const Color &colCenter, const Color &colEdge)
{
	va->Clear();
	va->Add(vector3f(0.f, 0.f, 0.f), colCenter);
	for (float ang=0; ang<float(M_PI)*2.f; ang+=0.1f) {
		va->Add(vector3f(radius*sin(ang), radius*cos(ang), 0.0f), colEdge);
	}
	va->Add(vector3f(0.f, radius, 0.f), colEdge);
}

void HyperspaceCloud::UpdateInterpTransform(double alpha)
{
	m_interpOrient = matrix3x3d::Identity();
	const vector3d oldPos = GetPosition() - m_vel*Pi::game->GetTimeStep();
	m_interpPos = alpha*GetPosition() + (1.0-alpha)*oldPos;
}

void HyperspaceCloud::Render(Renderer *renderer, Camera *camera, const vector3d &viewCoords, const matrix4x4d &viewTransform)
{
	if( m_isBeingKilled )
		return;

	matrix4x4f trans = matrix4x4f::Identity();
	trans.Translate(float(viewCoords.x), float(viewCoords.y), float(viewCoords.z));

	// face the camera dammit
	vector3f zaxis = vector3f(viewCoords.NormalizedSafe());
	vector3f xaxis = vector3f(0,1,0).Cross(zaxis).Normalized();
	vector3f yaxis = zaxis.Cross(xaxis);
	matrix4x4f rot = matrix4x4f::MakeRotMatrix(xaxis, yaxis, zaxis).Inverse();
	//renderer->SetTransform(trans * rot);

	// precise to the rendered frame (better than PHYSICS_HZ granularity)
	const double preciseTime = Pi::game->GetTime() + Pi::GetGameTickAlpha()*Pi::game->GetTimeStep();

	// Flickering gradient circle, departure clouds are red and arrival clouds blue
	// XXX could just alter the scale instead of recreating the model
	const float radius = 1000.0f + 200.0f*float(noise(vector3d(10.0*preciseTime, 0, 0)));
	
	Color outerColor = m_isArrival ? Color::BLUE : Color::RED;
	outerColor.a = 0;
	make_circle_thing(m_graphic->GetVertexArray(), radius, Color::WHITE, outerColor);
	//renderer->DrawTriangles(vertices.get(), renderState, material.get(), TRIANGLE_FAN);

	m_graphic->SetTransform(trans * rot);
	camera->GetCollector().AddTransparent(m_graphic.get());
}
