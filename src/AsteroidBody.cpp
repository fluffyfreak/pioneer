// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Ship.h"
#include "AsteroidBody.h"
#include "Game.h"
#include "Pi.h"
#include "Serializer.h"
#include "Sfx.h"
#include "Space.h"
#include "EnumStrings.h"
#include "LuaTable.h"
#include "collider/collider.h"
#include "scenegraph/SceneGraph.h"
#include "scenegraph/ModelSkin.h"

void AsteroidBody::SaveToJson(Json::Value &jsonObj, Space *space)
{
	DynamicBody::SaveToJson(jsonObj, space);

	Json::Value asteroidBodyObj(Json::objectValue); // Create JSON object to contain dynamic body data.

	asteroidBodyObj["hit_points"] = FloatToStr(m_hitpoints);

	jsonObj["asteroid_body"] = asteroidBodyObj; // Add asteroid body object to supplied object.
}

void AsteroidBody::LoadFromJson(const Json::Value &jsonObj, Space *space)
{
	DynamicBody::LoadFromJson(jsonObj, space);
	Init();

	if (!jsonObj.isMember("asteroid_body")) throw SavedGameCorruptException();
	Json::Value asteroidBodyObj = jsonObj["asteroid_body"];

	if (!asteroidBodyObj.isMember("hit_points")) throw SavedGameCorruptException();
	m_hitpoints = StrToFloat(asteroidBodyObj["hit_points"].asString());
}

void AsteroidBody::Init()
{
	m_hitpoints = 1.0f;
	//SetLabel(ScopedTable(m_cargo).CallMethod<std::string>("GetName"));
	SetMassDistributionFromModel();

	std::vector<Color> colors;
	//metallic blue-orangeish color scheme
	colors.push_back(Color(255, 198, 64));
	colors.push_back(Color(0, 222, 255));
	colors.push_back(Color(255, 255, 255));

	SceneGraph::ModelSkin skin;
	skin.SetColors(colors);
	skin.SetDecal("pioneer");
	skin.Apply(GetModel());
	GetModel()->SetColors(colors);
}

AsteroidBody::AsteroidBody()
{
	SetModel("cargo");
	Init();
	SetMass(1.0);
}

void AsteroidBody::TimeStepUpdate(const float timeStep)
{
	// Suggestion: since asteroids don't need thrust or AI, it could be
	// converted into an idle object on orbital rails, set up to only take
	// memory & save file space (not CPU power) when far from the player.
	
	DynamicBody::TimeStepUpdate(timeStep);
}

bool AsteroidBody::OnDamage(Object *attacker, float kgDamage, const CollisionContact& contactData)
{
	m_hitpoints -= kgDamage*0.001f;
	if (m_hitpoints < 0) {
		Pi::game->GetSpace()->KillBody(this);
		Sfx::Add(this, Sfx::TYPE_EXPLOSION);
	}
	return true;
}

bool AsteroidBody::OnCollision(Object *b, Uint32 flags, double relVel)
{
	// ignore collision if its about to be scooped
	if (b->IsType(Object::SHIP)) {
		int cargoscoop_cap = 0;
		static_cast<Ship*>(b)->Properties().Get("cargo_scoop_cap", cargoscoop_cap);
		if (cargoscoop_cap > 0)
			return true;
	}

	return DynamicBody::OnCollision(b, flags, relVel);
}

void AsteroidBody::Render(Graphics::Renderer *r, const Camera *camera, const vector3d &viewCoords, const matrix4x4d &viewTransform)
{
	RenderModel(r, camera, viewCoords, viewTransform);
}

void AsteroidBody::SetLabel(const std::string &label)
{
	assert(GetModel());
	GetModel()->SetLabel(label);
	Body::SetLabel(label);
}
