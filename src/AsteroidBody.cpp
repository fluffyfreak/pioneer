// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Ship.h"
#include "AsteroidBody.h"
#include "Game.h"
#include "GameSaveError.h"
#include "Pi.h"
#include "Sfx.h"
#include "Space.h"
#include "EnumStrings.h"
#include "collider/collider.h"
#include "scenegraph/SceneGraph.h"
#include "scenegraph/ModelSkin.h"

AsteroidBody::AsteroidBody()
{
	SetModel("cargo");
	Init();
	SetMass(1.0);
}

AsteroidBody::AsteroidBody(const Json &jsonObj, Space *space) :
	DynamicBody(jsonObj, space)
{
	Init();

	try {
		Json asteroidBodyObj = jsonObj["asteroid_body"];
		m_hitpoints = StrToFloat(asteroidBodyObj["hit_points"]);
	} catch (Json::type_error &) {
		throw SavedGameCorruptException();
	}

}

void AsteroidBody::SaveToJson(Json &jsonObj, Space *space)
{
	DynamicBody::SaveToJson(jsonObj, space);

	Json asteroidBodyObj = Json::object(); // Create JSON object to contain dynamic body data.

	asteroidBodyObj["hit_points"] = FloatToStr(m_hitpoints);

	jsonObj["asteroid_body"] = asteroidBodyObj; // Add asteroid body object to supplied object.
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
		SfxManager::Add(this, TYPE_EXPLOSION);
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
