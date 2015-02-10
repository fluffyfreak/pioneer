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

void AsteroidBody::Save(Serializer::Writer &wr, Space *space)
{
	DynamicBody::Save(wr, space);
}

void AsteroidBody::Load(Serializer::Reader &rd, Space *space)
{
	DynamicBody::Load(rd, space);
	Init();
}

void AsteroidBody::Init()
{
	SetMassDistributionFromModel();
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
