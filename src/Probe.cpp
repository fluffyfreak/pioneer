// Copyright © 2008-2014 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Probe.h"
#include "Serializer.h"
#include "Space.h"
#include "Sfx.h"
#include "ShipType.h"
#include "Lang.h"
#include "Pi.h"
#include "Game.h"
#include "LuaEvent.h"

Probe::Probe(ShipType::Id shipId, Body *owner): Ship(shipId)
{
	m_owner = owner;
	SetLabel(Lang::MISSILE);
}

void Probe::PostLoadFixup(Space *space)
{
	Ship::PostLoadFixup(space);
	m_owner = space->GetBodyByIndex(m_ownerIndex);
}

void Probe::Save(Serializer::Writer &wr, Space *space)
{
	Ship::Save(wr, space);
	wr.Int32(space->GetIndexForBody(m_owner));
}

void Probe::Load(Serializer::Reader &rd, Space *space)
{
	Ship::Load(rd, space);
	m_ownerIndex = rd.Int32();
}

void Probe::TimeStepUpdate(const float timeStep)
{
	Ship::TimeStepUpdate(timeStep);

	const float MISSILE_DETECTION_RADIUS = 100.0f;
	if (!m_owner) {
		Explode();
	}/* else if (m_armed) {
		Space::BodyNearList nearby;
		Pi::game->GetSpace()->GetBodiesMaybeNear(this, MISSILE_DETECTION_RADIUS, nearby);
		for (Space::BodyNearIterator i = nearby.begin(); i != nearby.end(); ++i) {
			if (*i == this) continue;
			double dist = ((*i)->GetPosition() - GetPosition()).Length();
			if (dist < MISSILE_DETECTION_RADIUS) {
				Explode();
				break;
			}
		}
	}*/
}

bool Probe::OnCollision(Object *o, Uint32 flags, double relVel)
{
	if (!IsDead()) {
		Explode();
	}
	return true;
}

bool Probe::OnDamage(Object *attacker, float kgDamage, const CollisionContact& contactData)
{
	if (!IsDead()) {
		Explode();
	}
	return true;
}

void Probe::Explode()
{
	Pi::game->GetSpace()->KillBody(this);

	const double damageRadius = 200.0;
	const double kgDamage = 10000.0;

	CollisionContact dummy;
	Space::BodyNearList nearby;
	Pi::game->GetSpace()->GetBodiesMaybeNear(this, damageRadius, nearby);
	for (Space::BodyNearIterator i = nearby.begin(); i != nearby.end(); ++i) {
		if ((*i)->GetFrame() != GetFrame()) continue;
		double dist = ((*i)->GetPosition() - GetPosition()).Length();
		if (dist < damageRadius) {
			// linear damage decay with distance
			(*i)->OnDamage(m_owner, kgDamage * (damageRadius - dist) / damageRadius, dummy);
			if ((*i)->IsType(Object::SHIP))
				LuaEvent::Queue("onShipHit", dynamic_cast<Ship*>(*i), m_owner);
		}
	}

	Sfx::Add(this, Sfx::TYPE_EXPLOSION);
}

void Probe::NotifyRemoved(const Body* const removedBody)
{
	if (m_owner == removedBody) {
		m_owner = 0;
	}
}
