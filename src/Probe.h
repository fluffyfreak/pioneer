// Copyright © 2008-2014 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#pragma once

#ifndef _PROBE_H
#define _PROBE_H

#include <list>
#include "libs.h"
#include "Ship.h"

class Probe: public Ship {
public:
	OBJDEF(Probe, Ship, PROBE);
	Probe(ShipType::Id type, Body *owner);
	Probe() {}
	virtual ~Probe() {}
	void TimeStepUpdate(const float timeStep);
	virtual bool OnCollision(Object *o, Uint32 flags, double relVel);
	virtual bool OnDamage(Object *attacker, float kgDamage, const CollisionContact& contactData);
	virtual void NotifyRemoved(const Body* const removedBody);
	virtual void PostLoadFixup(Space *space);
	Body *GetOwner() const { return m_owner; }

protected:
	virtual void Save(Serializer::Writer &wr, Space *space);
	virtual void Load(Serializer::Reader &rd, Space *space);
private:
	void Explode();

	Body *m_owner;

	int m_ownerIndex; // deserialisation
};

#endif /* _PROBE_H */
