// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _ASTEROIDBODY_H
#define _ASTEROIDBODY_H

#include "libs.h"
#include "DynamicBody.h"
#include "lua/LuaRef.h"

namespace Graphics { class Renderer; }

class AsteroidBody: public DynamicBody {
public:
	OBJDEF(AsteroidBody, DynamicBody, ASTEROIDBODY);
	AsteroidBody();
	AsteroidBody(const Json &jsonObj, Space *space);

	virtual void SetLabel(const std::string &label);

	virtual void Render(Graphics::Renderer *r, const Camera *camera, const vector3d &viewCoords, const matrix4x4d &viewTransform);
	virtual void TimeStepUpdate(const float timeStep);

	virtual bool OnCollision(Object *o, Uint32 flags, double relVel);
	virtual bool OnDamage(Object *attacker, float kgDamage, const CollisionContact& contactData);
protected:
	virtual void SaveToJson(Json &jsonObj, Space *space);
private:
	void Init();
	
	float m_hitpoints;
};

#endif /* _ASTEROIDBODY_H */