// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _ASTEROIDBODY_H
#define _ASTEROIDBODY_H

#include "libs.h"
#include "DynamicBody.h"
#include "LuaRef.h"

namespace Graphics { class Renderer; }

class AsteroidBody: public DynamicBody {
public:
	OBJDEF(AsteroidBody, DynamicBody, ASTEROIDBODY);
	AsteroidBody();

	virtual void SetLabel(const std::string &label);
	virtual void Render(Graphics::Renderer *r, const Camera *camera, const vector3d &viewCoords, const matrix4x4d &viewTransform);
	virtual void TimeStepUpdate(const float timeStep);

protected:
	virtual void Save(Serializer::Writer &wr, Space *space);
	virtual void Load(Serializer::Reader &rd, Space *space);
private:
	void Init();
};

#endif /* _ASTEROIDBODY_H */
