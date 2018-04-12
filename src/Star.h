// Copyright © 2008-2018 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _STAR_H
#define _STAR_H

#include "TerrainBody.h"
#include "graphics/RenderState.h"

namespace Graphics { class Renderer; }

class Star: public TerrainBody {
public:
	OBJDEF(Star, TerrainBody, STAR);
	Star(SystemBody *sbody);
	Star();
	virtual ~Star() {};

	virtual void Render(Graphics::Renderer *r, Camera *camera, const vector3d &viewCoords, const matrix4x4d &viewTransform) override final;
	virtual void SubRender(Graphics::Renderer *r, const matrix4x4d &modelView, const vector3d &camPos) override final {}
protected:
	void InitStar();
	virtual void LoadFromJson(const Json::Value &jsonObj, Space *space) override;

	Graphics::RenderState *m_haloState;
};

#endif /* _STAR_H */
