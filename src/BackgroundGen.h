// Copyright © 2008-2013 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _BACKGROUND_GEN_H
#define _BACKGROUND_GEN_H

#include "Background.h"

class BackgroundGen {
public:
	BackgroundGen(Graphics::Renderer *r, int width, int height);
	~BackgroundGen();
	void Draw();

private:
	int m_width, m_height;
	float m_aspectRatio;
	Graphics::Renderer *m_renderer;
	std::unique_ptr<Background::Container> m_background;
};

#endif
