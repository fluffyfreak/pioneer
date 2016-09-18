// Copyright © 2008-2013 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _GRAPHICCOLLECTOR_H
#define _GRAPHICCOLLECTOR_H

#include "libs.h"

class Renderable;

// Opaque graphics write to z-buffer and are preferably depth sorted front to back
// Transparent (or blended) graphics do not write to depth buffer and are drawn sorted back to front
// Additive transparent graphics do not write depth, do not need to be sorted, and are drawn last
class RenderableCollector {
public:
	typedef std::vector<Renderable*> GraphicList;

	void AddOpaque(Renderable *g) { m_opaque.push_back(g); }
	void AddTransparent(Renderable *g) { m_transparent.push_back(g); }
	void AddAdditive(Renderable *g) { m_additive.push_back(g); }

	GraphicList::const_iterator BeginOpaque() const { return m_opaque.begin(); }
	GraphicList::const_iterator EndOpaque() const { return m_opaque.end(); }

	GraphicList::const_iterator BeginTransparent() const { return m_transparent.begin(); }
	GraphicList::const_iterator EndTransparent() const { return m_transparent.end(); }

	GraphicList::const_iterator BeginAdditive() const { return m_additive.begin(); }
	GraphicList::const_iterator EndAdditive() const { return m_additive.end(); }

	void Reset() {
		m_opaque.clear();
		m_transparent.clear();
		m_additive.clear();
	}

private:
	GraphicList m_opaque;
	GraphicList m_transparent;
	GraphicList m_additive;
};

#endif
