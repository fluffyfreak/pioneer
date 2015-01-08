// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#pragma once

#ifndef _MARCHINGGEOM_H
#define _MARCHINGGEOM_H

//
// Marching Cubes Example Program 
// by Cory Bloyd (corysama@yahoo.com)
//
// A simple, portable and complete implementation of the Marching Cubes
// and Marching Tetrahedrons algorithms in a single source file.
// There are many ways that this code could be made faster, but the 
// intent is for the code to be easy to understand.
//
// For a description of the algorithm go to
// http://astronomy.swin.edu.au/pbourke/modelling/polygonise/
//
// This code is public domain.
//

#include "libs.h"

namespace Graphics {
	class Renderer;
	class Material;
	class RenderState;
	class VertexBuffer;
}

class MarchingGeometry
{
public:
	MarchingGeometry(Graphics::Renderer*);

	void Init();
	void Draw();
	void Update();

private:
	Graphics::Renderer *m_renderer;
	Graphics::RenderState *m_renderState;
	RefCountedPtr<Graphics::Material> m_material;
	RefCountedPtr<Graphics::VertexBuffer> m_vb;

	float m_fTime;
};

#endif
