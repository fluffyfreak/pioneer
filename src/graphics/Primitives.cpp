// Copyright © 2008-2012 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Primitives.h"
#include <GL/glew.h>
#include "utils.h"

namespace Graphics {

////////////////////////////////////////////////////////////////
// CGLquad related data and definitions
namespace QuadData {
	static const GLfloat s_vertex_buffer_data[18] = { 
		// triangle 1
		-1.0f,-1.0f, 0.0f, // v1
		 1.0f, 1.0f, 0.0f, // v3
		-1.0f, 1.0f, 0.0f, // v2
		// triangle 2
		 1.0f, 1.0f, 0.0f, // v1
		-1.0f,-1.0f, 0.0f, // v3
		 1.0f,-1.0f, 0.0f  // v2
	};
	static const GLfloat s_normal_buffer_data[18] = { 
		// triangle 1
		 0.0f, 0.0f, 1.0f, // v1
		 0.0f, 0.0f, 1.0f, // v3
		 0.0f, 0.0f, 1.0f, // v2
		// triangle 2
		 0.0f, 0.0f, 1.0f, // v1
		 0.0f, 0.0f, 1.0f, // v3
		 0.0f, 0.0f, 1.0f  // v2
	};
	static const GLfloat s_uv_buffer_data[12] = { 
		// triangle 1
		0.0f, 0.0f, // v1
		1.0f, 1.0f, // v3
		0.0f, 1.0f, // v2
		// triangle 2
		1.0f, 1.0f, // v1
		0.0f, 0.0f, // v3
		1.0f, 0.0f  // v2
	};
}; // namespace QuadData

CGLquad::CGLquad(const bool bNormals_, const bool bUVs_) 
	: mVBO(sizeof(QuadData::s_vertex_buffer_data)/sizeof(QuadData::s_vertex_buffer_data[0]), 
			&QuadData::s_vertex_buffer_data[0], 
			bNormals_	? &QuadData::s_normal_buffer_data[0]	: nullptr,	// normals are optional
			bUVs_		? &QuadData::s_uv_buffer_data[0]		: nullptr) 	// UVcoords are optional
{
}
CGLquad::~CGLquad() 
{
}

void CGLquad::Render() const
{
	mVBO.Bind();

	// Draw the triangle !
	glDrawArrays(GL_TRIANGLES, 0, 2*3); // From index 0 to 2*3 -> 2 triangles

	mVBO.Release();
}

};
