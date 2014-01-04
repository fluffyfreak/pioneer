// Copyright © 2008-2013 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _GRAPHICS_TYPES_H
#define _GRAPHICS_TYPES_H

namespace Graphics {

enum LineType {
	LINE_SINGLE = GL_LINES, //draw one line per two vertices
	LINE_STRIP = GL_LINE_STRIP,  //connect vertices
	LINE_LOOP = GL_LINE_LOOP    //connect vertices,  connect start & end
};

//how to treat vertices
enum PrimitiveType {
	TRIANGLES = GL_TRIANGLES,
	TRIANGLE_STRIP = GL_TRIANGLE_STRIP,
	TRIANGLE_FAN = GL_TRIANGLE_FAN,
	POINTS = GL_POINTS
};

enum BlendMode {
	BLEND_SOLID,
	BLEND_ADDITIVE,
	BLEND_ALPHA,
	BLEND_ALPHA_ONE, //"additive alpha"
	BLEND_ALPHA_PREMULT,
	BLEND_SET_ALPHA, // copy alpha channel
	BLEND_DEST_ALPHA // XXX maybe crappy name
};

}

#endif
