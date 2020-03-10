// Copyright © 2008-2020 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "attributes.glsl"
#include "logz.glsl"
#include "lib.glsl"

out vec4 vertexColor;
out vec3 varyingTexCoord0;

void main(void)
{
	gl_Position = logarithmicTransform();
	vertexColor = a_color;
	varyingTexCoord0 = a_normal.xyz;
}
