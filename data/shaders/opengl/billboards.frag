// Copyright © 2008-2017 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "attributes.glsl"
#include "logz.glsl"
#include "lib.glsl"

uniform sampler2D texture0;
#ifdef USE_SPRITE_ATLAS
uniform float coordDownScale;
in vec2 uv0;
in vec2 uv1;
in float tBlend;
#endif

out vec4 frag_color;

void main(void)
{
#ifdef USE_SPRITE_ATLAS
	frag_color = mix(	texture(texture0, (gl_PointCoord * coordDownScale) + uv0),
						texture(texture0, (gl_PointCoord * coordDownScale) + uv1),
						tBlend);
#else
	frag_color = texture(texture0, gl_PointCoord);
#endif
	SetFragDepth();
}
