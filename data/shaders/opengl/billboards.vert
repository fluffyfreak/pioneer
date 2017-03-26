// Copyright © 2008-2017 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "attributes.glsl"
#include "logz.glsl"
#include "lib.glsl"

#ifdef USE_SPRITE_ATLAS
out vec2 uv0;
out vec2 uv1;
out float tBlend;
#endif

void main(void)
{
	gl_Position = logarithmicTransform();
	gl_PointSize = a_normal.z;
#ifdef USE_SPRITE_ATLAS
	uv0 = a_uv0.xy;
	uv1 = a_normal.xy;
	tBlend = a_uv0.z;
#endif
}
