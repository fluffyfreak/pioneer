// Copyright © 2008-2016 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "attributes.glsl"
#include "logz.glsl"
#include "lib.glsl"

out vec3 varyingEyepos;
out vec3 varyingNormal;
out vec4 vertexColor;

uniform vec3 geosphereCenter;
uniform float geosphereRadius;

#ifdef DETAIL_MAPS
out vec2 texCoord0;
out float dist;
#endif // DETAIL_MAPS

void main(void)
{
	gl_Position = logarithmicTransform();
	vertexColor = a_color;
	varyingEyepos = vec3(uViewMatrix * a_vertex);
	varyingNormal = normalize(uNormalMatrix * a_normal);
	
#ifdef DETAIL_MAPS
	texCoord0 = a_uv0.xy;
	dist = abs(varyingEyepos.z);
#endif // DETAIL_MAPS
}
