// Copyright © 2008-2018 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "attributes.glsl"
#include "logz.glsl"
#include "lib.glsl"

out vec3 varyingEyepos;
out vec3 varyingNormal;

uniform vec3 geosphereCenter;
uniform float geosphereRadius;

out vec2 texCoord0;
out vec2 invTexCoord;
out float dist;

void main(void)
{
	gl_Position = logarithmicTransform();
	varyingEyepos = vec3(uViewMatrix * a_vertex);
	varyingNormal = normalize(uNormalMatrix * a_normal);
	
	texCoord0 = a_uv0.xy;
	invTexCoord = vec2(1.0 - a_uv0.x, a_uv0.y);	// nasty hack, needs X inverted for detail texturing but corrected for regular
	dist = abs(varyingEyepos.z);
}
