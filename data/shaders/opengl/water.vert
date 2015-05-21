// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

out vec3 varyingEyepos;
out vec3 varyingNormal;
out vec4 vertexColor;

uniform vec3 geosphereCenter;
uniform float geosphereScaledRadius;

void main(void)
{
	gl_Position = logarithmicTransform();
	vertexColor = a_color;
	varyingEyepos = vec3(uViewMatrix * a_vertex);
	varyingNormal = normalize(uNormalMatrix * a_normal);
}
