// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifdef VERTEXCOLOR
out vec4 vertexColor;
#endif
out vec2 texCoord0;
#if (NUM_LIGHTS > 0)
out vec3 eyePos;
out vec3 normal;
out vec3 wNormal;
out vec3 wCoords01;
out vec3 wCoords23;
#ifdef MAP_NORMAL
out vec3 tangent;
out vec3 bitangent;
#endif
#ifdef HEAT_COLOURING
uniform mat3 heatingMatrix;
uniform vec3 heatingNormal; // normalised
out vec3 heatingDir;
#endif // HEAT_COLOURING
#endif // (NUM_LIGHTS > 0)
uniform float texScale01;
uniform float texScale23;

void main(void)
{
	gl_Position = logarithmicTransform();
#ifdef VERTEXCOLOR
	vertexColor = a_color;
#endif
	texCoord0 = a_uv0.xy;
#if (NUM_LIGHTS > 0)
	eyePos = vec3(uViewMatrix * a_vertex);
	normal = normalize(uNormalMatrix * a_normal).xyz;
	wNormal = a_normal;
	wCoords01 = a_vertex.xyz * texScale01;
	wCoords23 = a_vertex.xyz * texScale23;
#ifdef MAP_NORMAL
	tangent = normalize(uNormalMatrix * a_tangent.xyz).xyz;
	bitangent = normalize(uNormalMatrix * cross(a_normal, a_tangent.xyz)).xyz;
#endif
#ifdef HEAT_COLOURING
	heatingDir = normalize(heatingMatrix * heatingNormal);
#endif
#endif
}