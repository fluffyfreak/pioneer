// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifdef VERTEXCOLOR
out vec4 vertexColor;
#endif
#if (NUM_LIGHTS > 0)
out vec3 eyePos;
out vec3 normal;
out vec3 wNormal;
out vec3 wCoords;
#ifdef MAP_NORMAL
out vec3 tangent;
out vec3 bitangent;
#endif //MAP_NORMAL
#ifdef HEAT_COLOURING
uniform mat3 heatingMatrix;
uniform vec3 heatingNormal; // normalised
out vec3 heatingDir;
#endif // HEAT_COLOURING
#endif // (NUM_LIGHTS > 0)
uniform float scale;

void main(void)
{
	gl_Position = logarithmicTransform();
#ifdef VERTEXCOLOR
	vertexColor = a_color;
#endif
#if (NUM_LIGHTS > 0)
	eyePos = vec3(uViewMatrix * a_vertex);
	normal = normalize(uNormalMatrix * vec4(a_normal, 1.0)).xyz;
	wNormal = a_normal;
	wCoords = a_vertex.xyz * scale;
#ifdef MAP_NORMAL
	tangent = normalize(uNormalMatrix * vec4(a_tangent, 1.0)).xyz;
	bitangent = normalize(uNormalMatrix * vec4(cross(a_normal, a_tangent), 1.0)).xyz;
#endif
#ifdef HEAT_COLOURING
	heatingDir = normalize(heatingMatrix * heatingNormal);
#endif
#endif
}
