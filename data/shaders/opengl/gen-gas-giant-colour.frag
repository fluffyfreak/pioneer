// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

uniform sampler2D texture0; // ???

uniform vec3 v0;
uniform vec3 v1;
uniform vec3 v2;
uniform vec3 v3;
uniform float fracStep;

uniform float time;

uniform vec3 frequency;

in vec3 vertex;
in vec2 uv;

out vec4 frag_color;

vec3 NoiseSample(in vec3 p)
{
	// fractal definition example:
	// float #fractalNoise#(vec3 position, int octaves, float frequency, float persistence, float time);
	// -------------------------------------------------------------------------------------------------
	float n1 = fbm      (p, 6, 10.0 * frequency[0], 0.8, time) * 0.01;
	float n2 = ridgedfbm(p, 5, 5.8 * frequency[1], 0.75, time) * 0.015 - 0.01;

	// Get the three threshold samples
	float s = 0.6;
	float t1 = snoise(vec4(p * 2.0, time)) - s;
	float t2 = snoise(vec4((p + 800.0) * 2.0, time)) - s;
	float t3 = snoise(vec4((p + 1600.0) * 2.0, time)) - s;
 
	// Intersect them and get rid of negatives
	float threshold = max(t1 * t2 * t3, 0.0);

	// Storms
	float n3 = snoise(vec4(p * frequency[2], time)) * threshold * 3.0;
	float n = n1 + n2 + n3;
	
	vec2 newUV = vec2(n,n);
	vec3 texColor = texture(texture0, vec2(0.0, (p.y + 1.0) * 0.5) + newUV).xyz;

	return texColor;
}

#ifdef GEN_JUPITER_ESQUE
vec4 GetColour(in vec3 p)
{
	vec3 texColor = NoiseSample(p);
	return vec4(texColor, 1.0);
}
#endif

#ifdef GEN_SATURN_ESQUE
vec4 GetColour(in vec3 p)
{
	vec3 texColor = NoiseSample(p);
	return vec4(texColor, 1.0);
}
#endif

#ifdef GEN_SATURN2_ESQUE
vec4 GetColour(in vec3 p)
{
	vec3 texColor = NoiseSample(p);
	return vec4(texColor, 1.0);
}
#endif // GEN_SATURN2_ESQUE

#ifdef GEN_NEPTUNE_ESQUE
vec4 GetColour(in vec3 p)
{
	vec3 texColor = NoiseSample(p);
	return vec4(texColor, 1.0);
}
#endif

#ifdef GEN_NEPTUNE2_ESQUE
vec4 GetColour(in vec3 p)
{
	vec3 texColor = NoiseSample(p);
	return vec4(texColor, 1.0);
}
#endif

#ifdef GEN_URANUS_ESQUE 
vec4 GetColour(in vec3 p)
{
	vec3 texColor = NoiseSample(p);
	return vec4(texColor, 1.0);
}
#endif

// in patch surface coords, [0,1]
// v[0] to v[3] are the corner vertices
vec3 GetSpherePoint(in float x, in float y) {
	return normalize(v0 + x*(1.0-y)*(v1-v0) + x*y*(v2-v0) + (1.0-x)*y*(v3-v0));
}

void main(void)
{
	float xfrac = (uv.x-0.5) * fracStep;
	float yfrac = (uv.y-0.5) * fracStep;
	vec3 p = GetSpherePoint(xfrac, yfrac);
	
	// call the GetColour function implemented for this shader type
	frag_color = GetColour(p);
	
	SetFragDepth();
}
