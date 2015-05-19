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

#ifdef GEN_JUPITER_ESQUE
vec4 GetColour(in vec3 p)
{	
	//float n = octavenoise(8, 0.5, 2.0, p, frequency[0], time);
	//float n = river_octavenoise(8, 0.5, 2.0, p, frequency[0], time);
	//float n = dunes_octavenoise(8, 0.5, 2.0, p, frequency[0], time);
	//float n = billow_octavenoise(8, 0.5, 2.0, p, frequency[0], time);
	//float n = ridged_octavenoise(8, 0.5, 2.0, p, frequency[0], time);
	//float n = combo_octavenoise(8, 0.5, 2.0, p, frequency[0], time);
	//float n = voronoiscam_octavenoise(8, 0.5, 2.0, p, frequency[0], time);
	//float n = megavolcano_function(8, frequency[0], amplitude[0], p);
	//vec2 F = cellular( p * 16.0 );
	//float n = F.y-F.x;
	//return mix(vec4(0.04, 0.05, 0.15, 1.0), vec4(0.80,0.94,0.96, 1.0), n);
	
	// spots
	//vec2 F = cellular( p * 4.0 );
	//float s = fwidth(F.x);
	//float n1 = smoothstep(0.4-s, 0.4+s, F.x);
	//float n2 = smoothstep(0.5-s, 0.5+s, F.x);
	//return vec4(n1, n2, n2, 1.0);
	
	// partial spot
	//vec2 F = cellular( p * 4.0 );
	//float s = fwidth(F.x);
	//float ss = smoothstep(0.5-s, 0.5+s, F.x);
	
	
	//float n1 = octavenoise(8, 0.5, 2.0, p, frequency[0], time);
	//float n2 = octavenoise(8, 0.5, 2.0, p * 3.14159, frequency[2], time);
	//vec4 color = vec4(texture(texture0, vec2(0.0, (p.y + 1.0) * 0.5) + vec2(n1*0.075,n2*0.075)).xyz, 1.0);
	//return color;

	// fractal definition example:
	// float #fractalNoise#(vec3 position, int octaves, float frequency, float persistence, float time);
	// -------------------------------------------------------------------------------------------------
	float n1 = fbm(p, 6, 0.1, 0.8, time) * 0.01;
	float n2 = ridgedfbm(p, 5, 5.8, 0.75, time) * 0.015 - 0.01;

	// Get the three threshold samples
	float s = 0.6;
	float t1 = snoise(vec4(p * 2.0, time)) - s;
	float t2 = snoise(vec4((p + 800.0) * 2.0, time)) - s;
	float t3 = snoise(vec4((p + 1600.0) * 2.0, time)) - s;
 
	// Intersect them and get rid of negatives
	float threshold = max(t1 * t2 * t3, 0.0);

	// Storms
	float n3 = snoise(vec4(p, time)) * threshold * 3.0;
	float n = n1 + n2 + n3;
	
	vec2 newUV = vec2(n,n);
	vec3 texColor = texture(texture0, vec2(0.0, (p.y + 1.0) * 0.5) + newUV).xyz;
 
	// Add to red color channel for debugging
	//return vec4(threshold * 3.0, 0.0, 0.0, 0.0) + vec4(texColor, 1.0);
	return vec4(texColor, 1.0);
}
#endif

#ifdef GEN_SATURN_ESQUE
vec4 GetColour(in vec3 p)
{
	float n1 = octavenoise(8, 0.5, 2.0, p, frequency[0], time);
	float n2 = octavenoise(8, 0.5, 2.0, p * 3.14159, frequency[2], time);
	vec4 color = vec4(texture(texture0, vec2(0.0, (p.y + 1.0) * 0.5) + vec2(n1*0.075,n2*0.075)).xyz, 1.0);
	return color;
}
#endif

#ifdef GEN_SATURN2_ESQUE
vec4 GetColour(in vec3 p)
{
	float n1 = octavenoise(8, 0.5, 2.0, p, frequency[0], time);
	float n2 = octavenoise(8, 0.5, 2.0, p * 3.14159, frequency[2], time);
	vec4 color = vec4(texture(texture0, vec2(0.0, (p.y + 1.0) * 0.5) + vec2(n1*0.075,n2*0.075)).xyz, 1.0);
	return color;
}
#endif // GEN_SATURN2_ESQUE

#ifdef GEN_NEPTUNE_ESQUE
vec4 GetColour(in vec3 p)
{
	float n1 = octavenoise(8, 0.5, 2.0, p, frequency[0], time);
	float n2 = octavenoise(8, 0.5, 2.0, p * 3.14159, frequency[2], time);
	vec4 color = vec4(texture(texture0, vec2(0.0, (p.y + 1.0) * 0.5) + vec2(n1*0.075,n2*0.075)).xyz, 1.0);
	return color;
}
#endif

#ifdef GEN_NEPTUNE2_ESQUE
vec4 GetColour(in vec3 p)
{
	float n1 = octavenoise(8, 0.5, 2.0, p, frequency[0], time);
	float n2 = octavenoise(8, 0.5, 2.0, p * 3.14159, frequency[2], time);
	vec4 color = vec4(texture(texture0, vec2(0.0, (p.y + 1.0) * 0.5) + vec2(n1*0.075,n2*0.075)).xyz, 1.0);
	return color;
}
#endif

#ifdef GEN_URANUS_ESQUE 
vec4 GetColour(in vec3 p)
{
	float n1 = octavenoise(8, 0.5, 2.0, p, frequency[0], time);
	float n2 = octavenoise(8, 0.5, 2.0, p * 3.14159, frequency[2], time);
	vec4 color = vec4(texture(texture0, vec2(0.0, (p.y + 1.0) * 0.5) + vec2(n1*0.075,n2*0.075)).xyz, 1.0);
	return color;
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
