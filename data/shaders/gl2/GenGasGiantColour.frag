// Copyright © 2008-2014 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

uniform vec3 v0;
uniform vec3 v1;
uniform vec3 v2;
uniform vec3 v3;
uniform float fracStep;

varying vec3 vertex;
varying vec2 uv;

#if defined( GEN_JUPITER_ESQUE )
vec4 GetColour(const vec3 p)
{
	vec4 colour = vec4(1.0, 0.0, 0.0, 1.0);
	return colour;
}

/*
uniform float time;

uniform int octaves[10];
uniform float lacunarity[10];
uniform float frequency[10];

uniform vec3 ggdarkColor[8];
uniform vec3 gglightColor[8];
uniform float entropy;
uniform float planetEarthRadii;

vec4 GetColour(const vec3 p)
{
	float n;
	float sn = snoise(vec4(p.x*8.0, p.y*32.0, p.z*8.0, time));
	float h = river_octavenoise(octaves[0], 0.5*entropy + 0.25, lacunarity[0], vec3(sn,sn,sn), frequency[0], 1.0) * 0.125;
	float equatorial_region_1 = billow_octavenoise(octaves[0], 0.8, lacunarity[0], p, frequency[0], 1.0) * p.y * p.x;
	float equatorial_region_2 = octavenoise(octaves[1], 0.8, lacunarity[1], p, frequency[1], 1.0) * p.x * p.x;
	vec3 col;
	col = mix(ggdarkColor[0], ggdarkColor[1], equatorial_region_1);
	col = mix(col, vec3(0.45, 0.3, 0.0), equatorial_region_2);
	vec3 pRadii03 = snoise(vec4(p.x, p.y*planetEarthRadii*0.3, p.z, time))*p;
	vec3 pRadii = snoise(vec4(p.x, p.y*planetEarthRadii, p.z, time))*p;

	// default stripes
	float IncLoop = 0.3;
	float entScale1 = 0.6;
	float entScale2 = 0.6;
	float entScale3 = 0.7;
	float compMin = 0.1+h;
	float compMax = -0.0+h;

	if (p.y < 0.5 && p.y > 0.1) {
		// top stripe ?
		IncLoop = 0.6;
		entScale1 = 0.7;
		entScale3 = 0.6;
		compMin = 0.15+h;
		compMax = -0.15+h;
	} else if (p.y < -0.1 && p.y > -0.5) {
		// bottom stripe ?
		IncLoop = 0.6;
		entScale2 = 0.7;
		entScale3 = 0.6;
		compMin = 0.15+h;
		compMax = -0.15+h;
	}

	vec3 colArrayStart[5];
	colArrayStart[0] = col;
	colArrayStart[1] = gglightColor[4];
	colArrayStart[2] = vec3(0.9, 0.89, 0.85);
	colArrayStart[3] = ggdarkColor[2];
	colArrayStart[4] = col;

	vec3 colArrayEnd[5];
	colArrayEnd[0] = ggdarkColor[7];
	colArrayEnd[1] = col;
	colArrayEnd[2] = gglightColor[4];
	colArrayEnd[3] = vec3(0.9, 0.89, 0.85);
	colArrayEnd[4] = ggdarkColor[2];

	for(float i=-1.0 ; i < 1.0; i+=IncLoop) {
		float temp = p.y - i;
		if ( temp < compMin && temp > compMax ) {
			n = billow_octavenoise(octaves[2], entScale1*entropy, lacunarity[2], pRadii03, frequency[2], time);
			n += 0.5*octavenoise(octaves[1], entScale2*entropy, lacunarity[1], pRadii, frequency[1], time);
			n += ridged_octavenoise(octaves[1], entScale3*entropy, lacunarity[1], pRadii03, frequency[1], time);
			n = (n<0.0 ? -n : n);
			n = (n>1.0 ? 2.0-n : n);

			float fltIdx = floor((n*0.5)*10.0); // 0 -> 4
			int idx = int(fltIdx);
			float subVal = ((fltIdx * 0.1) * 2.0);
			n -= subVal; n *= 5.0;
			col = mix( colArrayStart[idx], colArrayEnd[idx], n);
			return col;
		}
	}

	//if is not a stripe.
	n = octavenoise(octaves[1], 0.6*entropy + 0.25, lacunarity[1], snoise(vec4(p.x, p.y*planetEarthRadii*3.0, p.z, time))*p, frequency[1], time);
	n *= n*n;
	n = (n<0.0 ? -n : n);
	n = (n>1.0 ? 2.0-n : n);

	if (n>0.5) {
		n -= 0.5; n*= 2.0;
		col = mix(col, gglightColor[2], n);
		return col;
	} 

	n *= 2.0;
	col = mix(vec3(0.9, 0.89, 0.85), col, n);
	return col;
}*/
#elif defined( GEN_NEPTUNE_ESQUE )
vec4 GetColour(const vec3 p)
{
	vec4 colour = vec4(0.0, 1.0, 0.0, 1.0);
	return colour;
}
#elif defined( GEN_SATURN_ESQUE )
vec4 GetColour(const vec3 p)
{
	vec4 colour = vec4(0.0, 0.0, 1.0, 1.0);
	return colour;
}
#elif defined( GEN_URANUS_ESQUE )
vec4 GetColour(const vec3 p)
{
	vec4 colour = vec4(1.0, 0.0, 1.0, 1.0);
	return colour;
}
#endif

// in patch surface coords, [0,1]
// v[0] to v[3] are the corner vertices
vec3 GetSpherePoint(const float x, const float y) {
	return normalize(v0 + x*(1.0-y)*(v1-v0) + x*y*(v2-v0) + (1.0-x)*y*(v3-v0));
}

void main(void)
{
	float xfrac = (uv.x-0.5) * fracStep;
	float yfrac = (uv.y-0.5) * fracStep;
	vec3 p = GetSpherePoint(xfrac, yfrac);
	
	// call the GetColour function implemented for this shader type
	vec4 colour = GetColour(p);
	
	// success! (hopefully)
	gl_FragColor = vec4(1.0, 0.0, 0.0, 1.0);//colour;
}
