#ifdef FRAGMENT_SHADER
//scene uniform parameters
struct Scene {
	vec4 ambient;
};

//material uniform parameters
//struct Material {
//	...
//};

//Currently used by: planet ring shader, geosphere shaders
float findSphereEyeRayEntryDistance(in vec3 sphereCenter, in vec3 eyeTo, in float radius)
{
	vec3 v = -sphereCenter;
	vec3 dir = normalize(eyeTo);
	float b = -dot(v, dir);
	float det = (b * b) - dot(v, v) + (radius * radius);
	float entryDist = 0.0;
	if (det > 0.0) {
		det = sqrt(det);
		float i1 = b - det;
		float i2 = b + det;
		if (i2 > 0.0) {
			entryDist = max(i1, 0.0);
		}
	}
	return entryDist;
}

// Used by: geosphere shaders
// Calculate length*density product of a line through the atmosphere
// a - start coord (normalized relative to atmosphere radius)
// b - end coord " "
// centerDensity - atmospheric density at centre of sphere
// length - real length of line in meters
float AtmosLengthDensityProduct(vec3 a, vec3 b, float surfaceDensity, float len, float invScaleHeight)
{
	/* 4 samples */
	float ldprod = 0.0;
	vec3 dir = b-a;
	ldprod = surfaceDensity * (
			exp(-invScaleHeight*(length(a)-1.0)) +
			exp(-invScaleHeight*(length(a + 0.2*dir)-1.0)) +
			exp(-invScaleHeight*(length(a + 0.4*dir)-1.0)) +
			exp(-invScaleHeight*(length(a + 0.6*dir)-1.0)) +
			exp(-invScaleHeight*(length(a + 0.8*dir)-1.0)) +
			exp(-invScaleHeight*(length(b)-1.0)));
	ldprod *= len;
	return ldprod;
}


//
// Description : Array and textureless GLSL 2D/3D/4D simplex
//               noise functions.
//      Author : Ian McEwan, Ashima Arts.
//  Maintainer : ijm
//     Lastmod : 20110822 (ijm)
//     License : Copyright (C) 2011 Ashima Arts. All rights reserved.
//               Distributed under the MIT License. See LICENSE file.
//               https://github.com/ashima/webgl-noise
//

vec4 mod289(vec4 x) {
    return x - floor(x * (1.0 / 289.0)) * 289.0;
}

float mod289(float x) {
    return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec4 permute(vec4 x) {
    return mod289(((x*34.0)+1.0)*x);
}

float permute(float x) {
    return mod289(((x*34.0)+1.0)*x);
}

vec4 taylorInvSqrt(vec4 r)
{
    return 1.79284291400159 - 0.85373472095314 * r;
}

float taylorInvSqrt(float r)
{
    return 1.79284291400159 - 0.85373472095314 * r;
}

vec4 grad4(float j, vec4 ip)
{
    const vec4 ones = vec4(1.0, 1.0, 1.0, -1.0);
    vec4 p,s;

    p.xyz = floor( fract (vec3(j) * ip.xyz) * 7.0) * ip.z - 1.0;
    p.w = 1.5 - dot(abs(p.xyz), ones.xyz);
    s = vec4(lessThan(p, vec4(0.0)));
    p.xyz = p.xyz + (s.xyz*2.0 - 1.0) * s.www;

    return p;
}

float snoise(vec4 v)
{
    const vec4  C = vec4( 0.138196601125011,  // (5 - sqrt(5))/20  G4
                          0.276393202250021,  // 2 * G4
                          0.414589803375032,  // 3 * G4
                          -0.447213595499958); // -1 + 4 * G4

	// (sqrt(5) - 1)/4 = F4, used once below
	#define F4 0.309016994374947451

	// First corner
    vec4 i  = floor(v + dot(v, vec4(F4)) );
    vec4 x0 = v -   i + dot(i, C.xxxx);

	// Other corners

	// Rank sorting originally contributed by Bill Licea-Kane, AMD (formerly ATI)
    vec4 i0;
    vec3 isX = step( x0.yzw, x0.xxx );
    vec3 isYZ = step( x0.zww, x0.yyz );
	//  i0.x = dot( isX, vec3( 1.0 ) );
    i0.x = isX.x + isX.y + isX.z;
    i0.yzw = 1.0 - isX;
	//  i0.y += dot( isYZ.xy, vec2( 1.0 ) );
    i0.y += isYZ.x + isYZ.y;
    i0.zw += 1.0 - isYZ.xy;
    i0.z += isYZ.z;
    i0.w += 1.0 - isYZ.z;

    // i0 now contains the unique values 0,1,2,3 in each channel
    vec4 i3 = clamp( i0, 0.0, 1.0 );
    vec4 i2 = clamp( i0-1.0, 0.0, 1.0 );
    vec4 i1 = clamp( i0-2.0, 0.0, 1.0 );

    //  x0 = x0 - 0.0 + 0.0 * C.xxxx
    //  x1 = x0 - i1  + 1.0 * C.xxxx
    //  x2 = x0 - i2  + 2.0 * C.xxxx
    //  x3 = x0 - i3  + 3.0 * C.xxxx
    //  x4 = x0 - 1.0 + 4.0 * C.xxxx
    vec4 x1 = x0 - i1 + C.xxxx;
    vec4 x2 = x0 - i2 + C.yyyy;
    vec4 x3 = x0 - i3 + C.zzzz;
    vec4 x4 = x0 + C.wwww;

	// Permutations
    i = mod289(i);
    float j0 = permute( permute( permute( permute(i.w) + i.z) + i.y) + i.x);
    vec4 j1 = permute( permute( permute( permute (
            i.w + vec4(i1.w, i2.w, i3.w, 1.0 ))
                                         + i.z + vec4(i1.z, i2.z, i3.z, 1.0 ))
                                + i.y + vec4(i1.y, i2.y, i3.y, 1.0 ))
                       + i.x + vec4(i1.x, i2.x, i3.x, 1.0 ));

	// Gradients: 7x7x6 points over a cube, mapped onto a 4-cross polytope
	// 7*7*6 = 294, which is close to the ring size 17*17 = 289.
    vec4 ip = vec4(1.0/294.0, 1.0/49.0, 1.0/7.0, 0.0) ;

    vec4 p0 = grad4(j0,   ip);
    vec4 p1 = grad4(j1.x, ip);
    vec4 p2 = grad4(j1.y, ip);
    vec4 p3 = grad4(j1.z, ip);
    vec4 p4 = grad4(j1.w, ip);

	// Normalise gradients
    vec4 norm = taylorInvSqrt(vec4(dot(p0,p0), dot(p1,p1), dot(p2, p2), dot(p3,p3)));
    p0 *= norm.x;
    p1 *= norm.y;
    p2 *= norm.z;
    p3 *= norm.w;
    p4 *= taylorInvSqrt(dot(p4,p4));

	// Mix contributions from the five corners
    vec3 m0 = max(0.6 - vec3(dot(x0,x0), dot(x1,x1), dot(x2,x2)), 0.0);
    vec2 m1 = max(0.6 - vec2(dot(x3,x3), dot(x4,x4)            ), 0.0);
    m0 = m0 * m0;
    m1 = m1 * m1;
    return 49.0 * ( dot(m0*m0, vec3( dot( p0, x0 ), dot( p1, x1 ), dot( p2, x2 )))
                    + dot(m1*m1, vec2( dot( p3, x3 ), dot( p4, x4 ) ) ) ) ;
}

float octavenoise(in int octaves, in float roughness, in float lacunarity, in vec3 p, in float freq, in float time)
{
    float n = 0.0;
    float octaveAmplitude = 1.0/(1.0-pow(roughness,(float(octaves))));
    for (int i = 0; i < octaves; i++) {
        n += octaveAmplitude * snoise(vec4(freq*p, time));
        octaveAmplitude *= roughness;
        freq *= lacunarity;
    }
    return (n+1.0)*0.5;
}

float combo_octavenoise(in int octaves, in float roughness, in float lacunarity, in vec3 p, in float freq, in float time)
{
    float n = 0.0;
    float n1 = 0.0;
    float octaveAmplitude = 1.0/(1.0-pow(roughness,(float(octaves))));
    for (int i = 0; i < octaves; i++) {
        n += octaveAmplitude * snoise(vec4(freq*p, time));
        octaveAmplitude *= roughness;
        freq *= lacunarity;
    }
    //ridged noise
    n1 = 1.0 - abs(n);
    n1 *= n1;
    //billow noise
    n1 *= (2.0 * abs(n) - 1.0)+1.0;
    //voronoiscam noise
    n1 *= sqrt(10.0 * abs(n));
    return n1;
}

float ridged_octavenoise(in int octaves, in float roughness, in float lacunarity, in vec3 p, in float freq, in float time)
{
    float n = 0.0;
    float octaveAmplitude = 1.0/(1.0-pow(roughness,(float(octaves))));
    for (int i = 0; i < octaves; i++) {
        n += octaveAmplitude * snoise(vec4(freq*p, time));
        octaveAmplitude *= roughness;
        freq *= lacunarity;
    }
    //ridged noise
    n = 1.0 - abs(n);
    return(n*n);
}

float billow_octavenoise(in int octaves, in float roughness, in float lacunarity, in vec3 p, in float freq, in float time)
{
    float n = 0.0;
    float octaveAmplitude = 1.0/(1.0-pow(roughness,(float(octaves))));
    for (int i = 0; i < octaves; i++) {
        n += octaveAmplitude * snoise(vec4(freq*p, time));
        octaveAmplitude *= roughness;
        freq *= lacunarity;
    }
    //ridged noise
    n = (2.0 * abs(n) - 1.0)+1.0;
    return(n);
}

float dunes_octavenoise(in int octaves, in float roughness, in float lacunarity, in vec3 p, in float freq, in float time) {
	float n = 0.0;
	float octaveAmplitude = roughness;
	for (int i=0; i<octaves; i++) {
		n += octaveAmplitude * snoise(vec4(freq*p,1.0));
		octaveAmplitude *= roughness;
		freq *= lacunarity;
	}
	return 1.0 - abs(n);
}

float river_octavenoise(in int octaves, in float roughness, in float lacunarity, in vec3 p, in float freq, in float time) {
	float n = 0.0;
	float octaveAmplitude = roughness;
	for (int i=0; i<octaves; i++) {
		n += octaveAmplitude * abs(snoise(vec4(freq*p,1.0)));
		octaveAmplitude *= roughness;
		freq *= lacunarity;
	}
	return n;
}

float voronoiscam_octavenoise(in int octaves, in float roughness, in float lacunarity, in vec3 p, in float freq, in float time) {
	float n = 0.0;
	float octaveAmplitude = roughness;
	for (int i=0; i<octaves; i++) {
		n += octaveAmplitude * snoise(vec4(freq*p,1.0));
		octaveAmplitude *= roughness;
		freq *= lacunarity;
	}
	return sqrt(10.0 * abs(n));
}


// Creates small canyons.
float canyon_ridged_function(in int octaves, in float amplitude, in float frequency, in float lacunarity, in vec3 p)
{
	float h;
	float n = ridged_octavenoise(octaves, 0.54, 2.0, p, frequency, 1.0);
	float outer = 0.71;
	float inner = 0.715;
	float inner2 = 0.715;
	float outer2 = 0.72;
	if (n > outer2) {
		h = 1.0;
	} else if (n > inner2) {
		h = 0.0+1.0*(n-inner2)*(1.0/(outer2-inner2));
	} else if (n > inner) {
		h = 0.0;
	} else if (n > outer) {
		h = 1.0-1.0*(n-outer)*(1.0/(inner-outer));
	} else {
		h = 1.0;
	}
	return h * amplitude;
}

// Larger canyon.
float canyon2_ridged_function(in int octaves, in float amplitude, in float frequency, in float lacunarity, in vec3 p)
{
	float h;
	float n = ridged_octavenoise(octaves, 0.56, 2.0, p, frequency, 1.0);
	float outer = 0.7;
	float inner = 0.71;
	float inner2 = 0.72;
	float outer2 = 0.73;
	if (n > outer2) {
		h = 1.0;
	} else if (n > inner2) {
		h = 0.0+1.0*(n-inner2)*(1.0/(outer2-inner2));
	} else if (n > inner) {
		h = 0.0;
	} else if (n > outer) {
		h = 1.0-1.0*(n-outer)*(1.0/(inner-outer));
	} else {
		h = 1.0;
	}
	return h * amplitude;
}

// Largest and best looking canyon, combine them together for best results.
float canyon3_ridged_function(in int octaves, in float amplitude, in float frequency, in float lacunarity, in vec3 p)
{
	float h;
	float n = ridged_octavenoise(octaves, 0.585, 2.0, p, frequency, 1.0);
	float outer = 0.7;
	float inner = 0.71;
	float inner2 = 0.72;
	float outer2 = 0.73;
	if (n > outer2) {
		h = 1.0;
	} else if (n > inner2) {
		h = 0.0+1.0*(n-inner2)*(1.0/(outer2-inner2));
	} else if (n > inner) {
		h = 0.0;
	} else if (n > outer) {
		h = 1.0-1.0*(n-outer)*(1.0/(inner-outer));
	} else {
		h = 1.0;
	}
	return h * amplitude;
}

float canyon_normal_function(in int octaves, in float amplitude, in float frequency, in float lacunarity, in vec3 p)
{
	float h;
	float n = octavenoise(octaves, 0.54, 2.0, p, frequency, 1.0);
	float outer = 0.71;
	float inner = 0.715;
	float inner2 = 0.715;
	float outer2 = 0.72;
	if (n > outer2) {
		h = 1.0;
	} else if (n > inner2) {
		h = 0.0+1.0*(n-inner2)*(1.0/(outer2-inner2));
	} else if (n > inner) {
		h = 0.0;
	} else if (n > outer) {
		h = 1.0-1.0*(n-outer)*(1.0/(inner-outer));
	} else {
		h = 1.0;
	}
	return h * amplitude;
}

float canyon2_normal_function(in int octaves, in float amplitude, in float frequency, in float lacunarity, in vec3 p)
{
	float h;
	float n = octavenoise(octaves, 0.56, 2.0, p, frequency, 1.0);
	float outer = 0.7;
	float inner = 0.71;
	float inner2 = 0.72;
	float outer2 = 0.73;
	if (n > outer2) {
		h = 1.0;
	} else if (n > inner2) {
		h = 0.0+1.0*(n-inner2)*(1.0/(outer2-inner2));
	} else if (n > inner) {
		h = 0.0;
	} else if (n > outer) {
		h = 1.0-1.0*(n-outer)*(1.0/(inner-outer));
	} else {
		h = 1.0;
	}
	return h * amplitude;
}

float canyon3_normal_function(in int octaves, in float amplitude, in float frequency, in float lacunarity, in vec3 p)
{
	float h;
	float n = octavenoise(octaves, 0.585, 2.0, p, frequency, 1.0);
	float outer = 0.7;
	float inner = 0.71;
	float inner2 = 0.72;
	float outer2 = 0.73;
	if (n > outer2) {
		h = 1.0;
	} else if (n > inner2) {
		h = 0.0+1.0*(n-inner2)*(1.0/(outer2-inner2));
	} else if (n > inner) {
		h = 0.0;
	} else if (n > outer) {
		h = 1.0-1.0*(n-outer)*(1.0/(inner-outer));
	} else {
		h = 1.0;
	}
	return h * amplitude;
}

float canyon_voronoi_function(in int octaves, in float amplitude, in float frequency, in float lacunarity, in vec3 p)
{
	float h;
	float n = octavenoise(octaves, 0.54, 2.0, p, frequency, 1.0);
	float outer = 0.71;
	float inner = 0.715;
	float inner2 = 0.715;
	float outer2 = 0.72;
	if (n > outer2) {
		h = 1.0;
	} else if (n > inner2) {
		h = 0.0+1.0*(n-inner2)*(1.0/(outer2-inner2));
	} else if (n > inner) {
		h = 0.0;
	} else if (n > outer) {
		h = 1.0-1.0*(n-outer)*(1.0/(inner-outer));
	} else {
		h = 1.0;
	}
	return h * amplitude;
}

float canyon2_voronoi_function(in int octaves, in float amplitude, in float frequency, in float lacunarity, in vec3 p)
{
	float h;
	float n = octavenoise(octaves, 0.56, 2.0, p, frequency, 1.0);
	float outer = 0.7;
	float inner = 0.71;
	float inner2 = 0.72;
	float outer2 = 0.73;
	if (n > outer2) {
		h = 1.0;
	} else if (n > inner2) {
		h = 0.0+1.0*(n-inner2)*(1.0/(outer2-inner2));
	} else if (n > inner) {
		h = 0.0;
	} else if (n > outer) {
		h = 1.0-1.0*(n-outer)*(1.0/(inner-outer));
	} else {
		h = 1.0;
	}
	return h * amplitude;
}

float canyon3_voronoi_function(in int octaves, in float amplitude, in float frequency, in float lacunarity, in vec3 p)
{
	float h;
	float n = octavenoise(octaves, 0.585, 2.0, p, frequency, 1.0);
	float outer = 0.7;
	float inner = 0.71;
	float inner2 = 0.72;
	float outer2 = 0.73;
	if (n > outer2) {
		h = 1.0;
	} else if (n > inner2) {
		h = 0.0+1.0*(n-inner2)*(1.0/(outer2-inner2));
	} else if (n > inner) {
		h = 0.0;
	} else if (n > outer) {
		h = 1.0-1.0*(n-outer)*(1.0/(inner-outer));
	} else {
		h = 1.0;
	}
	return h * amplitude;
}

float canyon_billow_function(in int octaves, in float amplitude, in float frequency, in float lacunarity, in vec3 p)
{
	float h;
	float n = octavenoise(octaves, 0.54, 2.0, p, frequency, 1.0);
	float outer = 0.71;
	float inner = 0.715;
	float inner2 = 0.715;
	float outer2 = 0.72;
	if (n > outer2) {
		h = 1.0;
	} else if (n > inner2) {
		h = 0.0+1.0*(n-inner2)*(1.0/(outer2-inner2));
	} else if (n > inner) {
		h = 0.0;
	} else if (n > outer) {
		h = 1.0-1.0*(n-outer)*(1.0/(inner-outer));
	} else {
		h = 1.0;
	}
	return h * amplitude;
}

float canyon2_billow_function(in int octaves, in float amplitude, in float frequency, in float lacunarity, in vec3 p)
{
	float h;
	float n = octavenoise(octaves, 0.56, 2.0, p, frequency, 1.0);
	float outer = 0.7;
	float inner = 0.71;
	float inner2 = 0.72;
	float outer2 = 0.73;
	if (n > outer2) {
		h = 1.0;
	} else if (n > inner2) {
		h = 0.0+1.0*(n-inner2)*(1.0/(outer2-inner2));
	} else if (n > inner) {
		h = 0.0;
	} else if (n > outer) {
		h = 1.0-1.0*(n-outer)*(1.0/(inner-outer));
	} else {
		h = 1.0;
	}
	return h * amplitude;
}

float canyon3_billow_function(in int octaves, in float amplitude, in float frequency, in float lacunarity, in vec3 p)
{
	float h;
	float n = octavenoise(octaves, 0.585, 2.0, p, frequency, 1.0);
	float outer = 0.7;
	float inner = 0.71;
	float inner2 = 0.72;
	float outer2 = 0.73;
	if (n > outer2) {
		h = 1.0;
	} else if (n > inner2) {
		h = 0.0+1.0*(n-inner2)*(1.0/(outer2-inner2));
	} else if (n > inner) {
		h = 0.0;
	} else if (n > outer) {
		h = 1.0-1.0*(n-outer)*(1.0/(inner-outer));
	} else {
		h = 1.0;
	}
	return h * amplitude;
}

float crater_function_1pass(in vec3 p, float inp, float height)
{
	float res = inp;
	float n = abs(snoise(vec4(p, 1.0)));
	float ejecta_outer = 0.6;
	float outer = 0.9;
	float inner = 0.94;
	float midrim = 0.93;
	if (n > inner) {
		//res = 0;
	} else if (n > midrim) {
		float hrim = inner - midrim;
		float descent = (hrim-(n-midrim))/hrim;
		res += height * descent * descent;
	} else if (n > outer) {
		float hrim = midrim - outer;
		float ascent = (n-outer)/hrim;
		res += height * ascent * ascent;
	} else if (n > ejecta_outer) {
		// blow down walls of other craters too near this one,
		// so we don't have sharp transition
		//res *= (outer-n)/-(ejecta_outer-outer);
	}
	return res;
}

// makes large and small craters across the entire planet.
float crater_function(in int octaves, in float amplitude, in float frequency, in float lacunarity, in vec3 p)
{
	float crater = 0.0;
	float sz = frequency;
	float max_h = amplitude;
	for (int i=0; i<octaves; i++) {
		crater = crater_function_1pass(sz*p, crater, max_h);
		sz *= 2.0;
		max_h *= 0.5;
	}
	return crater;
}

float impact_crater_function_1pass(in vec3 p, float inp, float height)
{
	float res = inp;
	float n = abs(snoise(vec4(p, 1.0)));
	float ejecta_outer = 0.6;
	float outer = 0.9;
	float midrim = 0.93;
	float hrim;
	float descent;
	if (n > midrim) {
		res -= height;
	} else if (n > outer) {
		hrim = midrim - outer;
		descent = (n-outer)/hrim;
		res -= height * descent * descent;
	} else if (n > ejecta_outer) {
		// blow down walls of other craters too near this one,
		// so we don't have sharp transition
		//res *= (outer-n)/-(ejecta_outer-outer);
	}
	return res;
}

// makes large and small craters across the entire planet.
float impact_crater_function(in int octaves, in float amplitude, in float frequency, in float lacunarity, in vec3 p)
{
	float crater = 0.0;
	float sz = frequency;
	float max_h = amplitude;
	for (int i=0; i<octaves; i++) {
		crater = impact_crater_function_1pass(sz*p, crater, max_h);
		sz *= 2.0;
		max_h *= 0.5;
	}
	return crater;
}

float volcano_function_1pass(in vec3 p, float inp, float height)
{
	float res = inp;
	float n = abs(snoise(vec4(p, 1.0)));
	float ejecta_outer = 0.6;
	float outer = 0.9;
	float inner = 0.975;
	float midrim = 0.971;
	if (n > inner) {
		//res = 0;
	} else if (n > midrim) {
		float hrim = inner - midrim;
		float descent = (hrim-(n-midrim))/hrim;
		res += height * descent;
	} else if (n > outer) {
		float hrim = midrim - outer;
		float ascent = (n-outer)/hrim;
		res += height * ascent * ascent;
	} else if (n > ejecta_outer) {
		// blow down walls of other craters too near this one,
		// so we don't have sharp transition
		res *= (outer-n)/-(ejecta_outer-outer);
	}
	return res;
}

float volcano_function(in int octaves, in float amplitude, in float frequency, in float lacunarity, in vec3 p)
{
	float crater = 0.0;
	float sz = frequency;
	float max_h = amplitude;
	for (int i=0; i<octaves; i++) {
		crater = volcano_function_1pass(sz*p, crater, max_h);
		sz *= 1.0;  //frequency?
		max_h *= 0.4; // height??
	}
	return 3.0 * crater;
}

float megavolcano_function_1pass(in vec3 p, float inp, float height)
{
	float res = inp;
	float n = abs(snoise(vec4(p, 1.0)));
	float ejecta_outer = 0.6;
	float outer = 0.76;  //Radius
	float inner = 0.98;
	float midrim = 0.964;
	if (n > inner) {
		//res = 0;
	} else if (n > midrim) {
		float hrim = inner - midrim;
		float descent = (hrim-(n-midrim))/hrim;
		res += height * descent;
	} else if (n > outer) {
		float hrim = midrim - outer;
		float ascent = (n-outer)/hrim;
		res += height * ascent * ascent;
	} else if (n > ejecta_outer) {
		// blow down walls of other craters too near this one,
		// so we don't have sharp transition
		res *= (outer-n)/-(ejecta_outer-outer);
	}
	return res;
}

float megavolcano_function(in int octaves, in float amplitude, in float frequency, in float lacunarity, in vec3 p)
{
	float crater = 0.0;
	float sz = frequency;
	float max_h = amplitude;
	for (int i=0; i<octaves; i++) {
		crater = megavolcano_function_1pass(sz*p, crater, max_h);
		sz *= 1.0;  //frequency?
		max_h *= 0.15; // height??
	}
	return 4.0 * crater;
}

float river_function(in int octaves, in float amplitude, in float frequency, in float lacunarity, in vec3 p, int style)
{
	float h;
	float n = octavenoise(octaves, 0.585, 2.0, p*0.5, frequency, 1.0);

	float outer = 0.67;
	float inner = 0.715;
	float inner2 = 0.715;
	float outer2 = 0.76;
	if (1==style) {
		outer = 0.01;
		inner = 0.49;
		inner2 = 0.51;
		outer2 = 0.99;
	}

	if (n > outer2) {
		h = 1.0;
	} else if (n > inner2) {
		h = 0.0+1.0*(n-inner2)*(1.0/(outer2-inner2));
	} else if (n > inner) {
		h = 0.0;
	} else if (n > outer) {
		h = 1.0-1.0*(n-outer)*(1.0/(inner-outer));
	} else {
		h = 1.0;
	}
	return h * amplitude;
}

#endif
