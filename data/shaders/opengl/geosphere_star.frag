// Copyright © 2008-2020 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "attributes.glsl"
#include "logz.glsl"
#include "lib.glsl"

uniform samplerCube texture0; //diffuse

uniform Material material;

in vec4 vertexColor;
in vec3 varyingTexCoord0;

out vec4 frag_color;

mat3 rotationY( in float angle ) {
	return mat3(	cos(angle),		0,		sin(angle),
			 				0,		1.0,			 0,
					-sin(angle),	0,		cos(angle));
}

void main(void)
{
	float time = uTime * -1.0;
	
	vec3 newUv = varyingTexCoord0 * rotationY(time);
	newUv = normalize(newUv);
	
	// using the texture to offset a 2nd read produces these nice blobby effects,
	// moving it by some offset using time animates it vaguely like plasma on a star
	vec3 texSample 	= texture( texture0, newUv ).rgb;
	float uOff		= (texSample.g * 0.1 * 3.14 + time);
	vec3 starUV		= varyingTexCoord0 * rotationY(uOff) * rotationY(time);
	vec3 starTex	= texture( texture0, starUV ).rgb;
	
	// NUM_LIGHTS == 0 -- unlit rendering - stars
	//emission is used to boost colour of stars, which is a bit odd
	frag_color = vec4(starTex, 1.0) + (material.emission * 0.5);

	SetFragDepth();
}
