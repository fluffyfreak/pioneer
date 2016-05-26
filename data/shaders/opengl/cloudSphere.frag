// Copyright © 2008-2016 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "attributes.glsl"
#include "logz.glsl"
#include "lib.glsl"
#include "eclipse.glsl"
#include "noise.glsl"

uniform vec4 atmosColor;
// to keep distances sane we do a nearer, smaller scam. this is how many times
// smaller the geosphere has been made
uniform float geosphereRadius;
uniform float geosphereAtmosTopRad;
uniform vec3 geosphereCenter;
uniform float geosphereAtmosFogDensity;
uniform float geosphereAtmosInvScaleHeight;

uniform Material material;
uniform Scene scene;

in vec3 varyingEyepos;
in vec3 varyingNormal;

out vec4 frag_color;

void main(void)
{
	vec3 eyepos = varyingEyepos;
	vec3 eyenorm = normalize(eyepos);
	vec3 tnorm = normalize(varyingNormal);
	
	// generate some noise clouds
	const float nScale = 1.4; // Uniform?
	const float Density = 0.02;
		
	vec3 noisePosition = v_texCoord3D * 10.0;
	float noise = fbm(noisePosition, 8, 8.0, 0.5) * nScale;
	float rnoise = ridgedNoise(noisePosition, 4, 1.0, 0.5) * nScale;
	rnoise -= (1.0 - Density);
	
	float thickness = max(rnoise * 2.0 + noise, 0.0);
	thickness *= thickness;
	vec4 texColor = vec4(max(vec3(1.0, 1.0, 1.0) * thickness, 0.0) * 2.0, clamp(thickness, 0.0, 1.0));
	// end of noise clouds
	
	vec4 diff = texColor;
	float nDotVP=0.0;
	float nnDotVP=0.0;

	vec3 v = (eyepos - geosphereCenter)/geosphereRadius;
	float lenInvSq = 1.0/(dot(v,v));
	for (int i=0; i<NUM_LIGHTS; ++i) {
		float uneclipsed = clamp(calcUneclipsed(i, v, normalize(vec3(uLight[i].position))), 0.0, 1.0);
		nDotVP  = max(0.0, dot(tnorm, normalize(vec3(uLight[i].position))));
		nnDotVP = max(0.0, dot(tnorm, normalize(-vec3(uLight[i].position)))); //need backlight to increase horizon
		diff += uLight[i].diffuse * uneclipsed * 0.5*(nDotVP+0.5*clamp(1.0-nnDotVP*4.0,0.0,1.0) * INV_NUM_LIGHTS);
	}

	// when does the eye ray intersect atmosphere
	float atmosStart = findSphereEyeRayEntryDistance(geosphereCenter, eyepos, geosphereRadius * geosphereAtmosTopRad);
	float ldprod=0.0;
	float fogFactor=0.0;
	{
		float atmosDist = (length(eyepos) - atmosStart);
		
		// a&b scaled so length of 1.0 means planet surface.
		vec3 a = (atmosStart * eyenorm - geosphereCenter) / geosphereRadius;
		vec3 b = (eyepos - geosphereCenter) / geosphereRadius;
		ldprod = AtmosLengthDensityProduct(a, b, atmosColor.w*geosphereAtmosFogDensity, atmosDist, geosphereAtmosInvScaleHeight);
		fogFactor = clamp( 1.5 / exp(ldprod),0.0,1.0); 
	}

	//calculate sunset tone red when passing through more atmosphere, clamp everything.
	float atmpower = (diff.r+diff.g+diff.b)/3.0;
	vec4 sunset = vec4(0.8,clamp(pow(atmpower,0.8),0.0,1.0),clamp(pow(atmpower,1.2),0.0,1.0),1.0);
	
	frag_color =
		material.emission +
		fogFactor *
		(scene.ambient + diff) +
		(1.0-fogFactor)*(diff*atmosColor) +
		(0.02-clamp(fogFactor,0.0,0.01))*diff*ldprod*sunset +	      //increase fog scatter				
		(pow((1.0-pow(fogFactor,0.75)),256.0)*0.4*diff*atmosColor)*sunset;  //distant fog.
	
	// The alpha channel is a decoy! Use the red to determine alpha.
	frag_color.a = texColor.r;
	
	SetFragDepth();
}
