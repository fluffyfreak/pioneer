// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "attributes.glsl"
#include "logz.glsl"
#include "lib.glsl"

uniform sampler2D texture0; //diffuse + intensity
uniform sampler2D texture1; //normal(enc) + specular + AO
uniform sampler2D texture2; //diffuse + intensity
uniform sampler2D texture3; //normal(enc) + specular + AO
in vec2 texCoord0;

#if (NUM_LIGHTS > 0)
in vec3 eyePos;
in vec3 normal;
in vec3 wNormal;
in vec3 wCoords01;
in vec3 wCoords23;
in vec3 tangent;
in vec3 bitangent;
#endif // (NUM_LIGHTS > 0)

uniform Scene scene;
uniform Material material;

out vec4 frag_color;

vec3 blend(in vec4 texB, in vec4 texA, in float t)
{
    float depth = 0.2;
    float ma = max(texA.w + t, texB.w + (1.0 - t)) - depth;

    float b1 = max(texA.w + t - ma, 0);
    float b2 = max(texB.w + (1.0 - t) - ma, 0);

    return (texA.xyz * b1 + texB.xyz * b2) / (b1 + b2);
}

vec4 getTriplanarTex(in vec3 blending, in vec3 coords, in sampler2D sampler)
{
	vec4 xaxis = texture( sampler, coords.yz);
	vec4 yaxis = texture( sampler, coords.xz);
	vec4 zaxis = texture( sampler, coords.xy);
	// blend the results of the 3 planar projections.
	vec4 tex = (xaxis * blending.x) + (yaxis * blending.y) + (zaxis * blending.z);
	return tex;
}

#if (NUM_LIGHTS > 0)
//ambient, diffuse, specular
//would be a good idea to make specular optional
void ads(in float spec, in int lightNum, in vec3 pos, in vec3 n, inout vec4 light, inout vec4 specular)
{
	vec3 s = normalize(vec3(uLight[lightNum].position)); //directional light
	vec3 v = normalize(vec3(-pos));
	vec3 h = normalize(v + s);
	light += uLight[lightNum].diffuse * material.diffuse * max(dot(s, n), 0.0);
	specular += vec4(spec) * material.specular * uLight[lightNum].diffuse * pow(max(dot(h, n), 0.0), material.shininess);
	specular.a = 0.0;
	light.a = 1.0;
}
#endif

vec2 sincos(in float x)
{
	float s = sin(x);
	return vec2(s, sqrt(1.0 - (s * s)));
}

#define kPI 3.1415926536
vec3 decode(in vec2 enc)
{
	vec2 ang = (enc * vec2(2.0)) - vec2(1.0);
	vec2 scth = sincos(ang.x * kPI);
	vec2 scphi = vec2(sqrt(1.0 - ang.y*ang.y), ang.y);
	return normalize(vec3(scth.y*scphi.x, scth.x*scphi.x, scphi.y));
}

void main(void)
{
	vec4 color = material.diffuse;

	// in wNormal is the world-space normal of the fragment
	vec3 blending = abs( wNormal );
	blending = normalize(max(blending, 0.00001)); // Force weights to sum to 1.0
	float b = (blending.x + blending.y + blending.z);
	blending /= vec3(b, b, b);
	
	//diffuse + intensity
	vec4 tex0 = getTriplanarTex(blending, wCoords01, texture0);
	
	//normal(enc) + specular + AO
	vec4 tex1 = getTriplanarTex(blending, wCoords01, texture1);
	
	//diffuse + intensity
	vec4 tex2 = getTriplanarTex(blending, wCoords23, texture2);
	
	//normal(enc) + specular + AO
	vec4 tex3 = getTriplanarTex(blending, wCoords23, texture3);
	
	//mix(0.33, 0.66, texCoord0.x)
	//color *= vec4(blend(tex0, tex2, texCoord0.x), 1.0);
	
	color *= vec4(mix(tex0.xyz, tex2.xyz, texCoord0.x), 1.0);
	float spec = mix(tex1.z, tex3.z, texCoord0.x);
	float ambi = mix(tex1.w, tex3.w, texCoord0.x);

//directional lighting
#if (NUM_LIGHTS > 0)
	vec3 bump0 = (decode(tex1.xy) * 2.0) - vec3(1.0);
	vec3 bump1 = (decode(tex3.xy) * 2.0) - vec3(1.0);
	vec3 bump = mix(bump0, bump1, texCoord0.x);
	
	mat3 tangentFrame = mat3(tangent, bitangent, normal);
	vec3 v_normal = tangentFrame * bump;

	//v_normal = normal;
	//color = vec4(v_normal, 1.0);

	//ambient only make sense with lighting
	vec4 light = scene.ambient;
	vec4 specular = vec4(0.0);
	//color = uLight[0].diffuse;
	//ads(spec, 0, eyePos, v_normal, light, specular);
	#if (NUM_LIGHTS == 1)
		ads(spec, 0, eyePos, v_normal, light, specular);
	#else
		for (int i=0; i<NUM_LIGHTS; ++i) {
			ads(spec, i, eyePos, v_normal, light, specular);
		}
	#endif
	
	// this is crude "baked ambient occlusion" - basically multiply everything by the ambient texture
	// scaling whatever we've decided the lighting contribution is by 0.0 to 1.0 to account for sheltered/hidden surfaces
	light *= vec4(ambi, ambi, ambi, 1.0);

#endif //NUM_LIGHTS
	//light = vec4(1.0);
#if (NUM_LIGHTS > 0)
	frag_color = color; //* light;// + specular;
#else
	frag_color = color;
#endif
	SetFragDepth();
}
