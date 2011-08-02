uniform vec4 atmosColor;
// to keep distances sane we do a nearer, smaller scam. this is how many times
// smaller the georing has been made
uniform float georingScale;
uniform float georingAtmosTopRad;
uniform vec3 georingCenter;
uniform float georingAtmosFogDensity;

varying vec4 varyingEyepos;

void sphereEntryExitDist(out float near, out float far, in vec3 sphereCenter, in vec3 eyeTo, in float radius)
{
	vec3 v = -sphereCenter;
	vec3 dir = normalize(eyeTo);
	float b = -dot(v, dir);
	float det = (b * b) - dot(v, v) + (radius * radius);
	near = 0.0;
	far = 0.0;
	if (det > 0.0) {
		det = sqrt(det);
		float i1 = b - det;
		float i2 = b + det;
		if (i2 > 0.0) {
			near = max(i1, 0.0);
			far = i2;
		}
	}
}

void main(void)
{
	float skyNear, skyFar;
	vec3 eyepos = vec3(varyingEyepos);
	sphereEntryExitDist(skyNear, skyFar, georingCenter, eyepos, georingAtmosTopRad);
	float atmosDist = georingScale * (skyFar - skyNear);
	float ldprod;
	{
		vec3 dir = normalize(eyepos);
		vec3 a = (skyNear * dir - georingCenter) / georingAtmosTopRad;
		vec3 b = (skyFar * dir - georingCenter) / georingAtmosTopRad;
		ldprod = AtmosLengthDensityProduct(a, b, atmosColor.w*georingAtmosFogDensity, atmosDist);
	}
	float fogFactor = 1.0 / exp(ldprod);
	vec4 atmosDiffuse = vec4(0.0,0.0,0.0,1.0);
	{
		vec3 surfaceNorm = normalize(eyepos - georingCenter);
		for (int i=0; i<NUM_LIGHTS; ++i) {
			atmosDiffuse += gl_LightSource[i].diffuse * max(0.0, dot(surfaceNorm, normalize(vec3(gl_LightSource[i].position))));
		}
	}
	atmosDiffuse.a = 1.0;
	//float sun = max(0.0, dot(normalize(eyepos),normalize(vec3(gl_LightSource[0].position))));
	gl_FragColor = (1.0-fogFactor) * (atmosDiffuse*
		vec4(atmosColor.r, atmosColor.g, atmosColor.b, 1.0));

#ifdef ZHACK
	SetFragDepth(gl_TexCoord[6].z);
#endif
}
