uniform vec4 atmosColor;
// to keep distances sane we do a nearer, smaller scam. this is how many times
// smaller the geosphere has been made
uniform float geosphereScale;
uniform float geosphereScaledRadius;
uniform float geosphereAtmosTopRad;
uniform vec3 geosphereCenter;
uniform float geosphereAtmosFogDensity;
uniform float geosphereAtmosInvScaleHeight;

#ifdef ECLIPSE
uniform int shadows;
uniform ivec3 occultedLight;
uniform vec3 shadowCentreX;
uniform vec3 shadowCentreY;
uniform vec3 shadowCentreZ;
uniform vec3 srad;
uniform vec3 lrad;
uniform vec3 sdivlrad;
#endif // ECLIPSE

uniform Material material;
uniform Scene scene;

in vec3 varyingEyepos;
in vec3 varyingNormal;
in vec3 varyingTexCoord0;

uniform sampler2D texture0; //diffuse
uniform float time;

out vec4 frag_color;

// function taken from https://www.shadertoy.com/view/XslGRr
// Created by inigo quilez - iq/2013
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

float iqnoise( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);
	f = f*f*(3.0-2.0*f);
	vec2 uv = (p.xy+vec2(37.0,17.0)*p.z) + f.xy;
	vec2 rg = texture2D( texture0, (uv+ 0.5)/256.0, -100.0 ).yx;
	return -1.0+2.0*mix( rg.x, rg.y, f.z );
}

float fbm(vec3 position, int octaves, float frequency, float persistence) {
	float total = 0.0;
	float maxAmplitude = 0.0;
	float amplitude = 1.0;
	for (int i = 0; i < octaves; i++) {
		total += iqnoise(position * frequency) * amplitude;
		frequency *= 2.0;
		maxAmplitude += amplitude;
		amplitude *= persistence;
	}
	return total / maxAmplitude;
}

float absNoise(vec3 position, int octaves, float frequency, float persistence) {
	float total = 0.0;
	float maxAmplitude = 0.0;
	float amplitude = 1.0;
	for (int i = 0; i < octaves; i++) {
		total += abs(iqnoise(position * frequency)) * amplitude;
		frequency *= 2.0;
		maxAmplitude += amplitude;
		amplitude *= persistence;
	}
	return total / maxAmplitude;
}

float ridgedNoise(vec3 position, int octaves, float frequency, float persistence) {
	float total = 0.0;
	float maxAmplitude = 0.0;
	float amplitude = 1.0;
	for (int i = 0; i < octaves; i++) {
		total += ((1.0 - abs(iqnoise(position * frequency))) * 2.0 - 1.0) * amplitude;
		frequency *= 2.0;
		maxAmplitude += amplitude;
		amplitude *= persistence;
	}
	return total / maxAmplitude;
}

float squaredNoise(vec3 position, int octaves, float frequency, float persistence) {
	float total = 0.0;
	float maxAmplitude = 0.0;
	float amplitude = 1.0;
	for (int i = 0; i < octaves; i++) {
		float tmp = iqnoise(position * frequency);
		total += tmp * tmp * amplitude;
		frequency *= 2.0;
		maxAmplitude += amplitude;
		amplitude *= persistence;
	}
	return total / maxAmplitude;
}

float cubedNoise(vec3 position, int octaves, float frequency, float persistence) {
	float total = 0.0;
	float maxAmplitude = 0.0;
	float amplitude = 1.0;
	for (int i = 0; i < octaves; i++) {
		float tmp = iqnoise(position * frequency);
		total += tmp * tmp * tmp * amplitude;
		frequency *= 2.0;
		maxAmplitude += amplitude;
		amplitude *= persistence;
	}
	return total / maxAmplitude;
}


#ifdef ECLIPSE
#define PI 3.141592653589793

float discCovered(const in float dist, const in float rad) {
	// proportion of unit disc covered by a second disc of radius rad placed
	// dist from centre of first disc.
	//
	// XXX: same function is in Camera.cpp
	//
	// WLOG, the second disc is displaced horizontally to the right.
	// xl = rightwards distance to intersection of the two circles.
	// xs = normalised leftwards distance from centre of second disc to intersection.
	// d = vertical distance to an intersection point
	//
	// The clamps on xl,xs handle the cases where one disc contains the other.

	float radsq = rad*rad;

	float xl = clamp((dist*dist + 1.0 - radsq) / (2.0*max(0.001,dist)), -1.0, 1.0);
	float xs = clamp((dist - xl)/max(0.001,rad), -1.0, 1.0);
	float d = sqrt(max(0.0, 1.0 - xl*xl));

	float th = clamp(acos(xl), 0.0, PI);
	float th2 = clamp(acos(xs), 0.0, PI);

	// covered area can be calculated as the sum of segments from the two
	// discs plus/minus some triangles, and it works out as follows:
	return clamp((th + radsq*th2 - dist*d)/PI, 0.0, 1.0);
}
#endif // ECLIPSE

void main(void)
{
	vec3 eyepos = varyingEyepos;
	vec3 eyenorm = normalize(eyepos);
	vec3 tnorm = normalize(varyingNormal);
	
	// generate some noise clouds
	const float nScale = 1.4; // Uniform?
	const float Density = 0.02;
	
	
	vec3 noisePosition = varyingTexCoord0 + (vec3(0.0, 0.1, 1.0) * time);
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

	vec3 v = (eyepos - geosphereCenter)/geosphereScaledRadius;
	float lenInvSq = 1.0/(dot(v,v));
	for (int i=0; i<NUM_LIGHTS; ++i) {
		vec3 lightDir = normalize(vec3(uLight[i].position));
		float unshadowed = 1.0;
#ifdef ECLIPSE
		for (int j=0; j<shadows; j++) {
			if (i != occultedLight[j])
			continue;
			
			vec3 centre = vec3( shadowCentreX[j], shadowCentreY[j], shadowCentreZ[j] );
			
			// Apply eclipse:
			vec3 projectedPoint = v - dot(lightDir,v)*lightDir;
			// By our assumptions, the proportion of light blocked at this point by
			// this sphere is the proportion of the disc of radius lrad around
			// projectedPoint covered by the disc of radius srad around shadowCentre.
			float dist = length(projectedPoint - centre);
			unshadowed *= 1.0 - discCovered(dist/lrad[j], sdivlrad[j]);
		}
#endif // ECLIPSE
		unshadowed = clamp(unshadowed, 0.0, 1.0);
		nDotVP  = max(0.0, dot(tnorm, normalize(vec3(uLight[i].position))));
		nnDotVP = max(0.0, dot(tnorm, normalize(-vec3(uLight[i].position)))); //need backlight to increase horizon
		diff += uLight[i].diffuse * unshadowed * 0.5*(nDotVP+0.5*clamp(1.0-nnDotVP*4.0,0.0,1.0) * INV_NUM_LIGHTS);
	}

	// when does the eye ray intersect atmosphere
	float atmosStart = findSphereEyeRayEntryDistance(geosphereCenter, eyepos, geosphereScaledRadius * geosphereAtmosTopRad);
	float ldprod=0.0;
	float fogFactor=0.0;
	{
		float atmosDist = geosphereScale * (length(eyepos) - atmosStart);
		
		// a&b scaled so length of 1.0 means planet surface.
		vec3 a = (atmosStart * eyenorm - geosphereCenter) / geosphereScaledRadius;
		vec3 b = (eyepos - geosphereCenter) / geosphereScaledRadius;
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
