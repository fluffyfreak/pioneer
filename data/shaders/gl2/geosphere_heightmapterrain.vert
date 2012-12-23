varying vec3 varyingEyepos;
varying vec3 varyingNormal;
varying vec4 vertexColor;

uniform sampler2D texture0; //heightmap

uniform float radius;
uniform vec3 v0;
uniform vec3 v1;
uniform vec3 v2;
uniform vec3 v3;
uniform float fracStep;

// in patch surface coords, [0,1]
// v[0] to v[3] are the corner vertices
vec3 GetSpherePoint(float x, float y) {
	return normalize(v0 + x*(1.0-y)*(v1-v0) + x*y*(v2-v0) + (1.0-x)*y*(v3-v0));
}

void main(void)
{
#if 1
	const float heightscale = 0.05;
	//uv = gl_MultiTexCoord0;//vertexUVs;

	// get the uv offsets
	float x = gl_MultiTexCoord0.x;
	float y = gl_MultiTexCoord0.y;
	float xm1 = x-fracStep;
	float xp1 = x+fracStep;
	float ym1 = y-fracStep;
	float yp1 = y+fracStep;

	// get the heights
	float hxm1y = (texture2D(texture0, vec2(xm1,y)).x * heightscale) + 1.0;
	float hxp1y = (texture2D(texture0, vec2(xp1,y)).x * heightscale) + 1.0;
	float hxym1 = (texture2D(texture0, vec2(x,ym1)).x * heightscale) + 1.0;
	float hxyp1 = (texture2D(texture0, vec2(x,yp1)).x * heightscale) + 1.0;

	// normal
	vec3 x1 = GetSpherePoint(xm1,y) * hxm1y;
	vec3 x2 = GetSpherePoint(xp1,y) * hxm1y;
	vec3 y1 = GetSpherePoint(x,ym1) * hxm1y;
	vec3 y2 = GetSpherePoint(x,yp1) * hxm1y;

	vec3 xNormal = normalize(cross((x2-x1),(y2-y1)));

	// Normal of the the vertex, in camera space
	// Only correct if ModelMatrix does not scale the model ! Use its inverse transpose if not.
	varyingNormal = (gl_ModelViewMatrix * vec4(xNormal,0)).xyz;

	float height = (texture2D(texture0, gl_MultiTexCoord0.st).x * heightscale) + 1.0; 
	vec4 p = gl_Vertex * height * radius;
	p.w = 1.0;
	gl_Position = logarithmicTransformParam( p );
	varyingEyepos = vec3(gl_ModelViewMatrix * p);
	vertexColor = gl_Color;
#else
	gl_Position = logarithmicTransform();
	vertexColor = gl_Color;
	varyingEyepos = vec3(gl_ModelViewMatrix * gl_Vertex);
	varyingNormal = gl_NormalMatrix * gl_Normal;
#endif
}
