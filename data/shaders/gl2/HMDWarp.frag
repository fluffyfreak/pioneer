#ifdef TEXTURE0
uniform sampler2D texture0; //diffuse
varying vec2 texCoord0;
varying vec2 oTexCoord;
#endif

uniform vec2 LensCenter;
uniform vec2 ScreenCenter;
uniform vec2 Scale;
uniform vec2 ScaleIn;
uniform vec4 HmdWarpParam;

uniform Scene scene;
uniform Material material;

vec2 HmdWarp(vec2 in01)
{
   vec2  theta = (in01 - LensCenter) * ScaleIn; // Scales to [-1, 1]
   float rSq = theta.x * theta.x + theta.y * theta.y;
   vec2  theta1 = theta * (HmdWarpParam.x + HmdWarpParam.y * rSq + 
                           HmdWarpParam.z * rSq * rSq + HmdWarpParam.w * rSq * rSq * rSq);
   return LensCenter + Scale * theta1;
}

void main(void)
{
	vec4 color = material.diffuse;
	vec2 tc = HmdWarp(oTexCoord);
	if (!all(equal(clamp(tc, ScreenCenter-vec2(0.25,0.5), ScreenCenter+vec2(0.25,0.5)), tc)))
		color = vec4(0.0,0.0,0.0,1.0);//vec4(0.39,0.58,0.93,1.0);//vec4(0.1,0.0,0.1,1.0);
	else
		color = texture2D(texture0, tc);
	
	gl_FragColor = color;

	SetFragDepth();
}
