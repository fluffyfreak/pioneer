varying vec3 varyingEyepos;
varying vec3 varyingNormal;

void main(void)
{
	//vec4 p = vec4(gl_Normal,1.0);
	//float waveTime = 1.0;
	//float waveWidth = 0.0001;
	//float waveHeight = 0.003;
	//float theta = atan(p.z/pow(p.x*p.x+p.y*p.y, 0.5));
	//float phi = atan(p.y/p.x);
	//float height = (sin(waveWidth * theta + waveTime) * cos(waveWidth * phi + waveTime) * waveHeight);
	//p = gl_Vertex * (height + 1.0);
	
	//gl_Position = logarithmicTransformParam(p);
	gl_Position = logarithmicTransform();
	varyingEyepos = vec3(gl_ModelViewMatrix * gl_Vertex);
	varyingNormal = gl_NormalMatrix * gl_Normal;
}
