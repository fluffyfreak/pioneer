uniform mat4 Texm;
varying vec2 texCoord0;
varying vec2 oTexCoord;

void main(void)
{
	gl_Position = logarithmicTransform();
	texCoord0 = gl_MultiTexCoord0.xy;
	
	oTexCoord = vec2(Texm * vec4(texCoord0,0.0,1.0));
	//oTexCoord.y = 1.0-oTexCoord.y;
}