// Copyright © 2008-2014 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

varying vec3 vertex;
varying vec2 uv;

void main(void)
{
	gl_Position = vec4(vertex = vec3(gl_ModelViewMatrix * gl_Vertex), 1.0);
	uv = gl_MultiTexCoord0.xy;
}
