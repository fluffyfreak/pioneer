// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

uniform sampler2D texture0; //depth map
in vec2 texCoord0;

out vec4 frag_color;

void main(void)
{
	float depthValue = texture(texture0, texCoord0).r;
    frag_color = vec4(vec3(depthValue), 1.0);
}
