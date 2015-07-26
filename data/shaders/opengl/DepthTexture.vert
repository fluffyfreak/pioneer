// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

out vec2 texCoord0;

void main(void)
{
    gl_Position = uViewProjectionMatrix * a_vertex;
	texCoord0 = a_uv0.xy;
}
