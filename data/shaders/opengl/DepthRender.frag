// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

// Ouput data
layout(location = 0) out float fragmentdepth;

void main(void)
{
	fragmentdepth = gl_FragCoord.z;
}
