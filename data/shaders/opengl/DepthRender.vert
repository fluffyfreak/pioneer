// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

void main(void)
{
    gl_Position = uViewProjectionMatrix * a_vertex;
}
