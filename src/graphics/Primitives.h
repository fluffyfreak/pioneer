// Copyright © 2008-2012 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef __glprimitives_h__
#define __glprimitives_h__

#include <cstdint>
#include "VertexBuffer.h"

namespace Graphics {

class CGLquad
{
private:
	VertexBuffer mVBO;

public:
	CGLquad(const bool bNormals_, const bool bUVs_);
	~CGLquad();

	void Render() const;
};

};

#endif // __glprimitives_h__