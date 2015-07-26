// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _GL_DEPTHTEXTUREMATERIAL_H
#define _GL_DEPTHTEXTUREMATERIAL_H
/*
 * A generic material & program for simple uses
 * textured/untextured, vertex colors or no...
 *
 */
#include "MaterialGL.h"
#include "Program.h"

namespace Graphics {

	namespace OGL {
		class DepthTextureProgram : public Program {
		public:
			DepthTextureProgram(const MaterialDescriptor &);
		};

		class DepthTextureMaterial : public Material {
		public:
			virtual Program *CreateProgram(const MaterialDescriptor &);
			virtual void Apply();
			virtual void Unapply();
		};
	}
}

#endif
