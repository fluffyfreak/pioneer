// Copyright © 2008-2014 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#pragma once

#ifndef _GL2_GENGASGIANTCOLOURMATERIAL_H
#define _GL2_GENGASGIANTCOLOURMATERIAL_H
/*
 * Material(s) used to generate 
 *
 */
#include "GL2Material.h"
#include "Program.h"

namespace Graphics {

	namespace GL2 {
		class GenGasGiantColourProgram : public Program {
		public:
			GenGasGiantColourProgram(const MaterialDescriptor &);

			Uniform v0, v1, v2, v3;
			Uniform fracStep;

		protected:
			virtual void InitUniforms();
		};

		class GenGasGiantColourMaterial : public Material { //unlit
		public:
			virtual Program *CreateProgram(const MaterialDescriptor &);
			virtual void Apply();
			virtual void Unapply();
		};
	}
}

#endif
