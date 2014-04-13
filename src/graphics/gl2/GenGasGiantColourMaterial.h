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
		enum GasGiantQuality {
			GEN_JUPITER_TEXTURE = 0,
			GEN_SATURN_TEXTURE,
			GEN_SATURN2_TEXTURE,
			GEN_NEPTUNE_TEXTURE,
			GEN_NEPTUNE2_TEXTURE,
			GEN_URANUS_TEXTURE
		};

		class GenGasGiantColourProgram : public Program {
		public:
			GenGasGiantColourProgram(const MaterialDescriptor &);

			Uniform v0, v1, v2, v3;
			Uniform fracStep;

			Uniform time;

			Uniform octaves;
			Uniform lacunarity;
			Uniform frequency;

			Uniform ggdarkColor;
			Uniform gglightColor;
			Uniform entropy;
			Uniform planetEarthRadii;

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
