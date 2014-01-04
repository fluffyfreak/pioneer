// Copyright © 2008-2013 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _GL2_GEOSPHEREMATERIAL_H
#define _GL2_GEOSPHEREMATEIRAL_H
/*
 * Programs & Materials used by terrain
 */
#include "libs.h"
#include "GL2Material.h"
#include "Program.h"
#include "galaxy/StarSystem.h"

namespace Graphics {
	namespace GL2 {
		class HMDWarpProgram : public Program {
		public:
			HMDWarpProgram(const std::string &filename, const std::string &defines);

			Uniform LensCenter;
			Uniform ScreenCenter;
			Uniform Scale;
			Uniform ScaleIn;
			Uniform HmdWarpParam;
			Uniform ChromAbParam;
			Uniform Texm;

		protected:
			virtual void InitUniforms();
		};

		class HMDWarpMaterial : public Material {
			virtual Program *CreateProgram(const MaterialDescriptor &);
			virtual void Apply();
			virtual void Unapply();

		protected:
			void SetGSUniforms();
		};
	}
}
#endif
