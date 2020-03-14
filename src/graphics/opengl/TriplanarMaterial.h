// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _OGL_TRIPLANARMATERIAL_H
#define _OGL_TRIPLANARMATERIAL_H
/*
 * Triplanar texture mapping material
 *
 */
#include "MaterialGL.h"
#include "Program.h"

namespace Graphics {

	namespace OGL {
		class TriplanarProgram : public Program {
		public:
			TriplanarProgram(const MaterialDescriptor &, int numLights = 0);
			Uniform texScale01;
			Uniform texScale23;
		protected:
			virtual void InitUniforms();
		};

		/*
		* This material uses up to four variations of the Program,
		* one for each directional light
		*/
		class TriplanarMaterial : public Material {
		public:
			TriplanarMaterial();
			virtual Program *CreateProgram(const MaterialDescriptor &);
			virtual void SetProgram(Program *p);
			virtual void Apply();
			virtual void Unapply();

		private:
			Program* m_programs[5];
			Uint32 m_curNumLights;
		};
	}
}

#endif
