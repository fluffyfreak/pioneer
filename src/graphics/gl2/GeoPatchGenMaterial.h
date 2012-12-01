// Copyright © 2008-2012 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _GL2_GEOPATCHGENMATERIAL_H
#define _GL2_GEOPATCHGENMATEIRAL_H
/*
 * Programs & Materials used by terrain
 */
#include "libs.h"
#include "GL2Material.h"
#include "Program.h"
#include "galaxy/StarSystem.h"

namespace Graphics {
	namespace GL2 {
		class GeoPatchGenProgram : public Program {
		public:
			GeoPatchGenProgram(const std::string &filename, const std::string& defines);

			Uniform v0;
			Uniform v1;
			Uniform v2;
			Uniform v3;
			Uniform fracStep;
			
			Uniform maxHeight;
			Uniform seaLevel;
			Uniform fracnum;
			
			Uniform octaves;
			Uniform amplitude;
			Uniform lacunarity;
			Uniform frequency;
			
			Uniform heightmap;
			bool usesHeightmap;
		protected:
			virtual void InitUniforms();
		};

		enum EShaderBinding {
			eBothShaders = 0,
			eVertShader,
			eFragShader
		};

		typedef std::pair<std::string,EShaderBinding> ShaderBindPair;
		typedef std::vector<ShaderBindPair> vecBindings;
		static const vecBindings s_nullBindings;

		class GeoPatchGenMaterial : public Material {
			virtual Program *CreateProgram(const MaterialDescriptor &);
			virtual void Apply();

		protected:
			void SetGSUniforms();

		public:
			Program *CreateProgram(const std::string &vertstr, const std::string &fragstr, const vecBindings &includePaths = s_nullBindings);
		};
	}
}
#endif
