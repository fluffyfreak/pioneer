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
		class OceanSphereProgram : public Program {
		public:
			OceanSphereProgram(const std::string &filename, const std::string &defines);

			Uniform atmosColor;
			Uniform oceansphereAtmosFogDensity;
			Uniform oceansphereAtmosInvScaleHeight;
			Uniform oceansphereAtmosTopRad; // in planet radii
			Uniform oceansphereCenter;
			Uniform oceansphereScale;
			Uniform oceansphereScaledRadius; // (planet radius) / scale

		protected:
			virtual void InitUniforms();
		};

		class OceanSphereSurfaceMaterial : public Material {
			virtual Program *CreateProgram(const MaterialDescriptor &);
			virtual void Apply();

		protected:
			void SetGSUniforms();
		};
	}
}
#endif
