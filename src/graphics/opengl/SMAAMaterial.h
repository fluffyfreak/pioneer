// Copyright © 2008-2014 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _OGL_SMAA_MATERIAL_H_
#define _OGL_SMAA_MATERIAL_H_

/*
 * Performs a gaussian horizontal blur on a fullscreen quad.
 */
#include "libs.h"
#include "MaterialGL.h"
#include "Program.h"

namespace Graphics {
	namespace OGL {
		// edge detection material
		class SMAAEdgeMaterial : public Material {
		public:
			SMAAEdgeMaterial() {
				texture0 = nullptr;
			}

			Program *CreateProgram(const MaterialDescriptor &) {
				return new Program("postprocessing", "smaa_edge", "", false);
			}

			virtual void Apply() {
				m_program->Use();
				if(texture0) {
					m_program->texture0.Set(texture0, 0);
				}
			}

			virtual void Unapply() {
				m_program->Unuse();
			}
		};
		
		// blend material
		class SMAABlendMaterial : public Material {
		public:
			SMAABlendMaterial() {
				texture0 = nullptr;
			}

			Program *CreateProgram(const MaterialDescriptor &) {
				return new Program("postprocessing", "smaa_blend", "", false);
			}

			virtual void Apply() {
				m_program->Use();
				if(texture0) {
					m_program->texture0.Set(texture0, 0);
				}
			}

			virtual void Unapply() {
				m_program->Unuse();
			}
		};
		
		// neighbourhood material
		class SMAANeighbourhoodMaterial : public Material {
		public:
			SMAANeighbourhoodMaterial() {
				texture0 = nullptr;
			}

			Program *CreateProgram(const MaterialDescriptor &) {
				return new Program("postprocessing", "smaa_neighbourhood", "", false);
			}

			virtual void Apply() {
				m_program->Use();
				if(texture0) {
					m_program->texture0.Set(texture0, 0);
				}
			}

			virtual void Unapply() {
				m_program->Unuse();
			}
		};
	}
}

#endif