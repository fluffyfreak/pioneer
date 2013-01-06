// Copyright © 2008-2012 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "GeoPatchGenMaterial.h"
#include "GeoSphere.h"
#include "GeoContext.h"
#include "StringF.h"
#include "FileSystem.h"
#include "graphics/Graphics.h"
#include "graphics/TextureGL.h"
#include "graphics/RendererGL2.h"
#include "graphics/GL2/GL2Debug.h"
#include <sstream>

namespace Graphics {
namespace GL2 {

GeoPatchGenProgram::GeoPatchGenProgram(const std::string &filename, const std::string& vsdefines, const std::string& fsdefines)
{
	m_name = filename;
	m_defines = vsdefines;
	m_frag_defines = fsdefines;
	LoadShaders(filename, vsdefines, fsdefines);
	InitUniforms();
	Graphics::GL2::CheckGLError();
}

//virtual 
void GeoPatchGenProgram::InitUniforms()
{
	v0.Init("v0", m_program);
	v1.Init("v1", m_program);
	v2.Init("v2", m_program);
	v3.Init("v3", m_program);
	fracStep.Init("fracStep", m_program);

	maxHeight.Init("maxHeight", m_program);
	seaLevel.Init("seaLevel", m_program);
	fracnum.Init("fracnum", m_program);

	octaves.Init("octaves", m_program);
	amplitude.Init("amplitude", m_program);
	lacunarity.Init("lacunarity", m_program);
	frequency.Init("frequency", m_program);

	Graphics::GL2::CheckGLError();
}

Program *GeoPatchGenMaterial::CreateProgram(const MaterialDescriptor &desc)
{
	assert(false);
	return NULL;
}

#pragma optimize( "", off )
Program *GeoPatchGenMaterial::CreateProgram( const std::string &vertstr, const std::string &fragstr, const vecBindings &includePaths /*= s_nullBindings*/)
{
	std::string vs;
	std::string fs;
	vecBindings::const_iterator iter = includePaths.begin();
	while (iter!=includePaths.end())
	{
		const std::string libpath( "shaders/gl2/" + (*iter).first );
		RefCountedPtr<FileSystem::FileData> libCode = FileSystem::gameDataFiles.ReadFile(libpath);
		assert(libCode);
		const std::string lib = libCode.Get()->AsStringRange().ToString();
		switch ((*iter).second)
		{
		case eBothShaders:
			vs += lib;
			fs += lib;
			break;
		case eVertShader:
			vs += lib;
			break;
		case eFragShader:
			fs += lib;
			break;
		}
		// Next!
		++iter;
	}

	vs += stringf("#define NUM_LIGHTS 0\n");
	fs += stringf("#define NUM_LIGHTS 0\n");
	m_program = new Graphics::GL2::GeoPatchGenProgram(vertstr, vs, fs);
	Graphics::GL2::CheckGLError();
	return m_program;
}

void GeoPatchGenMaterial::Apply()
{
	//XXX replace with actual material parameter
	glMaterialfv (GL_FRONT, GL_EMISSION, &emissive[0]);

	SetGSUniforms();
}

void GeoPatchGenMaterial::SetGSUniforms()
{
	GeoPatchGenProgram *p = static_cast<GeoPatchGenProgram*>(m_program);
	const GeoPatchContext::PatchGenData& hg = *static_cast<GeoPatchContext::PatchGenData*>(this->specialParameter0);

	p->Use();

	p->v0.Set(hg.v0);
	p->v1.Set(hg.v1);
	p->v2.Set(hg.v2);
	p->v3.Set(hg.v3);
	p->fracStep.Set(hg.fracStep);

	p->maxHeight.Set(hg.maxHeight);
	p->seaLevel.Set(hg.seaLevel);
	p->fracnum.Set(hg.fracnum);

	p->octaves.Set(hg.octaves, 10);
	p->amplitude.Set(hg.amplitude, 10);
	p->lacunarity.Set(hg.lacunarity, 10);
	p->frequency.Set(hg.frequency, 10);

	if (hg.usesHeightmap) {
		p->texture0.Set(hg.pTex, 0);
	}
}

}
}
