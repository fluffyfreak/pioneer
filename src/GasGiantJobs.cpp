// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "libs.h"
#include "GasGiant.h"
#include "GasGiantJobs.h"
#include "perlin.h"
#include "Pi.h"
#include "Game.h"
#include "RefCounted.h"
#include "graphics/Material.h"
#include "graphics/opengl/GenGasGiantColourMaterial.h"
#include "graphics/Renderer.h"
#include "graphics/Frustum.h"
#include "graphics/Graphics.h"
#include "graphics/Texture.h"
#include "graphics/TextureBuilder.h"
#include "graphics/VertexArray.h"
#include "graphics/Types.h"
#include "vcacheopt/vcacheopt.h"
#include <deque>
#include <algorithm>

namespace GasGiantJobs
{

	STextureFaceRequest::STextureFaceRequest(const vector3d *v_, const SystemPath &sysPath_, const Sint32 face_, const Sint32 uvDIMs_, Terrain *pTerrain_) :
		corners(v_), sysPath(sysPath_), face(face_), uvDIMs(uvDIMs_), pTerrain(pTerrain_)
	{
		colors = new Color[NumTexels()];
	}

	// RUNS IN ANOTHER THREAD!! MUST BE THREAD SAFE!
	// Use only data local to this object
	void STextureFaceRequest::OnRun()
	{
		PROFILE_SCOPED()
		MsgTimer timey;

		assert( corners != nullptr );
		double fracStep = 1.0 / double(UVDims()-1);
		for( Sint32 v=0; v<UVDims(); v++ ) {
			for( Sint32 u=0; u<UVDims(); u++ ) {
				// where in this row & colum are we now.
				const double ustep = double(u) * fracStep;
				const double vstep = double(v) * fracStep;

				// get point on the surface of the sphere
				const vector3d p = GetSpherePoint(ustep, vstep);
				// get colour using `p`
				const vector3d colour = pTerrain->GetColor(p, 0.0, p);

				// convert to ubyte and store
				Color* col = colors + (u + (v * UVDims()));
				col[0].r = Uint8(colour.x * 255.0);
				col[0].g = Uint8(colour.y * 255.0);
				col[0].b = Uint8(colour.z * 255.0);
				col[0].a = 255;
			}
		}

		timey.Mark("SingleTextureFaceCPUJob::OnRun");
	}

	// ********************************************************************************
	// Overloaded PureJob class to handle generating the mesh for each patch
	// ********************************************************************************
	SingleTextureFaceJob::~SingleTextureFaceJob()
	{
		PROFILE_SCOPED()
		if(mpResults) {
			mpResults->OnCancel();
			delete mpResults;
			mpResults = nullptr;
		}
	}

	void SingleTextureFaceJob::OnRun() // RUNS IN ANOTHER THREAD!! MUST BE THREAD SAFE!
	{
		PROFILE_SCOPED()
		mData->OnRun();

		// add this patches data
		STextureFaceResult *sr = new STextureFaceResult(mData->Face());
		sr->addResult(mData->Colors(), mData->UVDims());

		// store the result
		mpResults = sr;
	}

	void SingleTextureFaceJob::OnFinish() // runs in primary thread of the context
	{
		PROFILE_SCOPED()
		GasGiant::OnAddTextureFaceResult( mData->SysPath(), mpResults );
		mpResults = nullptr;
	}

	// ********************************************************************************
	SGPUGenRequest::SGPUGenRequest(const SystemPath &sysPath_, const Sint32 uvDIMs_, Terrain *pTerrain_, const float planetRadius_, GenFaceQuad* pQuad_, Graphics::Texture *pTex_) :
		m_texture(pTex_), sysPath(sysPath_), uvDIMs(uvDIMs_), pTerrain(pTerrain_), planetRadius(planetRadius_), pQuad(pQuad_)
	{
		PROFILE_SCOPED()
		assert(m_texture.Valid());
	}

	// ********************************************************************************
	void SGPUGenResult::addResult(Graphics::Texture *t_, Sint32 uvDims_) {
		PROFILE_SCOPED()
		mData = SGPUGenData(t_, uvDims_);
	}

	void SGPUGenResult::OnCancel()	{
		if( mData.texture ) { mData.texture.Reset(); }
	}

	//static matrix3x3f LookAt (const vector3f &eye, const vector3f &target, const vector3f &up) {
	//	const vector3f viewDir  = (target - eye).Normalized();
	//	const vector3f viewUp   = (up - up.Dot(viewDir) * viewDir).Normalized();
	//	const vector3f viewSide = viewDir.Cross(viewUp);
	//
	//	return matrix3x3f::FromVectors(viewSide, viewUp, -viewDir);
	//}

	// ********************************************************************************
	// Overloaded JobGPU class to handle generating the mesh for each patch
	// ********************************************************************************
	SingleGPUGenJob::SingleGPUGenJob(SGPUGenRequest *data) : mData(data), mpResults(nullptr) { /* empty */ }
	SingleGPUGenJob::~SingleGPUGenJob()
	{
		PROFILE_SCOPED()
		if(mpResults) {
			mpResults->OnCancel();
			delete mpResults;
			mpResults = nullptr;
		}
	}

	void SingleGPUGenJob::OnRun() // Runs in the main thread, may trash the GPU state
	{
		PROFILE_SCOPED()
		MsgTimer timey;

		for( Uint32 iFace=0; iFace<NUM_PATCHES; iFace++ )
		{
			// render the scene
			GasGiant::SetRenderTargetCubemap( iFace, mData->Texture() );
			GasGiant::BeginRenderTarget();
			Pi::renderer->BeginFrame(); 
			Pi::renderer->SetViewport(0, 0, mData->UVDims(), mData->UVDims());
			Pi::renderer->SetTransform(matrix4x4f::Identity());

			// enter ortho
			{
				Pi::renderer->SetMatrixMode(Graphics::MatrixMode::PROJECTION);
				Pi::renderer->PushMatrix();
				Pi::renderer->SetOrthographicProjection(0, mData->UVDims(), mData->UVDims(), 0, -1, 1);
				Pi::renderer->SetMatrixMode(Graphics::MatrixMode::MODELVIEW);
				Pi::renderer->PushMatrix();
				Pi::renderer->LoadIdentity();
			}

			// draw to the texture here
			mData->SetupMaterialParams( iFace );
			mData->Quad()->Draw( Pi::renderer );

			// leave ortho?
			{
				Pi::renderer->SetMatrixMode(Graphics::MatrixMode::PROJECTION);
				Pi::renderer->PopMatrix();
				Pi::renderer->SetMatrixMode(Graphics::MatrixMode::MODELVIEW);
				Pi::renderer->PopMatrix();
			}
			
			Pi::renderer->EndFrame();
			GasGiant::EndRenderTarget();
			GasGiant::SetRenderTargetCubemap( iFace, nullptr );
		}

		// add this patches data
		SGPUGenResult *sr = new SGPUGenResult();
		sr->addResult( mData->Texture(), mData->UVDims() );

		// store the result
		mpResults = sr;

		timey.Mark("SingleGPUGenJob::OnRun");
	}

	void SingleGPUGenJob::OnFinish() // runs in primary thread of the context
	{
		PROFILE_SCOPED()
		GasGiant::OnAddGPUGenResult( mData->SysPath(), mpResults );
		mpResults = nullptr;
	}
};

