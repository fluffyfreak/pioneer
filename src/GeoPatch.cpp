// Copyright © 2008-2012 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include <cassert>
#include "utils.h"
#include "GeoMesh.h"
#include "GeoPatch.h"
#include "GeoContext.h"
#include "Pi.h"

// constructor
GeoPatch::GeoPatch(const GeoPatchContext &context_, GeoMesh *pGeoSphere_, 
	const vector3f &v0_, const vector3f &v1_, const vector3f &v2_, const vector3f &v3_, 
	const uint32_t depth_, const GeoPatchID &ID_)
	: mContext(context_), mpGeoSphere(pGeoSphere_), mV0(v0_), mV1(v1_), mV2(v2_), mV3(v3_), 
	mClipCentroid((v0_+v1_+v2_+v3_) * 0.25f), mCentroid(mClipCentroid.Normalized()), mDepth(depth_), mClipRadius(0.0f), mRoughLength(0.0f), 
	mPatchID(ID_), mHeightmap(0), mVBO(nullptr), mHasSplitRequest(false), parent(nullptr)
{
	for (int i=0; i<NUM_KIDS; i++) {
		edgeFriend[i]	= nullptr;
		kids[i]			= nullptr;
	}

	mClipRadius = std::max(mClipRadius, (mV0-mClipCentroid).Length());
	mClipRadius = std::max(mClipRadius, (mV1-mClipCentroid).Length());
	mClipRadius = std::max(mClipRadius, (mV2-mClipCentroid).Length());
	mClipRadius = std::max(mClipRadius, (mV3-mClipCentroid).Length());

	const float fDepth = float(Clamp<uint32_t>(mDepth, 1, 5));
 	const float distMult = 5.0f / fDepth;
	mRoughLength = GEOPATCH_SUBDIVIDE_AT_CAMDIST / powf(2.0f, fDepth) * distMult;
}

// destructor
GeoPatch::~GeoPatch() 
{
	for (int i=0; i<4; i++) {
		if (edgeFriend[i]) {
			edgeFriend[i]->NotifyEdgeFriendDeleted(this);
		}
	}
	for (int i=0; i<NUM_KIDS; i++) {
		edgeFriend[i] = nullptr;
		if(kids[i]) {
			delete kids[i];
			kids[i] = nullptr;
		}
	}

	if(glIsTexture(mHeightmap)==GL_TRUE) {
		glDeleteTextures(1, &mHeightmap);
		//checkGLError();
	}
	if( mVBO ) {
		delete mVBO;
		mVBO = nullptr;
	}
}

// Generates full-detail vertices, and also non-edge normals and colors
void GeoPatch::GenerateMesh() {
	////////////////////////////////////////////////////////////////
	// Create the base mesh that the heightmap will modify
	vector3f *vts = mContext.vertexs();
	GLfloat *pUV = mContext.uvs();
	assert(nullptr!=vts);
	const float fracStep = mContext.meshLerpStep();
	for (uint32_t y=0; y<mContext.edgeLen(); y++) {
		for (uint32_t x=0; x<mContext.edgeLen(); x++) {
			const float xfrac = float(x) * fracStep;
			const float yfrac = float(y) * fracStep;
			assert(xfrac<=1.0 && yfrac<=1.0);
			const vector3f p = GetSpherePoint(xfrac, yfrac);
			*(vts++) = p;

			*(pUV++) = xfrac;
			*(pUV++) = yfrac;
		}
	}
	assert(vts == &mContext.vertexs()[mContext.NUM_MESH_VERTS()]);
	assert(pUV == &mContext.uvs()[mContext.NUM_MESH_VERTS()*2]);

	// now create the VBO
	assert(nullptr==mVBO);
	mVBO = new Graphics::VertexBuffer( mContext.NUM_MESH_VERTS(), &mContext.vertexs()[0], nullptr, mContext.uvs() );
}

void GeoPatch::ReceiveHeightmaps(const SSplitResult *psr)
{
	if (mDepth<psr->depth) {
		// this should work because each depth should have a common history
		const uint32_t kidIdx = psr->data[0].patchID.GetPatchIdx(mDepth+1);
		kids[kidIdx]->ReceiveHeightmaps(psr);
	} else {
		const int nD = mDepth+1;
		for (int i=0; i<4; i++)
		{
			kids[i] = new GeoPatch(mContext, mpGeoSphere, 
				psr->data[i].v0, psr->data[i].v1, psr->data[i].v2, psr->data[i].v3, 
				nD, mPatchID.NextPatchID(nD,i));
		}

		for (int i=0; i<4; i++)
		{
			kids[i]->ReceiveHeightmapTex(psr->data[i].texID);
		}

		// hm.. edges. Not right to pass this
		// edgeFriend...
		kids[0]->edgeFriend[0] = GetEdgeFriendForKid(0, 0);
		kids[0]->edgeFriend[1] = kids[1];
		kids[0]->edgeFriend[2] = kids[3];
		kids[0]->edgeFriend[3] = GetEdgeFriendForKid(0, 3);
		kids[1]->edgeFriend[0] = GetEdgeFriendForKid(1, 0);
		kids[1]->edgeFriend[1] = GetEdgeFriendForKid(1, 1);
		kids[1]->edgeFriend[2] = kids[2];
		kids[1]->edgeFriend[3] = kids[0];
		kids[2]->edgeFriend[0] = kids[1];
		kids[2]->edgeFriend[1] = GetEdgeFriendForKid(2, 1);
		kids[2]->edgeFriend[2] = GetEdgeFriendForKid(2, 2);
		kids[2]->edgeFriend[3] = kids[3];
		kids[3]->edgeFriend[0] = kids[0];
		kids[3]->edgeFriend[1] = kids[2];
		kids[3]->edgeFriend[2] = GetEdgeFriendForKid(3, 2);
		kids[3]->edgeFriend[3] = GetEdgeFriendForKid(3, 3);
		kids[0]->parent = kids[1]->parent = kids[2]->parent = kids[3]->parent = this;
		for (int i=0; i<4; i++) {
			kids[i]->GenerateMesh();
		}
		for (int i=0; i<4; i++) {
			if(edgeFriend[i]) {
				edgeFriend[i]->NotifyEdgeFriendSplit(this);
			}
		}
		mHasSplitRequest = false;
	}
}

void GeoPatch::ReceiveHeightmapTex(const GLuint tex)
{
	mHeightmap = tex;
}

void GeoPatch::Render(const vector3f &campos, const Graphics::Frustum &frustum)
{
	if (kids[0]) {
		for (int i=0; i<4; i++) {
			kids[i]->Render(campos, frustum);
		}
	} 
	else if(parent)
	{
		if (!frustum.TestPoint(vector3d(mClipCentroid), mClipRadius))
			return;

		const vector3f relpos = mClipCentroid - campos;
		glPushMatrix();
		glTranslated(relpos.x, relpos.y, relpos.z);

		Pi::statSceneTris += 2*(mContext.edgeLen()-1)*(mContext.edgeLen()-1);

		float patchColour[4] = {
			float(mDepth+1) * (1.0f/float(GEOPATCH_MAX_DEPTH)), 
			0.0f, 
			float(GEOPATCH_MAX_DEPTH-mDepth) * (1.0f/float(GEOPATCH_MAX_DEPTH)), 
			1.0f
		};
		//glm::vec4 patchColour(1.0f, 1.0f, 1.0f, 1.0f);
		glUniform4fv(mContext.patchColourID(), 1, &patchColour[0]);
		glUniform3fv(mContext.patchV0ID(), 1, &mV0[0]);
		glUniform3fv(mContext.patchV1ID(), 1, &mV1[0]);
		glUniform3fv(mContext.patchV2ID(), 1, &mV2[0]);
		glUniform3fv(mContext.patchV3ID(), 1, &mV3[0]);
		glUniform1i(mContext.patchTexHeightmapID(), 0); //Texture unit 0 is for base images.

		// bind the heightmap texture
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, mHeightmap);

		assert(nullptr != mVBO);
		(*mVBO).Bind();
			
		// Index buffer
		const GLuint iBufIndex = determineIndexbuffer();
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mContext.elementBuffers(iBufIndex));

		// Draw the triangles !
		glDrawElements(
			GL_TRIANGLES,									// mode
			mContext.elementBuffersIndexCount(iBufIndex),	// count
			GL_UNSIGNED_SHORT,								// type
			nullptr											// element array buffer offset
		);

		(*mVBO).Release();
	}
}

void GeoPatch::LODUpdate(const vector3f &campos) {
	// there should be no LODUpdate'ing when we have active split requests
	if(mHasSplitRequest)
		return;

	bool canSplit = true;
	bool canMerge = (kids[0]!=NULL);

	// always split at first level
	if (parent) {
		for (int i=0; i<4; i++) {
			if (!edgeFriend[i]) { 
				canSplit = false; 
				break; 
			} else if (edgeFriend[i]->mDepth < mDepth) {
				canSplit = false;
				break;
			}
		}
		const float centroidDist = (campos - mCentroid).Length();
		const bool errorSplit = (centroidDist < mRoughLength);
		if( !(canSplit && (mDepth < GEOPATCH_MAX_DEPTH) && errorSplit) ) 
		{
			canSplit = false;
		}
	}

	if (canSplit) {
		if (!kids[0]) {
			mHasSplitRequest = true;
			SSplitRequestDescription *desc = new SSplitRequestDescription(mV0, mV1, mV2, mV3, mCentroid, mDepth, mPatchID);
			mpGeoSphere->AddSplitRequest(desc);
		} else {
			for (int i=0; i<4; i++) {
				kids[i]->LODUpdate(campos);
			}
		}
	} else if (canMerge) {
		for (int i=0; i<4; i++) { 
			delete kids[i]; 
			kids[i] = nullptr; 
		}
	}
}
