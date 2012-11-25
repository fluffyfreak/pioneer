// Copyright © 2008-2012 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef __GEOCONTEXT_H__
#define __GEOCONTEXT_H__

#include "graphics/FrameBuffer.h"
#include "graphics/Primitives.h"

#define GEOPATCH_SUBDIVIDE_AT_CAMDIST	2.0f	//1.5f
#define GEOPATCH_MAX_DEPTH  10

class GeoPatchContext : public RefCounted
{
private:
	////////////////////////////////////////////////////////////////
	// class static private members

	////////////////////////////////////////////////////////////////
	// private members
	const uint32_t mEdgeLen;
	const uint32_t mHalfEdgeLen;

	Graphics::CGLquad	mQuad;
	//Graphics::VertexBuffer *mVBO;
	Graphics::FrameBuffer	mFBO;
	
	vector3f * mVertexs;
	vector3f * mNormals;
	GLfloat * mUVs;

	static const uint32_t NUM_INDEX_LISTS = 16;
	GLuint mElementBuffers[NUM_INDEX_LISTS];
	GLuint mElementBuffersTriangleCounts[NUM_INDEX_LISTS];

	inline int VBO_COUNT_LO_EDGE() const { return 3*(mHalfEdgeLen); }
	inline int VBO_COUNT_HI_EDGE() const { return 3*(mEdgeLen-1); }
	inline int VBO_COUNT_MID_IDX() const { return (4*3*(mEdgeLen-3)) + 2*(mEdgeLen-3)*(mEdgeLen-3)*3; }

	// private copy ctor, done to prevent copying/cloning of GeoPatchContext.
	GeoPatchContext(const GeoPatchContext &ref) : mEdgeLen(ref.mEdgeLen), mHalfEdgeLen(ref.mHalfEdgeLen), mQuad(false, true), mFBO(mEdgeLen,mEdgeLen) {};

public:
	////////////////////////////////////////////////////////////////
	// public members
	
	////////////////////////////////////////////////////////////////
	// public methods

	inline uint32_t NUM_MESH_VERTS() const { return mEdgeLen*mEdgeLen; }
	inline uint32_t NUM_INDICES() const { return NUM_MESH_VERTS()*2*3; }

	int getIndices(std::vector<unsigned short> &pl, const unsigned int edge_hi_flags,
		unsigned short *midIndices, unsigned short *loEdgeIndices[4], unsigned short *hiEdgeIndices[4]) const;

	// constructor
	GeoPatchContext(const uint32_t edgeLen);

	// destructor
	~GeoPatchContext();

	vector3f * vertexs()	const { return mVertexs; }
	vector3f * normals()	const { return mNormals; }
	GLfloat * uvs()			const { return mUVs; }

	inline GLuint elementBuffers(const GLuint iBufIndex) const { 
		assert(iBufIndex>=0 && iBufIndex<NUM_INDEX_LISTS); 
		return mElementBuffers[iBufIndex]; 
	}
	inline GLuint elementBuffersIndexCount(const GLuint iBufIndex) const { 
		assert(iBufIndex>=0 && iBufIndex<NUM_INDEX_LISTS); 
		return mElementBuffersTriangleCounts[iBufIndex]*3; 
	}

	inline uint32_t edgeLen() const { return mEdgeLen; }
	inline uint32_t halfEdgeLen() const { return mHalfEdgeLen; }
	
	inline float meshLerpStep() const		{ return 1.0f/float(mEdgeLen-1); }
	inline float textureLerpStep() const	{ return 1.0f/float(mEdgeLen-1); }

	inline uint32_t fboWidth() const		{ return mFBO.Width(); }

private:
	struct SHeightmapGen{
		GLuint prog;
		GLint v0;
		GLint v1;
		GLint v2;
		GLint v3;
		GLint fracStep;

		GLint maxHeight;
		GLint seaLevel;
		GLint fracnum;

		GLint octaves;
		GLint amplitude;
		GLint lacunarity;
		GLint frequency;

		bool usesHeightmap;
		GLint heightmap;
	};
	std::vector<SHeightmapGen> mHeightmapProgs;
	size_t mCurrentHeightmapProg;
public:
	const SHeightmapGen& GetHeightmapGenData()	const { return mHeightmapProgs[mCurrentHeightmapProg]; }
	void renderHeightmap(const vector3f &v0, const vector3f &v1, const vector3f &v2, const vector3f &v3, const uint32_t targetTex) const;

private:
	GLuint patch_prog;
	GLuint patch_MatrixID;
	GLuint patch_ViewMatrixID;
	GLuint patch_ModelMatrixID;
	GLuint patch_radius;
	GLuint patch_v0;
	GLuint patch_v1;
	GLuint patch_v2;
	GLuint patch_v3;
	GLuint patch_fracStep;
	GLuint patch_colour;
	GLuint patch_texHeightmap;
public:
	GLuint patchV0ID()				const { return patch_v0; }
	GLuint patchV1ID()				const { return patch_v1; }
	GLuint patchV2ID()				const { return patch_v2; }
	GLuint patchV3ID()				const { return patch_v3; }
	GLuint patchColourID()			const { return patch_colour; }
	GLuint patchTexHeightmapID()	const { return patch_texHeightmap; }

	void UsePatchShader(const matrix4x4f &ViewMatrix, const matrix4x4f &ModelMatrix, const matrix4x4f &MVP) const;
};

#endif //__GEOCONTEXT_H__