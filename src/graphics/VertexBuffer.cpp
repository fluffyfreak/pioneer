// Copyright © 2008-2012 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include <cassert>
#include "VertexBuffer.h"
#include "utils.h"

namespace Graphics {

GLvbo::GLvbo(const int nElements, const void *pVertexBuf, const void *pNormalBuf, const void *pUVBuf, GLenum nUsage)
{
	assert(nullptr!=pVertexBuf);
	assert(0<nElements);
	assert( nUsage>=0x88E0 && nUsage<=0x88EA );

	mVAO = UINT_MAX;
	mVertObjId = UINT_MAX;
	mNormObjId = UINT_MAX;
	mUVObjId = UINT_MAX;

	Init(nElements, pVertexBuf, pNormalBuf, pUVBuf, nUsage);
}
GLvbo::~GLvbo()
{
	Cleanup();
}

// private
void GLvbo::Init(const int nElements, const void *pVertexBuf, const void *pNormalBuf, const void *pUVBuf, GLenum nUsage)
{
	assert(nullptr!=pVertexBuf);
	assert(0<nElements);
	assert( nUsage>=0x88E0 && nUsage<=0x88EA );

	Cleanup();

	glGenVertexArrays(1, &mVAO);
	glBindVertexArray(mVAO);

	// the vertex buffer is NOT optional
	glGenBuffers(1, &mVertObjId);
	glBindBuffer(GL_ARRAY_BUFFER, mVertObjId);
	glBufferData(GL_ARRAY_BUFFER, nElements * sizeof(vector3f), pVertexBuf, nUsage);

	if(pNormalBuf) {
		glGenBuffers(1, &mNormObjId);
		glBindBuffer(GL_ARRAY_BUFFER, mNormObjId);
		glBufferData(GL_ARRAY_BUFFER, nElements * sizeof(vector3f), pNormalBuf, nUsage);
	}

	if(pUVBuf) {
		glGenBuffers(1, &mUVObjId);
		glBindBuffer(GL_ARRAY_BUFFER, mUVObjId);
		glBufferData(GL_ARRAY_BUFFER, nElements * sizeof(vector2f), pUVBuf, nUsage);
	}

	// LOG_GLERROR();
}

void GLvbo::Cleanup()
{
	if(mVertObjId != UINT_MAX)
	{
		const GLboolean glbIsBuffer = glIsBuffer(mVertObjId);
		if(glbIsBuffer==GL_TRUE) {
			glDeleteBuffers(1, &mVertObjId);
			mVertObjId = UINT_MAX;
		}
	}
	if(mNormObjId != UINT_MAX)
	{
		const GLboolean glbIsBuffer = glIsBuffer(mNormObjId);
		if(glbIsBuffer==GL_TRUE) {
			glDeleteBuffers(1, &mNormObjId);
			mNormObjId = UINT_MAX;
		}
	}
	if(mUVObjId != UINT_MAX)
	{
		const GLboolean glbIsBuffer = glIsBuffer(mUVObjId);
		if(glbIsBuffer==GL_TRUE) {
			glDeleteBuffers(1, &mUVObjId);
			mUVObjId = UINT_MAX;
		}
	}
	if(mVAO != UINT_MAX)
	{
		const GLboolean glbIsVA = glIsVertexArray(mVAO);
		if(glbIsVA==GL_TRUE) {
			glDeleteVertexArrays(1, &mVAO);
			mVAO = UINT_MAX;
		}
	}
}

void GLvbo::Bind() const
{
	glBindVertexArray(mVAO);

	GLuint attributeIndex = 0;

	if(mVertObjId < UINT_MAX)
	{
		// 1rst attribute buffer : vertices
		glEnableVertexAttribArray(attributeIndex);
		glBindBuffer(GL_ARRAY_BUFFER, mVertObjId);
		glVertexAttribPointer(
			attributeIndex,     // attribute
			3,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			nullptr             // array buffer offset
		);
		++attributeIndex;
	}
	if(mNormObjId < UINT_MAX)
	{
		// 2nd attribute buffer : normals
		glEnableVertexAttribArray(attributeIndex);
		glBindBuffer(GL_ARRAY_BUFFER, mNormObjId);
		glVertexAttribPointer(
			attributeIndex,                   // attribute
			3,                                // size
			GL_FLOAT,                         // type
			GL_FALSE,                         // normalized?
			0,                                // stride
			nullptr                           // array buffer offset
		);
		++attributeIndex;
	}
	if(mUVObjId < UINT_MAX)
	{
		// 3rd attribute buffer : UVs
		glEnableVertexAttribArray(attributeIndex);
		glBindBuffer(GL_ARRAY_BUFFER, mUVObjId);
		glVertexAttribPointer(
			attributeIndex,                   // attribute
			2,                                // size
			GL_FLOAT,                         // type
			GL_FALSE,                         // normalized?
			0,                                // stride
			nullptr                           // array buffer offset
		);
		++attributeIndex;
	}
}

void GLvbo::Release() const
{
	GLuint attributeIndex = 0;
	if(mVertObjId != UINT_MAX) {
		glDisableVertexAttribArray(attributeIndex);
		++attributeIndex;
	}
	if(mNormObjId != UINT_MAX) {
		glDisableVertexAttribArray(attributeIndex);
		++attributeIndex;
	}
	if(mUVObjId != UINT_MAX) {
		glDisableVertexAttribArray(attributeIndex);
		++attributeIndex;
	}

	glBindVertexArray(0);
}

};
