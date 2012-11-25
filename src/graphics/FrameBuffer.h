// Copyright © 2008-2012 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef __framebuffer_h__
#define __framebuffer_h__

#include <cstdint>

namespace Graphics {

class FrameBuffer
{
private:
	uint32_t mFBO;
	uint32_t mTexture;
	const uint32_t mWidth;
	const uint32_t mHeight;

public:
	FrameBuffer(const uint32_t width, const uint32_t height);
	~FrameBuffer();

	void Bind() const;
	void Release() const;

	inline uint32_t Width()	const { return mWidth; }
	inline uint32_t Height()	const { return mHeight; }

	//void GetData(float *data) const;
	void CopyTexture(const uint32_t target) const;
};

};

#endif // __framebuffer_h__