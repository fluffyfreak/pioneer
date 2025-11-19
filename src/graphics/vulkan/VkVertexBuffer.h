// Copyright © 2008-2024 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#pragma once

#ifndef _VKVERTEXBUFFER_H
#define _VKVERTEXBUFFER_H

#include "graphics/Types.h"
#include "graphics/VertexBuffer.h"

#include <memory>

namespace Graphics {

	namespace Vulkan {

		class VertexBuffer : public Graphics::VertexBuffer {
		public:
			VertexBuffer(BufferUsage usage, uint32_t numVertices, uint32_t stride) :
				Graphics::VertexBuffer(usage, numVertices, stride),
				m_mapStart(0),
				m_mapLength(0)
			{
			}
			~VertexBuffer()
			{
				delete[] m_data;
			}

			void Unmap() final {}

			void UnmapRange(bool flush) final {}

			void FlushRange(size_t, size_t) final {}

			// change the buffer data without mapping
			void BufferData(const size_t, void *) final {}

			// release the underlying GPU memory and recreate
			void Reset() final {};

			void Bind() final {}
			void Release() final {}

		protected:
			Uint8 *MapInternal(BufferMapMode) final { return m_buffer.get(); }
			uint8_t *MapRangeInternal(BufferMapMode, size_t, size_t) override { return nullptr; }
			uint8_t *m_data;

			size_t m_mapStart;
			size_t m_mapLength;

		private:
			std::unique_ptr<Uint8[]> m_buffer;
		};

		class IndexBuffer : public Graphics::IndexBuffer {
		public:
			IndexBuffer(Uint32 size, BufferUsage bu, IndexBufferSize el) :
				Graphics::IndexBuffer(size, bu, el)
			{
				if (el == INDEX_BUFFER_32BIT)
					m_buffer.reset(new Uint32[size]);
				else
					m_buffer16.reset(new Uint16[size]);
			}

			Uint32 *Map(BufferMapMode) final { return m_buffer.get(); }
			Uint16 *Map16(BufferMapMode) final { return m_buffer16.get(); }
			void Unmap() final {}

			void BufferData(const size_t, void *) final {}

			void Bind() final {}
			void Release() final {}

		private:
			std::unique_ptr<Uint32[]> m_buffer;
			std::unique_ptr<Uint16[]> m_buffer16;
		};

		class MeshObject final : public Graphics::MeshObject {
		public:
			MeshObject(const Graphics::VertexFormatDesc &fmt, Graphics::VertexBuffer *vtx, Graphics::IndexBuffer *idx) :
				m_vtxBuffer(static_cast<Vulkan::VertexBuffer *>(vtx)),
				m_idxBuffer(static_cast<Vulkan::IndexBuffer *>(idx)),
				m_format(fmt)
			{
				assert(m_vtxBuffer.Valid());
			}
			~MeshObject() final {}

			Graphics::VertexBuffer *GetVertexBuffer() const final { return m_vtxBuffer.Get(); }
			Graphics::IndexBuffer *GetIndexBuffer() const final { return m_idxBuffer.Get(); }
			const Graphics::VertexFormatDesc &GetFormat() const override { return m_format; }

			void Bind() final {}
			void Release() final {}

		protected:
			RefCountedPtr<VertexBuffer> m_vtxBuffer;
			RefCountedPtr<IndexBuffer> m_idxBuffer;
			VertexFormatDesc m_format;
		};

	} // namespace Vulkan
} // namespace Graphics

#endif
