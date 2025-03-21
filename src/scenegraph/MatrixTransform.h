// Copyright © 2008-2025 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#ifndef _MATRIXTRANSFORM_H
#define _MATRIXTRANSFORM_H

#include "Group.h"
#include "matrix4x4.h"

namespace Graphics {
	class Renderer;
}

namespace SceneGraph {

	/*
	* Applies a matrix transform to child nodes
	*
	* Note: transforms are not automatically serialized when saving to disk;
	* they are derived from the original model and animations.
	* If you have programmatically positioned a MatrixTransform, it is your
	* responsibility to ensure the new position is properly serialized.
	*/
	class MatrixTransform : public Group {
	public:
		MatrixTransform(Graphics::Renderer *r, const matrix4x4f &m);
		MatrixTransform(const MatrixTransform &, NodeCopyCache *cache = 0);

		Node *Clone(NodeCopyCache *cache = 0) override;
		const char *GetTypeName() const override { return "MatrixTransform"; }
		void Accept(NodeVisitor &v) override;

		void Save(NodeDatabase &) override;
		static MatrixTransform *Load(NodeDatabase &);

		void Render(const matrix4x4f &trans, const RenderData *rd) override;
		void RenderInstanced(const std::vector<matrix4x4f> &trans, const RenderData *rd) override;

		const matrix4x4f &GetTransform() const { return m_transform; }
		void SetTransform(const matrix4x4f &m) { m_transform = m; }

		matrix4x4f CalcGlobalTransform() const override;

	protected:
		virtual ~MatrixTransform() {}

	private:
		matrix4x4f m_transform;
	};
} // namespace SceneGraph
#endif
