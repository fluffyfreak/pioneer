// Copyright © 2008-2012 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "GL2Debug.h"
#include "libs.h"
#include "OS.h"
#include <gl\GLU.h>

namespace Graphics {
	namespace GL2 {
		#pragma optimize( "", off )
		void CheckGLError()
		{
			GLenum err = glGetError();
			if( err ) {
				const char * errmsg = (const char *)gluErrorString(err);
#ifndef NDEBUG
				OS::Error("GL Error: %s\n-----------------\n\n", (errmsg!=NULL) ? errmsg : "UNKNOWN");
#else
				OS::Warning("GL Error: %s\n-----------------\n\n", (errmsg!=NULL) ? errmsg : "UNKNOWN");
#endif
			}
		}
	}
}
