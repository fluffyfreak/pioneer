#include "LuaOrbital.h"
#include "LuaUtils.h"

/*
 * Class: Orbital
 *
 * Class representing a Culture style Orbital (mini-Ringworld). Inherits from <Body>.
 */

static bool promotion_test(DeleteEmitter *o)
{
	return dynamic_cast<Orbital*>(o);
}

template <> const char *LuaObject<Orbital>::s_type = "Orbital";

template <> void LuaObject<Orbital>::RegisterClass()
{
	const char *l_parent = "Body";

	LuaObjectBase::CreateClass(s_type, l_parent, NULL, NULL, NULL);
	LuaObjectBase::RegisterPromotion(l_parent, s_type, promotion_test);
}
