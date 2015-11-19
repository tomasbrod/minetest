#ifndef L_DLL_H_
#define L_DLL_H_

#include "lua_api/l_base.h"

class ModApiDll : public ModApiBase {
private:
	static int l_sampletext(lua_State *L);
	static int l_loadobject(lua_State *L);
public:
	static void Initialize(lua_State *L, int top);
};

#endif /* L_DLL_H_ */
