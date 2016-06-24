#include "lua_api/l_dll.h"
#include "lua_api/l_internal.h"
#include "common/c_converter.h"
#include "common/c_content.h"
#include "server.h"


// sampletext()
int ModApiDll::l_sampletext(lua_State *L)
{
	NO_MAP_LOCK_REQUIRED;
	actionstream<<"sampletext"<<std::endl;
	return 0; /* number of results */
}
int ModApiDll::l_loadobject(lua_State *L)
{
	NO_MAP_LOCK_REQUIRED;
	const char * libname = luaL_checkstring(L, 1);
	actionstream<<"LoadObject: "<<libname<<std::endl;
	return 0; /* number of results */
}

void ModApiDll::Initialize(lua_State *L, int top)
{
	API_FCT(sampletext);
	API_FCT(loadobject);
}
