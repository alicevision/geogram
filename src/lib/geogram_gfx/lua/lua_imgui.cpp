/*
 *  Copyright (c) 2012-2016, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#include <geogram_gfx/lua/lua_imgui.h>
#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>
#include <geogram_gfx/ImGui_ext/imgui_ext.h>
#include <geogram/lua/lua_wrap.h>

extern void LoadImguiBindings();
extern lua_State* lState;

namespace {
    using namespace GEO;
    
    int wrapper_TextInput(lua_State* L) {
	
	if(
	    lua_gettop(L) != 2 &&
	    lua_gettop(L) != 3
	) {
	    return luaL_error(
		L, "'imgui.TextInput()' invalid number of arguments"
	    );
	}
	
	if(!lua_isstring(L,1)) {
	    return luaL_error(
		L, "'imgui.TextInput()' argument 1 should be a string"
	    );
	}
	
	if(!lua_isstring(L,2)) {
	    return luaL_error(
		L, "'imgui.TextInput()' argument 2 should be a string"
	    );
	}

	ImGuiInputTextFlags flags = 0;
	
	if(lua_gettop(L) == 3) {
	    if(!lua_isnumber(L,3)) {
		return luaL_error(
		    L, "'imgui.TextInput()' argument 3 should be a number"
		);
	    }
	    flags = ImGuiInputTextFlags(lua_tonumber(L,3));
	}
	
	const char* label  = lua_tostring(L,1);
	const char* str = lua_tostring(L,2);	
	static char buff[geo_imgui_string_length];
	strcpy(buff,str);
	bool result = ImGui::InputText(
	    label, buff, geo_imgui_string_length, flags
	);
	lua_pushboolean(L,result);
	lua_pushstring(L,buff);
	
	// Note: If I do not do that, then the system thinks that the
	// enter key is still pressed (it misses the key release event)
	// I do not know why...
	if(result) {
            const int key_index = ImGui::GetIO().KeyMap[ImGuiKey_Enter];
            ImGui::GetIO().KeysDown[key_index] = false;
	}
	
	return 2;
    }


    int wrapper_Combo(lua_State* L) {
	if(lua_gettop(L) != 3) {
	    return luaL_error(
		L, "'imgui.Combo()' invalid number of arguments"
	    );
	}
	
	if(!lua_isstring(L,1)) {
	    return luaL_error(
		L, "'imgui.Combo()' argument should be a string"
	    );
	}
	
	if(!lua_isstring(L,2)) {
	    return luaL_error(
		L, "'imgui.Combo()' argument should be a string"
	    );
	}

	if(!lua_isstring(L,3)) {
	    return luaL_error(
		L, "'imgui.Combo()' argument should be a string"
	    );
	}
	
	const char* label = lua_tostring(L,1);
	const char* current_item = lua_tostring(L,2);
	const char* items = lua_tostring(L,3);		

	char* lua_items = (char*)alloca(strlen(items)+2);
	strcpy(lua_items,items);
	{
	    size_t n = strlen(lua_items);
	    lua_items[n] = ';';
	    lua_items[n+1] = '\0';
	}
	
	int lua_current_item=0;
	
	const char* prev_item = lua_items;
	int nb_items = 0;

	char* p = lua_items;
	while(*p != '\0') {
	    if(*p == ';') {
		*p = '\0';
		if(!strcmp(prev_item, current_item)) {
		    lua_current_item = nb_items;
		}
		prev_item = p+1;
		++nb_items;
	    }
	    ++p;
	}
	*p = '\0'; // Double '\0' to indicate end of item list to lua.

	bool result = ImGui::Combo(label, &lua_current_item, lua_items);

	current_item = lua_items;
	while(lua_current_item > 0) {
	    while(*current_item) {
		++current_item;
	    }
	    ++current_item;
	    --lua_current_item;
	}

	lua_pushboolean(L, result);
	lua_pushstring(L, current_item);

	return 2;
    }

    int wrapper_ColorEdit3WithPalette(
	lua_State* L	
    ) {
	if(lua_gettop(L) != 4) {
	    return luaL_error(
		L, "'imgui.ColorEdit3WithPalette()' invalid number of arguments"
	    );
	}
	
	if(!lua_isstring(L,1)) {
	    return luaL_error(
		L, "'imgui.ColorEdit3WithPalette()' argument 1 should be a string"
	    );
	}

	if(!lua_isnumber(L,2)) {
	    return luaL_error(
		L, "'imgui.ColorEdit3WithPalette()' argument 2 should be a number"
	    );
	}

	if(!lua_isnumber(L,3)) {
	    return luaL_error(
		L, "'imgui.ColorEdit3WithPalette()' argument 3 should be a number"
	    );
	}

	if(!lua_isnumber(L,4)) {
	    return luaL_error(
		L, "'imgui.ColorEdit3WithPalette()' argument 4 should be a number"
	    );
	}

	const char* label = lua_tostring(L,1);
	
	float rgb[3];
	rgb[0] = float(lua_tonumber(L,2));
	rgb[1] = float(lua_tonumber(L,3));
	rgb[2] = float(lua_tonumber(L,4));	

	bool sel = ImGui::ColorEdit3WithPalette(
	    label, rgb
	);

	lua_pushboolean(L,sel);
	lua_pushnumber(L,rgb[0]);
	lua_pushnumber(L,rgb[1]);
	lua_pushnumber(L,rgb[2]);

	return 4;
    }
    
    int wrapper_OpenFileDialog(
	lua_State* L
    ) {
	if(lua_gettop(L) != 4) {
	    return luaL_error(
		L, "'imgui.OpenFileDialog()' invalid number of arguments"
	    );
	}

	if(!lua_isstring(L,1)) {
	    return luaL_error(
		L, "'imgui.OpenFileDialog()' argument 1 should be a string"
	    );
	}

	if(!lua_isstring(L,2)) {
	    return luaL_error(
		L, "'imgui.OpenFileDialog()' argument 2 should be a string"
	    );
	}

	if(!lua_isstring(L,3)) {
	    return luaL_error(
		L, "'imgui.OpenFileDialog()' argument 3 should be a string"
	    );
	}

	if(!lua_isnumber(L,4)) {
	    return luaL_error(
		L, "'imgui.OpenFileDialog()' argument 4 should be a number"
	    );
	}

	const char* label      = lua_tostring(L,1);
	const char* extensions = lua_tostring(L,2);
	const char* filename   = lua_tostring(L,3);
	ImGuiExtFileDialogFlags flags = ImGuiExtFileDialogFlags(lua_tonumber(L,4));

	ImGui::OpenFileDialog(label, extensions, filename, flags);
	
	return 0;
    }

    int wrapper_FileDialog(
	lua_State* L
    ) {
	if(lua_gettop(L) != 2) {
	    return luaL_error(
		L, "'imgui.FileDialog()' invalid number of arguments"
	    );
	}

	if(!lua_isstring(L,1)) {
	    return luaL_error(
		L, "'imgui.OpenFileDialog()' argument 1 should be a string"
	    );
	}

	if(!lua_isstring(L,2)) {
	    return luaL_error(
		L, "'imgui.OpenFileDialog()' argument 2 should be a string"
	    );
	}

	const char* label      = lua_tostring(L,1);
	char filename[geo_imgui_string_length];

	const char* filename_in = lua_tostring(L,3);
	if(filename_in != nil) {
	    if(strlen(filename_in) > geo_imgui_string_length + 1) {
		Logger::err("ImGui") << "Max file name length exceeded"
				     << std::endl;
		return false;
	    }
	    strcpy(filename, filename_in);
	} else {
	    filename[0] = '\0';
	}
	
	bool result = ImGui::FileDialog(label, filename, geo_imgui_string_length);

	lua_pushboolean(L,result);
	lua_pushstring(L, filename);
	
	return 2;
    }

    int wrapper_SetNextWindowPos(lua_State* L) {
	if(
	    lua_gettop(L) != 2 &&
	    lua_gettop(L) != 3
	) {
	    return luaL_error(
		L, "'imgui.SetNextWindowPos()' invalid number of arguments"
	    );
	}

	if(!lua_isnumber(L,1)) {
	    return luaL_error(
		L, "'imgui.SetNextWindowPos()' argument 1 should be a number"
	    );
	}

	if(!lua_isnumber(L,2)) {
	    return luaL_error(
		L, "'imgui.SetNextWindowPos()' argument 2 should be a number"
	    );
	}

	ImGuiCond cond = 0;
	if(lua_gettop(L) == 3) {
	    if(!lua_isnumber(L,3)) {
		return luaL_error(
		    L, "'imgui.SetNextWindowPos()' argument 3 should be a number"
		);
	    }
	    cond = ImGuiCond(lua_tonumber(L,3));
	}

	
	ImGui::SetNextWindowPos(
	    ImVec2(float(lua_tonumber(L,1)), float(lua_tonumber(L,2))),
	    cond
	);

	return 0;
    }

    int wrapper_SetNextWindowSize(lua_State* L) {
	if(
	    lua_gettop(L) != 2 &&
	    lua_gettop(L) != 3
	) {
	    return luaL_error(
		L, "'imgui.SetNextWindowSize()' invalid number of arguments"
	    );
	}

	if(!lua_isnumber(L,1)) {
	    return luaL_error(
		L, "'imgui.SetNextWindowSize()' argument 1 should be a number"
	    );
	}

	if(!lua_isnumber(L,2)) {
	    return luaL_error(
		L, "'imgui.SetNextWindowSize()' argument 2 should be a number"
	    );
	}

	ImGuiCond cond = 0;
	if(lua_gettop(L) == 3) {
	    if(!lua_isnumber(L,3)) {
		return luaL_error(
		    L, "'imgui.SetNextWindowSize()' argument 3 should be a number"
		);
	    }
	    cond = ImGuiCond(lua_tonumber(L,3));
	}

	
	ImGui::SetNextWindowSize(
	    ImVec2(float(lua_tonumber(L,1)), float(lua_tonumber(L,2))),
	    cond
	);

	return 0;
    }

    
}

void init_lua_imgui(lua_State* L) {
    lState = L;
    LoadImguiBindings();

    lua_pushinteger(L, ImGuiExtFileDialogFlags_Load);
    lua_setglobal(L,"ImGuiExtFileDialogFlags_Load");

    lua_pushinteger(L, ImGuiExtFileDialogFlags_Save);
    lua_setglobal(L,"ImGuiExtFileDialogFlags_Save");

    lua_pushinteger(L, ImGuiCond_Always);
    lua_setglobal(L,"ImGuiCond_Always");

    lua_pushinteger(L, ImGuiCond_Once);
    lua_setglobal(L,"ImGuiCond_Once");

    lua_pushinteger(L, ImGuiCond_FirstUseEver);
    lua_setglobal(L,"ImGuiCond_FirstUseEver");
    
    lua_pushinteger(L, ImGuiCond_Appearing);
    lua_setglobal(L,"ImGuiCond_Appearing");

    lua_getglobal(L, "imgui");

    lua_pushliteral(L,"TextInput");
    lua_pushcfunction(L,wrapper_TextInput); 
    lua_settable(L,-3);

    lua_pushliteral(L,"Combo");
    lua_pushcfunction(L,wrapper_Combo); 
    lua_settable(L,-3);

    lua_pushliteral(L,"ColorEdit3WithPalette");
    lua_pushcfunction(L,wrapper_ColorEdit3WithPalette); 
    lua_settable(L,-3);

    lua_pushliteral(L,"OpenFileDialog");
    lua_pushcfunction(L,wrapper_OpenFileDialog);
    lua_settable(L,-3);

    lua_pushliteral(L,"FileDialog");
    lua_pushcfunction(L,wrapper_FileDialog);
    lua_settable(L,-3);

    lua_pushliteral(L,"SetNextWindowPos");
    lua_pushcfunction(L,wrapper_SetNextWindowPos);
    lua_settable(L,-3);

    lua_pushliteral(L,"SetNextWindowSize");
    lua_pushcfunction(L,wrapper_SetNextWindowSize);
    lua_settable(L,-3);
    
    lua_pop(L,1);
}

