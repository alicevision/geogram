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

#include <geogram_gfx/glup_viewer/glup_viewer_lua.h>
#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>
#include <geogram_gfx/glup_viewer/glup_viewer.h>
#include <geogram_gfx/GLUP/GLUP.h>
#include <geogram/basic/geometry.h>
#include <geogram/basic/stopwatch.h>

extern "C" {
#include <geogram/third_party/lua/lauxlib.h>
#include <geogram/third_party/lua/lualib.h>
}

#include <map>
#include <string>

extern void LoadImguiBindings();
extern lua_State* lState;

namespace {
    const GEO::index_t lua_glup_max_depth = 13;
    GEO::index_t lua_glup_modelview_depth = 0;
    GEO::index_t lua_glup_projection_depth = 0;
    GEO::index_t lua_glup_texture_depth = 0;
    bool lua_glup_primitive_active = false;
    std::map<std::string, const char*> embedded_files;
}

#define DECLARE_GLUP_CONSTANT(C) \
    lua_pushliteral(L,#C);       \
    lua_pushinteger(L,GLUP_##C); \
    lua_settable(L,1)

#define DECLARE_GLUP_FUNC(F)           \
    lua_pushliteral(L,#F);             \
    lua_pushcfunction(L,lua_glup_##F); \
    lua_settable(L,1)

static std::map<std::string, GEO::vec4> lua_glup_colormap;

static void DECLARE_GLUP_COLOR(
    const char* name, double r, double g, double b, double a=1.0
) {
    lua_glup_colormap[name] = GEO::vec4(r,g,b,a);
}

inline bool get_vec4(lua_State* L, double* xyzw, int pos=1) {
    int nargs = lua_gettop(L);
    xyzw[0] = 0.0;
    xyzw[1] = 0.0;
    xyzw[2] = 0.0;
    xyzw[3] = 1.0;

    if(nargs == pos && lua_isstring(L,pos)) {
	const char* name = lua_tostring(L,pos);
	std::map<std::string, GEO::vec4>::iterator it =
	    lua_glup_colormap.find(name);
	if(it == lua_glup_colormap.end()) {
	    return false;
	}
	xyzw[0] = it->second.x;
	xyzw[1] = it->second.y;
	xyzw[2] = it->second.z;
	xyzw[3] = it->second.w;	
    } else {
	if(nargs > pos+3 || nargs < pos) {
	    return false;
	}
	if(nargs == pos+3) {
	    if(!lua_isnumber(L,pos+3)) {
		return false;
	    }
	    xyzw[3] = lua_tonumber(L,pos+3);
	}
	if(nargs >= pos+2) {
	    if(!lua_isnumber(L,pos+2)) {
		return false;
	    }
	    xyzw[2] = lua_tonumber(L,pos+2);
	}
	if(nargs >= pos+1) {
	    if(!lua_isnumber(L,pos+1)) {
		return false;
	    }
	    xyzw[1] = lua_tonumber(L,pos+1);
	}
	if(nargs >= pos) {
	    if(!lua_isnumber(L,pos)) {
		return false;
	    }
	    xyzw[0] = lua_tonumber(L,pos);
	}
    }
    return true;
}

static int lua_glup_Enable(lua_State* L) {
    if(lua_gettop(L) != 1) {
	return luaL_error(
	    L, "'GLUP.Enable()' invalid number of arguments"
	);
    }
    if(!lua_isinteger(L,1)) {
	return luaL_error(
	    L, "'GLUP.Enable()' argument should be an integer"
	);
    }
    lua_Integer toggle = lua_tointeger(L,1);
    glupEnable(GLUPtoggle(toggle));
    return 0;
}

static int lua_glup_Disable(lua_State* L) {
    if(lua_gettop(L) != 1) {
	return luaL_error(
	    L, "'GLUP.Disable()' invalid number of arguments"
	);
    }
    if(!lua_isinteger(L,1)) {
	return luaL_error(
	    L, "'GLUP.Disable()' argument should be an integer"
	);
    }
    lua_Integer toggle = lua_tointeger(L,1);
    glupDisable(GLUPtoggle(toggle));
    return 0;
}

static int lua_glup_IsEnabled(lua_State* L) {
    if(lua_gettop(L) != 1) {
	return luaL_error(
	    L, "'GLUP.IsEnabled()' invalid number of arguments"
	);
    }
    if(!lua_isinteger(L,1)) {
	return luaL_error(
	    L, "'GLUP.IsEnabled()' argument should be an integer"
	);
    }
    lua_Integer toggle = lua_tointeger(L,1);
    lua_pushboolean(L,glupIsEnabled(GLUPtoggle(toggle)));
    return 1;
}

static int lua_glup_SetColor(lua_State* L) {
    if(lua_gettop(L) < 2) {
	return luaL_error(
	    L, "'GLUP.SetColor()' invalid number of arguments"
	);
    }
    if(!lua_isinteger(L,1)) {
	return luaL_error(
	    L, "'GLUP.SetColor()' argument should be an integer"
	);
    }
    lua_Integer color = lua_tointeger(L,1);
    double rgba[4];
    if(!get_vec4(L,rgba,2)) {
	return luaL_error(
	    L, "'GLUP.SetColor()' invalid arguments"
	);
    }
    glupSetColor4dv(GLUPcolor(color),rgba);
    return 0;

}

static int lua_glup_GetColor(lua_State* L) {
    if(lua_gettop(L) != 1) {
	return luaL_error(
	    L, "'GLUP.GetColor()' invalid number of arguments"
	);
    }
    if(!lua_isinteger(L,1)) {
	return luaL_error(
	    L, "'GLUP.GetColor()' argument should be an integer"
	);
    }
    lua_Integer color = lua_tointeger(L,1);
    float rgba[4];
    glupGetColor4fv(GLUPcolor(color),rgba);
    lua_pushnumber(L,lua_Number(rgba[0]));
    lua_pushnumber(L,lua_Number(rgba[1]));
    lua_pushnumber(L,lua_Number(rgba[2]));
    lua_pushnumber(L,lua_Number(rgba[3]));    
    return 4;
}

static int lua_glup_LightVector(lua_State* L) {
    if(lua_gettop(L) != 3) {
	return luaL_error(
	    L, "'GLUP.LightVector()' invalid number of arguments"
	);
    }
    if(
	!lua_isnumber(L,1) ||
	!lua_isnumber(L,2) ||
	!lua_isnumber(L,3)
    ) {
	return luaL_error(
	    L, "'GLUP.LightVector()' invalid arguments"
	);
    }
    glupLightVector3f(
	float(lua_tonumber(L,1)),
	float(lua_tonumber(L,2)),
	float(lua_tonumber(L,3))
    );
    return 0;
}

static int lua_glup_SetPointSize(lua_State* L) {
    if(lua_gettop(L) != 1) {
	return luaL_error(
	    L, "'GLUP.SetPointSize()' invalid number of arguments"
	);
    }
    if(!lua_isnumber(L,1)) {
	return luaL_error(
	    L, "'GLUP.SetPointSize()' invalid arguments"
	);
    }
    glupSetPointSize(float(lua_tonumber(L,1)));
    return 0;
}

static int lua_glup_GetPointSize(lua_State* L) {
    if(lua_gettop(L) != 0) {
	return luaL_error(
	    L, "'GLUP.GetPointSize()' invalid number of arguments"
	);
    }
    lua_pushnumber(L,lua_Number(glupGetPointSize()));
    return 1;
}

static int lua_glup_SetMeshWidth(lua_State* L) {
    if(lua_gettop(L) != 1) {
	return luaL_error(
	    L, "'GLUP.SetMeshWidth()' invalid number of arguments"
	);
    }
    if(!lua_isinteger(L,1)) {
	return luaL_error(
	    L, "'GLUP.SetMeshWidth()' invalid arguments"
	);
    }
    glupSetMeshWidth(int(lua_tointeger(L,1)));
    return 0;
}

static int lua_glup_GetMeshWidth(lua_State* L) {
    if(lua_gettop(L) != 0) {
	return luaL_error(
	    L, "'GLUP.GetMeshWidth()' invalid number of arguments"
	);
    }
    lua_pushinteger(L,glupGetMeshWidth());
    return 1;
}


static int lua_glup_SetCellsShrink(lua_State* L) {
    if(lua_gettop(L) != 1) {
	return luaL_error(
	    L, "'GLUP.SetCellsShrink()' invalid number of arguments"
	);
    }
    if(!lua_isnumber(L,1)) {
	return luaL_error(
	    L, "'GLUP.SetCellsShrink()' invalid arguments"
	);
    }
    glupSetCellsShrink(float(lua_tonumber(L,1)));
    return 0;
}

static int lua_glup_GetCellsShrink(lua_State* L) {
    if(lua_gettop(L) != 0) {
	return luaL_error(
	    L, "'GLUP.GetCellsShrink()' invalid number of arguments"
	);
    }
    lua_pushnumber(L,lua_Number(glupGetCellsShrink()));
    return 1;
}


static int lua_glup_SetAlphaThreshold(lua_State* L) {
    if(lua_gettop(L) != 1) {
	return luaL_error(
	    L, "'GLUP.SetAlphaThreshold()' invalid number of arguments"
	);
    }
    if(!lua_isnumber(L,1)) {
	return luaL_error(
	    L, "'GLUP.SetAlphaThreshold()' invalid arguments"
	);
    }
    glupSetAlphaThreshold(float(lua_tonumber(L,1)));
    return 0;
}

static int lua_glup_GetAlphaThreshold(lua_State* L) {
    if(lua_gettop(L) != 0) {
	return luaL_error(
	    L, "'GLUP.GetAlphaThreshold()' invalid number of arguments"
	);
    }
    lua_pushnumber(L,lua_Number(glupGetAlphaThreshold()));
    return 1;
}

static int lua_glup_PickingMode(lua_State* L) {
    if(lua_gettop(L) != 1) {
	return luaL_error(
	    L, "'GLUP.PickingMode()' invalid number of arguments"
	);
    }
    if(!lua_isinteger(L,1)) {
	return luaL_error(
	    L, "'GLUP.PickingMode()' invalid arguments"
	);
    }
    glupPickingMode(GLUPpickingMode(lua_tointeger(L,1)));
    return 0;
}

static int lua_glup_GetPickingMode(lua_State* L) {
    if(lua_gettop(L) != 0) {
	return luaL_error(
	    L, "'GLUP.SetPickingMode()' invalid number of arguments"
	);
    }
    lua_pushinteger(L,glupGetPickingMode());
    return 1;
}

static int lua_glup_PickingId(lua_State* L) {
    if(lua_gettop(L) != 1) {
	return luaL_error(
	    L, "'GLUP.PickingId()' invalid number of arguments"
	);
    }
    if(!lua_isinteger(L,1)) {
	return luaL_error(
	    L, "'GLUP.PickingId()' invalid arguments"
	);
    }
    glupPickingId(GLUPuint64(lua_tointeger(L,1)));
    return 0;
}

static int lua_glup_GetPickingId(lua_State* L) {
    if(lua_gettop(L) != 0) {
	return luaL_error(
	    L, "'GLUP.GetPickingId()' invalid number of arguments"
	);
    }
    lua_pushinteger(L,lua_Integer(glupGetPickingId()));
    return 1;
}

static int lua_glup_BasePickingId(lua_State* L) {
    if(lua_gettop(L) != 1) {
	return luaL_error(
	    L, "'GLUP.BasePickingId()' invalid number of arguments"
	);
    }
    if(!lua_isinteger(L,1)) {
	return luaL_error(
	    L, "'GLUP.BasePickingId()' invalid arguments"
	);
    }
    glupBasePickingId(GLUPuint64(lua_tointeger(L,1)));
    return 0;
}

static int lua_glup_GetBasePickingId(lua_State* L) {
    if(lua_gettop(L) != 0) {
	return luaL_error(
	    L, "'GLUP.GetBasePickingId()' invalid number of arguments"
	);
    }
    lua_pushinteger(L,lua_Integer(glupGetBasePickingId()));
    return 1;
}


static int lua_glup_ClipMode(lua_State* L) {
    if(lua_gettop(L) != 1) {
	return luaL_error(
	    L, "'GLUP.ClipMode()' invalid number of arguments"
	);
    }
    if(!lua_isinteger(L,1)) {
	return luaL_error(
	    L, "'GLUP.ClipMode()' invalid arguments"
	);
    }
    glupClipMode(GLUPclipMode(lua_tointeger(L,1)));
    return 0;
}

static int lua_glup_GetClipMode(lua_State* L) {
    if(lua_gettop(L) != 0) {
	return luaL_error(
	    L, "'GLUP.GetClipMode()' invalid number of arguments"
	);
    }
    lua_pushinteger(L,lua_Integer(glupGetClipMode()));
    return 1;
}

static int lua_glup_ClipPlane(lua_State* L) {
    if(lua_gettop(L) != 4) {
	return luaL_error(
	    L, "'GLUP.ClipPlane()' invalid number of arguments"
	);
    }
    if(
	!lua_isnumber(L,1) ||
	!lua_isnumber(L,2) ||
	!lua_isnumber(L,3) ||
	!lua_isnumber(L,4) 	
    ) {
	return luaL_error(
	    L, "'GLUP.ClipPlane()' invalid arguments"
	);
    }
    double eqn[4];
    eqn[0] = lua_tonumber(L,1);
    eqn[1] = lua_tonumber(L,2);
    eqn[2] = lua_tonumber(L,3);
    eqn[3] = lua_tonumber(L,4);
    glupClipPlane(eqn);
    return 0;
}

static int lua_glup_GetClipPlane(lua_State* L) {
    if(lua_gettop(L) != 0) {
	return luaL_error(
	    L, "'GLUP.GetClipPlane()' invalid number of arguments"
	);
    }
    double eqn[4];
    glupGetClipPlane(eqn);
    lua_pushnumber(L,eqn[0]);
    lua_pushnumber(L,eqn[1]);
    lua_pushnumber(L,eqn[2]);
    lua_pushnumber(L,eqn[3]);    
    return 4;
}

static int lua_glup_MatrixMode(lua_State* L) {
    if(lua_gettop(L) != 1) {
	return luaL_error(
	    L, "'GLUP.MatrixMode()' invalid number of arguments"
	);
    }
    if(!lua_isinteger(L,1)) {
	return luaL_error(
	    L, "'GLUP.MatrixMode()' invalid argument"
	);
    }
    glupMatrixMode(GLUPmatrix(lua_tointeger(L,1)));
    return 0;
}

static int lua_glup_GetMatrixMode(lua_State* L) {
    if(lua_gettop(L) != 0) {
	return luaL_error(
	    L, "'GLUP.MatrixMode()' invalid number of arguments"
	);
    }
    lua_pushinteger(L,glupGetMatrixMode());
    return 1;
}

static int lua_glup_PushMatrix(lua_State* L) {
    if(lua_gettop(L) != 0) {
	return luaL_error(
	    L, "'GLUP.PushMatrix()' invalid number of arguments"
	);
    }
    GEO::index_t* matrix_depth = nil;
    switch(glupGetMatrixMode()) {
	case GLUP_MODELVIEW_MATRIX:
	    matrix_depth = &lua_glup_modelview_depth;
	    break;
	case GLUP_PROJECTION_MATRIX:
	    matrix_depth = &lua_glup_projection_depth;	    
	    break;
	case GLUP_TEXTURE_MATRIX:
	    matrix_depth = &lua_glup_texture_depth;	    
	    break;
    }
    if(*matrix_depth >= lua_glup_max_depth) {
	return luaL_error(
	    L, "'GLUP.PushMatrix()' stack overflow"
	);
    }
    ++(*matrix_depth);
    glupPushMatrix();
    return 0;
}

static int lua_glup_PopMatrix(lua_State* L) {
    if(lua_gettop(L) != 0) {
	return luaL_error(
	    L, "'GLUP.PopMatrix()' invalid number of arguments"
	);
    }
    GEO::index_t* matrix_depth = nil;
    switch(glupGetMatrixMode()) {
	case GLUP_MODELVIEW_MATRIX:
	    matrix_depth = &lua_glup_modelview_depth;
	    break;
	case GLUP_PROJECTION_MATRIX:
	    matrix_depth = &lua_glup_projection_depth;	    
	    break;
	case GLUP_TEXTURE_MATRIX:
	    matrix_depth = &lua_glup_texture_depth;	    
	    break;
    }
    if(*matrix_depth == 0) {
	return luaL_error(
	    L, "'GLUP.PopMatrix()' stack underflow"
	);
    }
    --(*matrix_depth);
    glupPopMatrix();
    return 0;
}

static int lua_glup_GetMatrix(lua_State* L) {
    if(lua_gettop(L) != 1) {
	return luaL_error(
	    L, "'GLUP.GetMatrix()' invalid number of arguments"
	);
    }
    if(!lua_isinteger(L,1)) {
	return luaL_error(
	    L, "'GLUP.GetMatrix()' invalid argument"
	);
    }
    double M[16];
    glupGetMatrixdv(GLUPmatrix(lua_tointeger(L,1)),M);
    for(int i=0; i<16; ++i) {
	lua_pushnumber(L,M[i]);
    }
    return 16;
}

static int lua_glup_LoadIdentity(lua_State* L) {
    if(lua_gettop(L) != 0) {
	return luaL_error(
	    L, "'GLUP.GetMatrix()' invalid number of arguments"
	);
    }
    glupLoadIdentity();
    return 0;
}

static int lua_glup_LoadMatrix(lua_State* L) {
    if(lua_gettop(L) != 16) {
	return luaL_error(
	    L, "'GLUP.LoadMatrix()' invalid number of arguments"
	);
    }
    double M[16];
    for(int i=0; i<16; ++i) {
	if(!lua_isnumber(L,i+1)) {
	    return luaL_error(
		L, "'GLUP.LoadMatrix()' invalid argument"
	    );
	}
	M[i] = lua_tonumber(L,i+1);
    }
    glupLoadMatrixd(M);
    return 0;
}

static int lua_glup_MultMatrix(lua_State* L) {
    if(lua_gettop(L) != 16) {
	return luaL_error(
	    L, "'GLUP.MultMatrix()' invalid number of arguments"
	);
    }
    double M[16];
    for(int i=0; i<16; ++i) {
	if(!lua_isnumber(L,i+1)) {
	    return luaL_error(
		L, "'GLUP.MultMatrix()' invalid argument"
	    );
	}
	M[i] = lua_tonumber(L,i+1);
    }
    glupMultMatrixd(M);
    return 0;
}

static int lua_glup_Translate(lua_State* L) {
    if(lua_gettop(L) != 3) {
	return luaL_error(
	    L, "'GLUP.Translate()' invalid number of arguments"
	);
    }
    double T[3];
    for(int i=0; i<3; ++i) {
	if(!lua_isnumber(L,i+1)) {
	    return luaL_error(
		L, "'GLUP.Translate()' invalid argument"
	    );
	}
	T[i] = lua_tonumber(L,i+1);
    }
    glupTranslated(T[0],T[1],T[2]);
    return 0;
}

static int lua_glup_Scale(lua_State* L) {
    if(lua_gettop(L) != 3) {
	return luaL_error(
	    L, "'GLUP.Scale()' invalid number of arguments"
	);
    }
    double S[3];
    for(int i=0; i<3; ++i) {
	if(!lua_isnumber(L,i+1)) {
	    return luaL_error(
		L, "'GLUP.Scale()' invalid argument"
	    );
	}
	S[i] = lua_tonumber(L,i+1);
    }
    glupScaled(S[0],S[1],S[2]);
    return 0;
}

static int lua_glup_Rotate(lua_State* L) {
    if(lua_gettop(L) != 4) {
	return luaL_error(
	    L, "'GLUP.Rotate()' invalid number of arguments"
	);
    }
    double R[4];
    for(int i=0; i<4; ++i) {
	if(!lua_isnumber(L,i+1)) {
	    return luaL_error(
		L, "'GLUP.Rotate()' invalid argument"
	    );
	}
	R[i] = lua_tonumber(L,i+1);
    }
    glupRotated(R[0],R[1],R[2],R[3]);
    return 0;
}

static int lua_glup_Begin(lua_State* L) {
    if(lua_gettop(L) != 1) {
	return luaL_error(
	    L, "'GLUP.Begin()' invalid number of arguments"
	);
    }
    if(!lua_isinteger(L,1)) {
	return luaL_error(
	    L, "'GLUP.Begin()' argument should be an integer"
	);
    }
    if(lua_glup_primitive_active) {
	return luaL_error(
	    L, "'GLUP.Begin()' called twice without GLUP.End()"
	);
    }
    lua_glup_primitive_active = true;
    lua_Integer prim = lua_tointeger(L,1);
    glupBegin(GLUPprimitive(prim));
    return 0;
}

static int lua_glup_End(lua_State* L) {
    if(lua_gettop(L) != 0) {
	return luaL_error(
	    L, "'GLUP.End()' invalid number of arguments"
	);
    }
    if(!lua_glup_primitive_active) {
	return luaL_error(
	    L, "'GLUP.End()' called without GLUP.Begin()"
	);
    }
    lua_glup_primitive_active = false;
    glupEnd();
    return 0;
}

static int lua_glup_Vertex(lua_State* L) {
    double xyzw[4];
    if(!get_vec4(L,xyzw)) {
	return luaL_error(
	    L, "'GLUP.Vertex()' invalid arguments"
	);
    }
    glupVertex4dv(xyzw);
    return 0;
}

static int lua_glup_Color(lua_State* L) {
    double rgba[4];
    if(!get_vec4(L,rgba)) {
	return luaL_error(
	    L, "'GLUP.Color()' invalid arguments"
	);
    }
    glupColor4dv(rgba);
    return 0;
}

static int lua_glup_TexCoord(lua_State* L) {
    double xyzw[4];
    if(!get_vec4(L,xyzw)) {
	return luaL_error(
	    L, "'GLUP.TexCoord()' invalid arguments"
	);
    }
    glupTexCoord4dv(xyzw);
    return 0;
}

static int lua_glup_Normal(lua_State* L) {
    double xyzw[4];
    if(!get_vec4(L,xyzw)) {
	return luaL_error(
	    L, "'GLUP.Normal()' invalid arguments"
	);
    }
    glupNormal3dv(xyzw);
    return 0;
}

static int lua_glup_SetRegionOfInterest(lua_State* L) {
    if(lua_gettop(L) != 6) {
	return luaL_error(
	    L, "'GLUP.SetRegionOfInterest()' invalid number of arguments"
	);
    }
    if(
	!lua_isnumber(L,1) ||
	!lua_isnumber(L,2) ||
	!lua_isnumber(L,3) ||
	!lua_isnumber(L,4) ||	
	!lua_isnumber(L,5) ||
	!lua_isnumber(L,6) 
    ) {
	return luaL_error(
	    L, "'GLUP.SetRegionOfInterest()' arguments should be numbers"
	);
    }
    glup_viewer_set_region_of_interest(
	float(lua_tonumber(L,1)),
	float(lua_tonumber(L,2)),
	float(lua_tonumber(L,3)),
	float(lua_tonumber(L,4)),
	float(lua_tonumber(L,5)),
	float(lua_tonumber(L,6))	
    );
    return 0;
}

static int lua_glup_GetRegionOfInterest(lua_State* L) {
    if(lua_gettop(L) != 0) {
	return luaL_error(
	    L, "'GLUP.GetRegionOfInterest()' invalid number of arguments"
	);
    }
    float xm,ym,zm,xM,yM,zM;
    glup_viewer_get_region_of_interest(
	&xm, &ym, &zm, &xM, &yM, &zM
    );
    lua_pushnumber(L,double(xm));
    lua_pushnumber(L,double(ym));
    lua_pushnumber(L,double(zm));
    lua_pushnumber(L,double(xM));
    lua_pushnumber(L,double(yM));
    lua_pushnumber(L,double(zM));            
    return 6;
}

static int lua_glup_import(lua_State* L) {
    if(lua_gettop(L) != 1) {
	return luaL_error(
	    L, "'import()' invalid number of arguments"
	);
    }
    std::string k =
	std::string("lib/") + std::string(lua_tostring(L,1)) + ".lua";
    std::map<std::string, const char*>::iterator it =
	embedded_files.find(k);
    if(it == embedded_files.end()) {
	return luaL_error(
	    L, "'import()' unknown embedded file"
	);
    }
    const char* source = it->second;
    if(luaL_dostring(L,source)) {
	const char* msg = lua_tostring(L,-1);
	GEO::Logger::err("LUA") << msg << std::endl;
    }
    return 0;
}


static double t0 = 0.0;

static int lua_glup_ElapsedTime(lua_State* L) {
    if(lua_gettop(L) != 0) {
	return luaL_error(
	    L, "'GLUP.ElapsedTime()' invalid number of arguments"
	);
    }
    double result = 0.0;
    result = GEO::SystemStopwatch::now() - t0;
    lua_pushnumber(L,double(result));    
    return 1;
}

static int lua_glup_ResetViewer(lua_State* L) {
    if(lua_gettop(L) != 0) {
	return luaL_error(
	    L, "'GLUP.ResetViewer()' invalid number of arguments"
	);
    }
    glup_viewer_home();
    glup_viewer_disable(GLUP_VIEWER_CLIP);
    glup_viewer_disable(GLUP_VIEWER_EDIT_CLIP);
    glup_viewer_disable(GLUP_VIEWER_FIXED_CLIP);    
    glup_viewer_disable(GLUP_VIEWER_ROTATE_LIGHT);
    glup_viewer_enable(GLUP_VIEWER_BACKGROUND);    
    GEO::Application::instance()->set_lighting(true);
    GEO::Application::instance()->set_white_bg(true);
    t0 = GEO::SystemStopwatch::now();
    return 0;
}

static int lua_glup_ArcadeStyle(lua_State* L) {
    if(lua_gettop(L) != 0) {
	return luaL_error(
	    L, "'GLUP.ArcadeStyle()' invalid number of arguments"
	);
    }
    glup_viewer_home();
    glup_viewer_disable(GLUP_VIEWER_CLIP);
    glup_viewer_disable(GLUP_VIEWER_EDIT_CLIP);
    glup_viewer_disable(GLUP_VIEWER_FIXED_CLIP);    
    glup_viewer_disable(GLUP_VIEWER_ROTATE_LIGHT);
    glup_viewer_disable(GLUP_VIEWER_BACKGROUND);    
    GEO::Application::instance()->set_lighting(false);
    GEO::Application::instance()->set_white_bg(false);
    return 0;    
}

void init_lua_glup(lua_State* L) {
  lua_newtable(L);

  // Enable/Disable constants
  DECLARE_GLUP_CONSTANT(LIGHTING);
  DECLARE_GLUP_CONSTANT(VERTEX_COLORS);
  DECLARE_GLUP_CONSTANT(DRAW_MESH);
  DECLARE_GLUP_CONSTANT(CLIPPING);
  DECLARE_GLUP_CONSTANT(INDIRECT_TEXTURING);
  DECLARE_GLUP_CONSTANT(VERTEX_NORMALS);
  DECLARE_GLUP_CONSTANT(PICKING);
  DECLARE_GLUP_CONSTANT(ALPHA_DISCARD);

  // Colors
  DECLARE_GLUP_CONSTANT(FRONT_COLOR);
  DECLARE_GLUP_CONSTANT(BACK_COLOR);
  DECLARE_GLUP_CONSTANT(MESH_COLOR);
  DECLARE_GLUP_CONSTANT(FRONT_AND_BACK_COLOR);

  // Picking
  DECLARE_GLUP_CONSTANT(PICK_PRIMITIVE);
  DECLARE_GLUP_CONSTANT(PICK_CONSTANT);

  // Clipping
  DECLARE_GLUP_CONSTANT(CLIP_STANDARD);
  DECLARE_GLUP_CONSTANT(CLIP_WHOLE_CELLS);
  DECLARE_GLUP_CONSTANT(CLIP_STRADDLING_CELLS);
  DECLARE_GLUP_CONSTANT(CLIP_SLICE_CELLS);  

  // Matrices
  DECLARE_GLUP_CONSTANT(MODELVIEW_MATRIX);
  DECLARE_GLUP_CONSTANT(PROJECTION_MATRIX);
  DECLARE_GLUP_CONSTANT(TEXTURE_MATRIX);  

  // Primitives
  DECLARE_GLUP_CONSTANT(POINTS);
  DECLARE_GLUP_CONSTANT(LINES);
  DECLARE_GLUP_CONSTANT(TRIANGLES);
  DECLARE_GLUP_CONSTANT(QUADS);
  DECLARE_GLUP_CONSTANT(TETRAHEDRA);
  DECLARE_GLUP_CONSTANT(HEXAHEDRA);
  DECLARE_GLUP_CONSTANT(PRISMS);
  DECLARE_GLUP_CONSTANT(PYRAMIDS);
  DECLARE_GLUP_CONSTANT(CONNECTORS);
  DECLARE_GLUP_CONSTANT(SPHERES);  

  // Functions
  DECLARE_GLUP_FUNC(Enable);
  DECLARE_GLUP_FUNC(Disable);
  DECLARE_GLUP_FUNC(IsEnabled);
  DECLARE_GLUP_FUNC(SetColor);
  DECLARE_GLUP_FUNC(GetColor);
  DECLARE_GLUP_FUNC(LightVector);
  DECLARE_GLUP_FUNC(SetPointSize);
  DECLARE_GLUP_FUNC(GetPointSize);
  DECLARE_GLUP_FUNC(SetMeshWidth);
  DECLARE_GLUP_FUNC(GetMeshWidth);
  DECLARE_GLUP_FUNC(SetCellsShrink);
  DECLARE_GLUP_FUNC(GetCellsShrink);
  DECLARE_GLUP_FUNC(SetAlphaThreshold);
  DECLARE_GLUP_FUNC(GetAlphaThreshold);
  DECLARE_GLUP_FUNC(PickingMode);
  DECLARE_GLUP_FUNC(GetPickingMode);
  DECLARE_GLUP_FUNC(PickingId);
  DECLARE_GLUP_FUNC(GetPickingId);
  DECLARE_GLUP_FUNC(BasePickingId);
  DECLARE_GLUP_FUNC(GetBasePickingId);
  DECLARE_GLUP_FUNC(ClipMode);
  DECLARE_GLUP_FUNC(GetClipMode);
  DECLARE_GLUP_FUNC(ClipPlane);
  DECLARE_GLUP_FUNC(GetClipPlane);
  DECLARE_GLUP_FUNC(MatrixMode);
  DECLARE_GLUP_FUNC(GetMatrixMode);
  DECLARE_GLUP_FUNC(PushMatrix);
  DECLARE_GLUP_FUNC(PopMatrix);
  DECLARE_GLUP_FUNC(GetMatrix);
  DECLARE_GLUP_FUNC(LoadIdentity);
  DECLARE_GLUP_FUNC(LoadMatrix);
  DECLARE_GLUP_FUNC(MultMatrix);
  DECLARE_GLUP_FUNC(Translate);
  DECLARE_GLUP_FUNC(Scale);
  DECLARE_GLUP_FUNC(Rotate);      
  DECLARE_GLUP_FUNC(Begin);
  DECLARE_GLUP_FUNC(End);
  DECLARE_GLUP_FUNC(Vertex);
  DECLARE_GLUP_FUNC(Color);
  DECLARE_GLUP_FUNC(TexCoord);
  DECLARE_GLUP_FUNC(Normal);

  DECLARE_GLUP_FUNC(ElapsedTime);
  DECLARE_GLUP_FUNC(ResetViewer);
  DECLARE_GLUP_FUNC(ArcadeStyle);            

  DECLARE_GLUP_FUNC(SetRegionOfInterest);
  DECLARE_GLUP_FUNC(GetRegionOfInterest);  

  DECLARE_GLUP_COLOR("black", 0.0, 0.0, 0.0);
  DECLARE_GLUP_COLOR("white", 1.0, 1.0, 1.0);
  DECLARE_GLUP_COLOR("gray",  0.5, 0.5, 0.5);
  DECLARE_GLUP_COLOR("red", 1.0, 0.0, 0.0);
  DECLARE_GLUP_COLOR("green", 0.0, 1.0, 0.0);
  DECLARE_GLUP_COLOR("blue", 0.0, 0.0, 1.0);
  DECLARE_GLUP_COLOR("yellow", 1.0, 1.0, 0.0);
  DECLARE_GLUP_COLOR("pink", 1.0, 0.75, 0.793);            
  DECLARE_GLUP_COLOR("brown",0.6445,0.164,0.164);
  
  lua_setglobal(L, "GLUP");

  lState = L;
  LoadImguiBindings();

  lua_register(L, "import", lua_glup_import);

}

void adjust_lua_glup_state(lua_State* L) {
    GEO::geo_argused(L);

    if(glupCurrentContext() == nil) {
	return;
    }
    
    if(lua_glup_primitive_active) {
	glupEnd();
	lua_glup_primitive_active = false;
    }
    GLUPmatrix mode = glupGetMatrixMode();
    
    if(lua_glup_modelview_depth != 0) {
	glupMatrixMode(GLUP_MODELVIEW_MATRIX);
	while(lua_glup_modelview_depth != 0) {
	    glupPopMatrix();
	    --lua_glup_modelview_depth;
	}
    }
    
    if(lua_glup_projection_depth != 0) {
	glupMatrixMode(GLUP_PROJECTION_MATRIX);
	while(lua_glup_projection_depth != 0) {
	    glupPopMatrix();
	    --lua_glup_projection_depth;
	}
    }
    
    if(lua_glup_texture_depth != 0) {
	glupMatrixMode(GLUP_TEXTURE_MATRIX);
	while(lua_glup_texture_depth != 0) {
	    glupPopMatrix();
	    --lua_glup_texture_depth;
	}
    }
    
    glupMatrixMode(mode);
}

void register_embedded_lua_file(
   const char* filename, const char* data
) {
   embedded_files[std::string(filename)] = data;
}

void list_embedded_lua_files(
    std::vector<std::string>& filenames
) {
    filenames.clear();
    for(std::map<std::string, const char*>::iterator it =
	    embedded_files.begin();
	it != embedded_files.end(); ++it) {
	filenames.push_back(it->first);
    }
}

void get_embedded_lua_file(
    const std::string& filename, const char** data
) {
    std::map<std::string, const char*>::iterator it =
	embedded_files.find(filename);
    *data = (it == embedded_files.end()) ? nil : it->second;
}
