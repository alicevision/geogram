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

#ifndef GEOGRAM_GFX_GLUP_VIEWER_GLUP_VIEWER_LUA
#define GEOGRAM_GFX_GLUP_VIEWER_GLUP_VIEWER_LUA

#include <geogram_gfx/basic/common.h>
extern "C" {
#include <geogram/third_party/lua/lua.h>
}

/**
 * \brief Registers all LUA extension functions.
 * \details This registers LUA wrappers for GLUP and
 *  ImGUI.
 * \param[in] L a pointer to the LUA state.
 */
void GEOGRAM_GFX_API init_lua_glup(lua_State* L);


/**
 * \brief Makes sure GLUP is in a valid state.
 * \details Restores the previous depth of matrix stacks
 *  and terminates pending GLUP primitives. This makes sure
 *  GLUP is in a valid state even if there was an error in
 *  LUA code.
 */
void GEOGRAM_GFX_API adjust_lua_glup_state(lua_State* L);


void GEOGRAM_GFX_API register_embedded_lua_file(
   const char* filename, const char* data
);

void GEOGRAM_GFX_API list_embedded_lua_files(
    std::vector<std::string>& filenames
);

void GEOGRAM_GFX_API get_embedded_lua_file(
    const std::string& filename, const char** data
);


#endif

