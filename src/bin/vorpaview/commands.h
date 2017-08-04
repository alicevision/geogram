/*
 *  Copyright (c) 2012-2014, Bruno Levy
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

#ifndef VORPAVIEW_COMMANDS_H
#define VORPAVIEW_COMMANDS_H

namespace GEO {
    class Mesh;
    class MeshGfx;
}

/**
 * \brief Initializes Vorpaview commands.
 * \details Should be called once, at initialization.
 * \param[in] M a reference to the Mesh displayed by Vorpaview
 * \param[in] M_gfx a reference to the MeshGfx used to display
 *  the Mesh
 * \param[in] show_vertices_in a reference to the vertices 
 *  visibility flag of Vorpaview
 * \param[in] show_surface_in a reference to the surface
 *  visibility flag of Vorpaview
 * \param[in] show_volume_in a reference to the volume
 *  visibility flag of Vorpaview
 * \param[in] show_console_in a reference to the visibility
 *  flag of the console
 */
void vorpaview_commands_init(
    GEO::Mesh& M, GEO::MeshGfx& M_gfx,
    bool& show_vertices_in, bool& show_surface_in, bool& show_volume_in,
    bool& show_console_in
);

/**
 * \brief Manages the menus for the commands.
 * \details Needs to be called in the function that manages
 *  Vorpaview's menubar.
 */
void vorpaview_commands_menus();

/**
 * \brief Manages the GUI for command's arguments.
 */
void vorpaview_commands_gui();

#endif
