/*
 *  OGF/Graphite: Geometry and Graphics Programming Library + Utilities
 *  Copyright (C) 2000-2015 INRIA - Project ALICE
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact for Graphite: Bruno Levy - Bruno.Levy@inria.fr
 *  Contact for this Plugin: Nicolas Ray - nicolas.ray@inria.fr
 *
 *     Project ALICE
 *     LORIA, INRIA Lorraine,
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs.
 *
 * As an exception to the GPL, Graphite can be linked with the following
 * (non-GPL) libraries:
 *     Qt, tetgen, SuperLU, WildMagic and CGAL
 */

#ifndef H_HEXDOM_ALGO_PGP_EXPORT_H
#define H_HEXDOM_ALGO_PGP_EXPORT_H

#include <exploragram/basic/common.h>
#include <exploragram/hexdom/PGP.h>


namespace GEO{

        // TODO put attributes into parameters list
        
        void EXPLORAGRAM_API split_edges_by_iso_uvs(Mesh* hex, Attribute<vec2>& uv, Attribute<index_t> &singular);
	
        void EXPLORAGRAM_API facets_split_ca_va(Mesh* m, Attribute<vec2>& uv, Attribute<index_t> &singular);
	
        void EXPLORAGRAM_API mark_charts(Mesh* m, Attribute<vec2>& uv, Attribute<int>& charts, Attribute<index_t> &singular);
	
        void EXPLORAGRAM_API hex_set_2_hex_mesh(Mesh* hex,Mesh* quadtri=NULL);
	
        void EXPLORAGRAM_API merge_hex_boundary_and_quadtri(Mesh* hex, Mesh* quad);
	
        void EXPLORAGRAM_API export_quadtri_from_charts(Mesh* hex);

        void EXPLORAGRAM_API fill_cavity_with_tetgen(Mesh* input, Mesh* tri, bool propagate_hexvid, bool with_pyramid);

        void EXPLORAGRAM_API simplify_quad_charts(Mesh* m);

        void EXPLORAGRAM_API kill_intersecting_hexes(Mesh* hex);

        void EXPLORAGRAM_API add_hexes_to_tetmesh(Mesh* hex, Mesh* tet_mesh);

        // for debug
        void EXPLORAGRAM_API export_tb_faces_DEBUG(Mesh* hex, Mesh* tbfaces);
}

#endif
