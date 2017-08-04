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

#include <exploragram/hexdom/hexdom_pipeline.h>
#include <exploragram/hexdom/preprocessing.h>
#include <exploragram/hexdom/FF.h>
#include <exploragram/hexdom/PGP.h>
#include <exploragram/hexdom/PGP_export.h>
#include <exploragram/hexdom/time_log.h>

namespace GEO {
    
    namespace HexdomPipeline {
	
	bool SetConstraints(Mesh*m, std::string& msg) {
            try {
                produce_hexdom_input(m, msg);
            }
            catch (const char* s) {
                plop(s);
                msg = std::string(s);
                logt.add_string("fail", msg);
                return false;
            }
            return true;
        }
	
	void FrameField(Mesh*m, bool smooth, bool brush) {
            GEO::FFopt ffopt(m);
            ffopt.FF_init();
            if (smooth) ffopt.FF_smooth();
            ffopt.compute_Bid_norm();
            if (brush)  ffopt.brush_frame();
        }

        //{algo} = {0: PGP, 1 : PGP with correction, 2 CubeCover}
	void Parameterization(Mesh*m, int algo, double PGP_max_scale_corr) {
            {// scope to allows destrction of attributes that are members of pgp
                PGPopt pgp(m);
                if (algo < 2) {
                    pgp.optimize_corr(PGP_max_scale_corr);
                    logt.add_step("pgp.optimize_PGP");
                    pgp.optimize_PGP();
                }
                else if (algo == 2) {
                    logt.add_step("cubcover");
                    pgp.cubcover();
                }
                else geo_assert_not_reached;
            }
            m->edges.attributes().delete_attribute_store("corr"); // ? can we really do that ?
            if (algo>0) { // we have launch cubecover
                m->edges.attributes().delete_attribute_store("ccgrp");
                m->edges.attributes().delete_attribute_store("ccid");
            }
        }
	
	void HexCandidates(Mesh*m, Mesh* result) {
            PGPopt pgp(m);
            pgp.export_hexes(result);
        }

	bool QuadDominant(Mesh*m, Mesh* chartmesh) {
            {
                PGPopt pgp(m);
                Attribute<GEO::vec2> uv(chartmesh->facet_corners.attributes(), "uv");
                Attribute<index_t> singtri(chartmesh->facets.attributes(), "singular");

                pgp.export_boundary_with_uv(chartmesh, uv, singtri);
                get_facet_stats(chartmesh, "export with uv");


                STEP(split_edges_by_iso_uvs(chartmesh, uv, singtri));
                get_facet_stats(chartmesh, "split edges");

                STEP(facets_split_ca_va(chartmesh, uv, singtri));
                get_facet_stats(chartmesh, "facets split");

                Attribute<int> chart(chartmesh->facets.attributes(), "chart");
                STEP(mark_charts(chartmesh, uv, chart, singtri));
                get_facet_stats(chartmesh, "mark charts");

            } // attention attribute references are invalidated in these steps

            STEP(simplify_quad_charts(chartmesh));
            get_facet_stats(chartmesh, "remove vertices");

            STEP(export_quadtri_from_charts(chartmesh));

            if (!surface_is_tetgenifiable(chartmesh)) {
                logt.add_string("fail", " tetgen is not able to remesh the quadtri");
                return false;
            }

            chartmesh->facet_corners.attributes().delete_attribute_store("uv");
            chartmesh->facets.attributes().delete_attribute_store("chart");
            chartmesh->facets.attributes().delete_attribute_store("quadelement");
            chartmesh->facets.attributes().delete_attribute_store("singular");
            return true;
        }

	void Hexahedrons(Mesh* quaddominant, Mesh* hexcandidates, Mesh* result) {
            result->copy(*hexcandidates);
            STEP(hex_set_2_hex_mesh(result, quaddominant));
        }
	
	bool Cavity(Mesh* quaddominant, Mesh* hexahedrons, Mesh* result) {
            result->copy(*quaddominant);
            STEP(merge_hex_boundary_and_quadtri(hexahedrons, result));
            if (result->facets.nb() > 0 && !surface_is_tetgenifiable(result))
                return false;
            return true;
        }

	void HexDominant(Mesh* cavity, Mesh* hexahedrons, Mesh* result, bool with_pyramid) {
            Mesh tets;
            tets.copy(*cavity);
            STEP(fill_cavity_with_tetgen(cavity, &tets, true, with_pyramid));
            result->copy(tets);
            result->facets.clear();
            STEP(add_hexes_to_tetmesh(hexahedrons, result));
            result->cells.connect();
            result->cells.compute_borders();
        }
    }
    
}

