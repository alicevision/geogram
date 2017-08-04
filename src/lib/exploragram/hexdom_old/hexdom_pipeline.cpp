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
#include <exploragram/hexdom/hex_candidates.h> 
#include <exploragram/hexdom/quad_dominant.h> 
#include <exploragram/hexdom/hex.h> 
#include <exploragram/hexdom/cavity.h> 
#include <exploragram/hexdom/hex_dominant.h>

#include <exploragram/hexdom/PGP_export.h>
#include <exploragram/hexdom/time_log.h>

#include <geogram/points/colocate.h>

// [BEGIN_HEXMEN_ONLY]
#ifdef GEOGRAM_WITH_VORPALINE
#include <vorpalib/mesh/mesh_tet2hex.h>
#define HAS_TET2HEX
#endif
// [END_HEXMEN_ONLY]

namespace GEO {

	namespace HexdomPipeline {

#define STEP(funcname,args) { \
		logt.add_step(#funcname); \
		funcname args; \
}

		bool SetConstraints(Mesh*m, std::string& msg,bool hilbert_sort ) {
			try {
				STEP(produce_hexdom_input,(m, msg, hilbert_sort ));
			}
			catch (const char* s) {
				plop(s);
				msg = std::string(s);
				logt.add_string("fail", msg);
				return false;
			}
			return true;
		}

		void FrameField(Mesh*m, bool smooth) {
			GEO::FFopt ffopt(m);
			STEP(ffopt.FF_init,(true));
			if (smooth) 
				STEP(ffopt.FF_smooth,());
			STEP(ffopt.compute_Bid_norm,());
			//uncomment may ease debugging...  STEP(ffopt.brush_frame,());
		}

		//{algo} = {0: CubeCover, 1 : PGP with correction, 2 PGP}
		void Parameterization(Mesh*m, int algo, double PGP_max_scale_corr) {
			{// scope to allows destruction of attributes that are members of pgp
				PGPopt pgp(m);
#ifdef WITH_CUBECOVER				
				if (algo > 0) {
					STEP(pgp.optimize_corr,(PGP_max_scale_corr));
					STEP(pgp.optimize_PGP,());
				}
				else if (algo == 0) {
					STEP(pgp.cubcover,());
				}
				else geo_assert_not_reached;
#else
				geo_argused(PGP_max_scale_corr);
				if(algo != 2) {
				    Logger::warn("PGP")
					<< "cubecover/scale correction not available in the public version"
					<< std::endl;
				    Logger::warn("PGP")
					<< "falling back with PGP without scale correction"
					<< std::endl;
				}
				algo = 2;
				STEP(pgp.optimize_corr,(0.0));
				STEP(pgp.optimize_PGP,());				
#endif				
			}

			m->edges.attributes().delete_attribute_store("corr"); // ? can we really do that ?
			if (algo <2) { // we have launch cubecover
				m->edges.attributes().delete_attribute_store("ccgrp");
				m->edges.attributes().delete_attribute_store("ccid");
			}
			m->vertices.attributes().delete_attribute_store("U");
			m->edges.clear();
		}

		void HexCandidates(Mesh*m, Mesh* result) {
			STEP(export_hexes,(m,result));
		}

		bool QuadDominant(Mesh*m, Mesh* chartmesh) {

			STEP(export_boundary_with_uv,(m, chartmesh, "uv", "singular"));
			//get_facet_stats(chartmesh, "export with uv");

			STEP(split_edges_by_iso_uvs, (chartmesh, "uv", "singular"));

			STEP(facets_split,(chartmesh, "uv", "singular"));

			STEP(mark_charts,(chartmesh, "uv", "chart", "singular"));

			STEP(simplify_quad_charts,(chartmesh));

			STEP(export_quadtri_from_charts,(chartmesh));

			if (!surface_is_tetgenifiable(chartmesh)) {
				logt.add_string("fail", " tetgen is not able to remesh the quadtri");
				return false;
			}

			chartmesh->facet_corners.attributes().delete_attribute_store("uv");
			chartmesh->facets.attributes().delete_attribute_store("chart");
			chartmesh->facets.attributes().delete_attribute_store("quadelement");
			chartmesh->facets.attributes().delete_attribute_store("singular");
			chartmesh->vertices.attributes().delete_attribute_store("quadcorners");
			chartmesh->facets.attributes().delete_attribute_store("orig_tri_fid");
			chartmesh->facet_corners.attributes().delete_attribute_store("isovalue");
			return true;
		}

		void Hexahedrons(Mesh* quaddominant, Mesh* hexcandidates, Mesh* result) {
			result->copy(*hexcandidates);
			STEP(hex_set_2_hex_mesh,(result, quaddominant));
		}

		bool Cavity(Mesh* quaddominant, Mesh* hexahedrons, Mesh* result) {
			result->copy(*quaddominant);
			STEP(merge_hex_boundary_and_quadtri,(hexahedrons, result));

			if (result->facets.nb() > 0 && !surface_is_tetgenifiable(result)) {
				logt.add_string("fail", "empty cavity, is it normal?");
				return false;
			}

			return true;
		}

		void HexDominant(Mesh* cavity, Mesh* hexahedrons, Mesh* result, bool with_pyramid,bool baudoin_carrier) {

// [BEGIN_HEXMEN_ONLY]
#ifdef HAS_TET2HEX
			if (baudoin_carrier) {
				Mesh tets;
				tets.copy(*cavity);
				Baudoin_mesher(&tets);
				STEP(fill_cavity_with_tetgen, (cavity, &tets, true, with_pyramid));
				tets.facets.clear();
				tets.edges.clear();
				mesh_tet2hex(tets, *result, 0.5, 0.6, false, false);
				//STEP(add_hexes_to_tetmesh, (hexahedrons, result)); without  check
				index_t off_v = result->vertices.create_vertices(hexahedrons->vertices.nb());
				FOR(v, hexahedrons->vertices.nb()) X(result)[off_v + v] = X(hexahedrons)[v];

				index_t off_c = result->cells.create_hexes(hexahedrons->cells.nb());

				//vector<index_t> old2new(hexahedrons->vertices.nb());
				//Geom::colocate(hexahedrons->vertices.point_ptr(0), 3, hexahedrons->vertices.nb(), old2new, 1e-10);

				FOR(c, hexahedrons->cells.nb())FOR(cv, 8) {
					index_t vid = hexahedrons->cells.vertex(c, cv);
					result->cells.set_vertex(off_c + c, cv, off_v+vid);
				}

				result->facets.clear();
				result->cells.connect();
				result->cells.compute_borders();
			} else
#endif			    
// [END_HEXMEN_ONLY]
			{
			    #ifndef HAS_TET2HEX
			    if(baudoin_carrier) {
				Logger::warn("hexdom") << "This version does not have Vorpaline" << std::endl;
				Logger::warn("hexdom") << "Ignored flag: Carrier-Baudouin algo." << std::endl;				
				Logger::warn("hexdom") << "(filling cavity with tets, no recombination)" << std::endl;				
			    }
			    #endif
				Mesh tets;
				tets.copy(*cavity);
				STEP(fill_cavity_with_tetgen,(cavity, &tets, true, with_pyramid));
				result->copy(tets);
				result->facets.clear();
				STEP(add_hexes_to_tetmesh,(hexahedrons, result));
				result->cells.connect();
				result->cells.compute_borders();
			}
		}
	}

}

