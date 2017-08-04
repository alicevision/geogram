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

#include <exploragram/hexdom/preprocessing.h>

namespace GEO {

    void compute_input_constraints(Mesh* m) {

        Attribute<mat3> B(m->vertices.attributes(), "B");
        Attribute<vec3> lockB(m->vertices.attributes(), "lockB");// how many vectors are locked
        Attribute<vec3> U(m->vertices.attributes(), "U");
        Attribute<vec3> lockU(m->vertices.attributes(), "lockU");// how many dimensions are locked
	// init all normal for each vertex
        vector<vector<vec3> > normals(m->vertices.nb());
        vector<vector<double> > weight(m->vertices.nb());

        FOR(c, m->cells.nb()) FOR(cf, 4) {
            if ((m->cells.adjacent(c, cf) != NO_CELL)) continue;
            vec3 n = tet_facet_cross(m, c, cf);
            if (n.length2() > 1e-10) {
                FOR(cfv, 3) {
                    normals[m->cells.facet_vertex(c, cf, cfv)].push_back(normalize(n));
                    weight[m->cells.facet_vertex(c, cf, cfv)].push_back(n.length());
                }
            }
        }



        FOR(v, m->vertices.nb()) {
            vector<vec3>& n = normals[v];
            vector<double>& w = weight[v];

            B[v].load_identity();
            lockB[v] = vec3(0, 0, 0);
            lockU[v] = vec3(0, 0, 0);
            if (n.size() > 0) {
                B[v] = Frame::representative_frame(n, w);// rot_to_B(representative_frame(n, w));
                AxisPermutation ap;
                ap.make_col2_equal_to_z(B[v], n[0]);
                B[v] = Frame(B[v]).apply_permutation(ap);
                FOR(i, n.size()) FOR(a, 3) 
                    if (std::abs(dot(col(B[v], a), n[i])) > .7) {
                        lockU[v][a] = 1;
                        lockB[v][a] = 1;
                    }
                if (lockB[v].length2() > 1) lockB[v] = vec3(1, 1, 1);
                U[v] = vec3(0, 0, 0);
            }
            //B[v].load_identity();
            //lockB[v] = vec3(0, 0, 0);
            //if (n.size() > 0) {
            //    vec3 pca_result;
            //    vec3 frame_result[3];

            //    // compute PCA and see how it explains the distribution
            //    UncenteredPCA3D pca;
            //    pca.begin_points();
            //    for (index_t i = 0; i < n.size(); i++) pca.point(n[i], w[i]);
            //    pca.end_points();
            //    pca_result = normalize(pca.axis[0]);

            //    // fast output if the PCA is very good
            //    bool PCA_is_good_enough = true;
            //    FOR(i, n.size()) PCA_is_good_enough = PCA_is_good_enough && (std::abs(dot(n[0], n[i])) > .7);
            //    if (PCA_is_good_enough) {
            //        Frame(B[v]).make_z_equal_to(pca_result);
            //        lockB[v] = vec3(0, 0, 1);
            //    }
            //    else {// PCA was not sufficient
            //        B[v] = Frame::representative_frame(n, w);// rot_to_B(representative_frame(n, w));
            //        
            //        lockB[v] = vec3(1, 1, 1);
            //    }
            //}
            //lockU[v] = vec3(0, 0, 0);
            //FOR(i, n.size()) FOR(a, 3)
            //    if (std::abs(dot(col(B[v], a), n[i])) > .7) lockU[v][a] = 1;
            //U[v] = vec3(0, 0, 0);
        }

    }






    void reorder_vertices_according_to_constraints(Mesh* m) {
        Attribute<vec3> lockB(m->vertices.attributes(), "lockB");// how many vectors are locked

        GEO::vector<index_t > ind_map(m->vertices.nb());
        index_t n = 0;
        index_t num_l_v = 0;
        index_t num_ln_v = 0;
        FOR(v, m->vertices.nb()) if (lockB[v][0] == 1) { ind_map[n] = v; n++; }
        num_l_v = n;
        FOR(v, m->vertices.nb()) if (lockB[v][0] == 0 && lockB[v][2] == 1) { ind_map[n] = v; n++; }
        num_ln_v = n;
        FOR(v, m->vertices.nb()) if (lockB[v][2] == 0) { ind_map[n] = v; n++; }

        for (index_t i = 0; i < num_l_v; i++)                   geo_assert(lockB[ind_map[i]][0] == 1);
        for (index_t i = num_l_v; i < num_ln_v; i++)            geo_assert(lockB[ind_map[i]][2] == 1);
        for (index_t i = num_ln_v; i < m->vertices.nb(); i++)   geo_assert(lockB[ind_map[i]][2] == 0);

        geo_assert(num_l_v <= num_ln_v && num_ln_v <= m->vertices.nb());
	
        compute_Hilbert_order(m->vertices.nb(), m->vertices.point_ptr(0), ind_map, 0, num_l_v);
        compute_Hilbert_order(m->vertices.nb(), m->vertices.point_ptr(0), ind_map, num_l_v, num_ln_v);
        compute_Hilbert_order(m->vertices.nb(), m->vertices.point_ptr(0), ind_map, num_ln_v, m->vertices.nb());
	
        m->vertices.permute_elements(ind_map);          // note: it also updates the cell_corners.vertex... and invert ind_map :(
    }


    void produce_hexdom_input(Mesh* m,std::string& error_msg) {

        // some check on the input mesh
        m->edges.clear();
        m->vertices.remove_isolated();
       


        if (have_negative_tet_volume(m)) {
            throw ("contains tets with negative volume");
        }

        if (!m->cells.are_simplices()) {
            throw ("cells contains non tet elements");
        }
   
        if (!volume_boundary_is_manifold(m, error_msg)) {
            throw (error_msg.c_str());
        }

        if (!volume_is_tetgenifiable(m)) {
            throw (" tetgen is not able to remesh the volume from its boundary");
        }

        // add some attributes

        compute_input_constraints(m);
        
        // compute scale
        double wanted_edge_length = get_cell_average_edge_size(m);
        
        Attribute<mat3> B(m->vertices.attributes(), "B"); 
        FOR(v, m->vertices.nb()) FOR(ij, 9) B[v].data()[ij] *= wanted_edge_length;

        reorder_vertices_according_to_constraints(m);
    }


}

