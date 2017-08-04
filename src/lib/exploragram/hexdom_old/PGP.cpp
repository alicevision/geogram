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

#include <exploragram/hexdom/PGP.h>
#include <exploragram/hexdom/frame.h> 
#include <exploragram/hexdom/extra_connectivity.h>
#include <exploragram/hexdom/geometry.h>
#include <exploragram/hexdom/mesh_utils.h>
#include <exploragram/hexdom/time_log.h>

// [BEGIN_HEXMEN_ONLY]    
#include <exploragram/hexdom/mixed_constrained_solver.h>
// [END_HEXMEN_ONLY]

#include <geogram/NL/nl.h>
#include <geogram/bibliography/bibliography.h>

#include <algorithm>
#include <cmath>
#include <queue>

namespace GEO {

    PGPopt::PGPopt(Mesh* p_m) : m(p_m) {
        geo_cite("DBLP:journals/tog/SokolovRUL16");
        // init fast acces to edges from vertices
        if (m->edges.nb() == 0) compute_tet_edge_graph(m, v2e,false);
        else restore_v2e(m, v2e);        
        v2eopp = vector<vector<index_t> >(m->edges.nb(), vector<index_t>());
        FOR(e, m->edges.nb()) v2eopp[m->edges.vertex(e, 1)].push_back(e);
       
        //bind attributes
        U.bind(m->vertices.attributes(), "U");
        B.bind(m->vertices.attributes(), "B");
        corr.bind(m->edges.attributes(), "corr");
        tij.bind(m->edges.attributes(), "tij");
    }

    bool PGPopt::is_PGP_singular(index_t c, index_t lf) {   
        index_t v[3];
        FOR(i, 3) v[i] = m->cells.facet_vertex(c, lf, i);
        vec3i t(0, 0, 0);
        mat3 R;
        R.load_identity();
        FOR(i, 3) {
            AxisPermutation r = Rij(m, B, v[i], v[(i + 1) % 3]);
            bool inv;
            index_t e = edge_from_vertices(v[i], v[(i + 1) % 3], inv);
            t += R*(inv ? -(r.inverse()*tij[e]) : tij[e]);
            R = R* r.inverse().get_mat();
        }

        geo_assert(R.is_identity());
        return t[0] || t[1] || t[2];
    }


    bool PGPopt::tet_is_PGP_singular_fct(index_t t) {
        bool is_sing = false;
        FOR(f, 4) is_sing = is_sing || is_PGP_singular(t, f);
        return is_sing;
    }


    bool PGPopt::face_is_resp(index_t c, index_t lf) {
        if (m->cells.adjacent(c, lf) == NO_CELL) return  true;
        return m->cells.adjacent(c, lf) < c;
    }

    struct TriangleEdgesInSameBasis {
        bool try_new_triangle(PGPopt* pgp, index_t c, index_t lf) {
            if (!pgp->face_is_resp(c, lf)) { return false; }                // avoid doing twice the same work
            if (triangle_is_frame_singular(pgp->m, pgp->B, c, lf)) { return false; }                                                        // curl correction on singular face is meaningless
            // without chain basis change
            index_t vid[3];                                             // vertices index
            // init vertex indices
            FOR(e, 3)               vid[e] = pgp->m->cells.facet_vertex(c, lf, e);
            // init edges
            FOR(e, 3) {
                edge[e] = pgp->edge_from_vertices(vid[e], vid[next_mod(e, 3)], inv[e]);
                geo_assert(edge[e] != NOT_AN_ID);
            }
            // express corr of all edges in a common basis with edge_ap[e]
            FOR(e, 3) {
                if (inv[e]) edge_ap[e] = Rij(pgp->m, pgp->B, vid[0], vid[next_mod(e, 3)]);
                else        edge_ap[e] = Rij(pgp->m, pgp->B, vid[0], vid[e]);
            }
            return true;
        }


        index_t edge[3];
        bool inv[3];
        AxisPermutation edge_ap[3];
    };




    //     ___              __
    //    / _ \ _ _  ___   / _|___ _ _ _ __  ___
    //   | (_) | ' \/ -_) |  _/ _ \ '_| '  \(_-<
    //    \___/|_||_\___| |_| \___/_| |_|_|_/__/
    //



    vec3 PGPopt::wish_angle_corr(index_t e, bool inv) {
        vec3 c(0, 0, 0);
        if (corr.is_bound()) { // CubeCover ne l'utilise pas forcement
            c = corr[e];
            if (inv) {
                AxisPermutation chg = Rij(m, B, m->edges.vertex(e, 0), m->edges.vertex(e, 1));
                c = -(chg.inverse()  *c);
            }
        }
        return c;
    }

    vec3 PGPopt::wish_angle_edge_geom(index_t e, bool inv) {
        index_t org = m->edges.vertex(e, 0);
        index_t dest = m->edges.vertex(e, 1);
        if (inv) std::swap(org, dest);

        AxisPermutation ap=Rij(m,B,org,dest);

        mat3 frame = B[org] + Frame(B[dest]).apply_permutation(ap);
        frame *= 0.5;
        frame = invert_columns_norm(frame);
        vec3 angle = frame.transpose() * (m->vertices.point(dest) - m->vertices.point(org));
        return 2.*M_PI *angle ; // one cycle length is edgelength_
    }

    vec3 PGPopt::wish_angle(index_t e, bool inv) {
        return  wish_angle_edge_geom(e, inv) + PGPopt::wish_angle_corr(e, inv);
    }





    //             _   _       _
    //    ___ _ __| |_(_)_ __ (_)______   __ ___ _ _ _ _
    //   / _ \ '_ \  _| | '  \| |_ / -_) / _/ _ \ '_| '_|
    //   \___/ .__/\__|_|_|_|_|_/__\___|_\__\___/_| |_|
    //       |_|


    ///////////////////////////////////////////////////
    //optimize_corr
    ///////////////////////////////////////////////////
    // problem: the objective one form derived from the frame field is not close, leading to PGP singularities (T_junctions)
    // solution: this function computes a one form such that :
    //           *its norm is minimal
    //                       *adding it to the original objective one form, make it close... everywhere but on frame field singularities
    // param : max_corr_prop typically in [0,.5] (0 => no correction, .5 => take initial field into account for only 50%
    //
    // note: it is the same energy as in cubecover, exept that there is no integer constraints, and the frame singularities are ignored


    void PGPopt::optimize_corr(double max_corr_prop) {
        FOR(e, m->edges.nb()) corr[e] = vec3(0, 0, 0);
        if (max_corr_prop == 0) {
            return;
        }
        cubcover(true); 
            // clamp the result: very naive way to avoid too large corrections. May be improved.
            FOR(e, m->edges.nb()) {
                //   BasisChg chg = edge_basis_change(e, false);
                vec3 geom_form = wish_angle_edge_geom(e, false);
                vec3 corr_form = wish_angle_corr(e, false);
                double scale = 1.;
                double cl = corr_form.length();
                double gl = geom_form.length();
                if (cl > max_corr_prop* gl)  scale = max_corr_prop* gl / cl;
                FOR(d, 3) corr[e][d] = scale * corr[e][d];
            }
        }



	void PGPopt::move_U_to_corner() {
		Attribute<vec3> UC(m->cell_corners.attributes(), "U");
		Attribute<bool> has_param(m->cell_facets.attributes(), "has_param");

		FOR(c, m->cells.nb()) {

			bool all_faces_have_param = true;
			// flag has param
			FOR(lf, 4) {
				index_t f = m->cells.facet(c, lf);
				has_param[f] = true;
				if (triangle_is_frame_singular(m, B, c, lf))has_param[f] = false;
				else if (is_PGP_singular(c, lf)) has_param[f] = false;
				all_faces_have_param = all_faces_have_param && has_param[f];
			}


			index_t org = 0;

			// find a seed that will properly reconstruct all valid triangle param
			if (!all_faces_have_param) {
				bool can_be_ref[4] = { true, true, true, true };
				FOR(lf, 4) {
					index_t f = m->cells.facet(c, lf);
					if (has_param[f]) {
						bool local_can_be_ref[4] = { false, false, false, false };
						FOR(lv, 3) local_can_be_ref[m->cells.descriptor(c).facet_vertex[lf][lv]] = true;
						FOR(i, 4) can_be_ref[i] = can_be_ref[i] && local_can_be_ref[i];
					}
				}
				org = NOT_AN_ID;
				FOR(i, 4) if (can_be_ref[i]) org = i;
			}


			if (org != NOT_AN_ID)
				FOR(i, 4) {
				index_t corner = m->cells.corner(c, i);
				UC[corner] = U[m->cells.vertex(c, i)];
				if (i != org) {
					AxisPermutation change = Rij(m, B, m->cells.vertex(c, org), m->cells.vertex(c, i));
					UC[corner] = change.inverse()  * U[m->cells.vertex(c, i)];
					bool inv;
					index_t e = edge_from_vertices(m->cells.vertex(c, org), m->cells.vertex(c, i), inv);
					geo_assert(e != NOT_AN_ID);
					vec3i t2 = !inv ? tij[e] : -(change.inverse()*tij[e]);
					UC[corner] += vec3(t2[0], t2[1], t2[2]);
				}
			}
			else geo_assert(!has_param[m->cells.facet(c, 0)] && !has_param[m->cells.facet(c, 1)] && !has_param[m->cells.facet(c, 2)] && !has_param[m->cells.facet(c, 3)]);

		}
	}

	/*
     *             _   _       _          ___  ___ ___
     *    ___ _ __| |_(_)_ __ (_)______  | _ \/ __| _ \
     *   / _ \ '_ \  _| | '  \| |_ / -_) |  _/ (_ |  _/
     *   \___/ .__/\__|_|_|_|_|_/__\___|_|_|  \___|_|
     *       |_|
     */

    void PGPopt::optimize_PGP() {
        Attribute<vec3> lockU(m->vertices.attributes(), "lockU");
        // Create and initialize OpenNL context
        nlNewContext();
        nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
        nlSolverParameteri(NL_NB_VARIABLES, NLint(6 * m->vertices.nb()));

        nlBegin(NL_SYSTEM);
        FOR(v, m->vertices.nb())FOR(d, 3) if (std::fabs(lockU[v][d]) > 0)FOR(c, 2) {
            nlSetVariable(6 * v + 2 * d + c, 1 - c);
            nlLockVariable(6 * v + 2 * d + c);
        }

        nlBegin(NL_MATRIX);
        FOR(e, m->edges.nb()) {
            AxisPermutation ap = Rij(m, B, m->edges.vertex(e, 0), m->edges.vertex(e, 1));

            vec3 theta = wish_angle(e, false);
            FOR(d, 3) {
                double c = cos(theta[d]);
                double s = sin(theta[d]);
                index_t off0 = 6 * m->edges.vertex(e, 0) + 2 * d;

                nlBegin(NL_ROW);
                FOR(dd,3)  if (ap.get_mat()(dd,d)!=0) 
                        nlCoefficient(6 * m->edges.vertex(e, 1) + 2 * dd, -1.);
                nlCoefficient(off0, c);
                nlCoefficient(off0 + 1, s);
                nlEnd(NL_ROW);
                nlBegin(NL_ROW);
                FOR(dd, 3) 
                        nlCoefficient(6 * m->edges.vertex(e, 1) + 2 * dd+1, -ap.get_mat()(dd, d));
                nlCoefficient(off0, -s);
                nlCoefficient(off0 + 1, c);
                nlEnd(NL_ROW);
            }
        }

        nlEnd(NL_MATRIX);
        nlEnd(NL_SYSTEM);
        // Solve and get solution
        nlSolve();

        FOR(v, m->vertices.nb())  FOR(d, 3)
            U[v][d] = (.5 / M_PI) *  atan2(nlGetVariable(6 * v + 2 * d + 1), nlGetVariable(6 * v + 2 * d));

        nlDeleteContext(nlGetCurrent());

        snap_U_to_round();

        FOR(e, m->edges.nb()) {
            index_t i = m->edges.vertex(e, 0);
            index_t j = m->edges.vertex(e, 1);
            AxisPermutation rij = Rij(m, B, i, j);

            vec3 gij = wish_angle(e, false) / (2.*M_PI);

            FOR(d, 3) {
                tij[e][d] = int(round(-(rij.inverse().get_mat()*U[j])[d] - gij[d] + U[i][d]));
            }

        }
		move_U_to_corner();
    }










    //     _____      _           _____
    //    / ____|    | |         / ____|
    //   | |    _   _| |__   ___| |     _____   _____ _ __
    //   | |   | | | | '_ \ / _ \ |    / _ \ \ / / _ \ '__|
    //   | |___| |_| | |_) |  __/ |___| (_) \ V /  __/ |
    //    \_____\__,_|_.__/ \___|\_____\___/ \_/ \___|_|
    //
    //
    struct TetHalfedge {
        TetHalfedge(index_t p_cell, index_t p_org, index_t p_dest) {
            cell = p_cell; org = p_org; dest = p_dest;
        }
        index_t cell;
        index_t org;
        index_t dest;
    };
    index_t PGPopt::get_the_only_non_zero_lc(index_t c, index_t cf, Attribute<index_t>& CCedgeid) {
        int nbzero = 0;
        FOR(cfv, 3) {
            bool inv;
            index_t e = edge_from_vertices(m->cells.facet_vertex(c, cf, cfv), m->cells.facet_vertex(c, cf, next_mod(cfv, 3)), inv);
            geo_assert(e != NOT_AN_ID);
            if (CCedgeid[e] == NOT_AN_ID)
                nbzero++;
        }
        if (nbzero == 2 && !triangle_is_frame_singular(m, B, c, cf)) {
            FOR(cfv, 3) {
                bool inv;
                index_t e = edge_from_vertices(m->cells.facet_vertex(c, cf, cfv), m->cells.facet_vertex(c, cf, next_mod(cfv, 3)), inv);
                geo_assert(e != NOT_AN_ID);
                if (CCedgeid[e] != NOT_AN_ID) {
                    CCedgeid[e] = NOT_AN_ID;
                    return cfv;
                }
            }
        }
        return NOT_AN_ID;
    }

    index_t PGPopt::number_edges_on_cut_graph(Attribute<index_t>& CCedgeid) {

        Attribute<vec3> lockU(m->vertices.attributes(), "lockU");

        index_t alot = std::numeric_limits<index_t>::max() - 10;




        FOR(e, m->edges.nb()) CCedgeid[e] = alot;

        vector<index_t> offset_from_org;
        vector<index_t> dest;

        // compute covering tree
        cell_edges_in_RCS(m, offset_from_org, dest);
        std::vector<index_t > dist(m->vertices.nb(), alot);


        {
            FOR(seed, m->vertices.nb()) {
                if (dist[seed] != alot) continue;
                dist[seed] = 0;
                std::queue<index_t> queue;
                queue.push(seed);
                while (!queue.empty()) {
                    index_t v = queue.front();
                    queue.pop();
                    for (index_t ed = offset_from_org[v]; ed < offset_from_org[v + 1]; ed++) {
                        index_t vopp = dest[ed];
                        if (dist[vopp] == alot) {
                            bool inv;
                            index_t e = edge_from_vertices(v, vopp, inv);
                            geo_assert(e != NOT_AN_ID);
                            CCedgeid[e] = NOT_AN_ID;
                            dist[vopp] = dist[v] + 1;
                            queue.push(vopp);
                           
                            AxisPermutation ap = Rij(m,B,v,vopp);
                            lockU[vopp] = ap.inverse().get_mat()*lockU[vopp];
                            B[vopp] = Frame(B[vopp]).apply_permutation(ap);
                        }
                    }
                }
            }
        }

        std::queue<TetHalfedge> queue;
        // init
        FOR(c, m->cells.nb()) FOR(cf, 4) {
            index_t cfe = get_the_only_non_zero_lc(c, cf, CCedgeid);
            if (cfe != NOT_AN_ID)
                queue.push(TetHalfedge(c, m->cells.facet_vertex(c, cf, cfe), m->cells.facet_vertex(c, cf, next_mod(cfe, 3))));
        }
        while (!queue.empty()) {
            TetHalfedge cur = queue.front();
            queue.pop();
            index_t cir = cur.cell;
            do {
                FOR(cf, 4) {
                    index_t cfe = get_the_only_non_zero_lc(cir, cf, CCedgeid);
                    if (cfe != NOT_AN_ID)
                        queue.push(TetHalfedge(cir, m->cells.facet_vertex(cir, cf, cfe), m->cells.facet_vertex(cir, cf, next_mod(cfe, 3))));
                }
                index_t next_cir = next_cell_around_oriented_edge(m, cir, cur.org, cur.dest);

                if (next_cir == NOT_AN_ID) {
                    next_cir = cir;
                    while (next_cell_around_oriented_edge(m, next_cir, cur.dest, cur.org) != NOT_AN_ID)
                        next_cir = next_cell_around_oriented_edge(m, next_cir, cur.dest, cur.org);
                }

                cir = next_cir;
            } while (cir != cur.cell);
        }
        // init m->CCedgeid
        index_t nb_edge_var = 0;
        FOR(e, m->edges.nb()) {
            if (CCedgeid[e] == alot) {
                CCedgeid[e] = nb_edge_var;
                nb_edge_var++;
            }
        }

        return nb_edge_var;

    }

// [BEGIN_HEXMEN_ONLY]    
#ifdef WITH_CUBECOVER
    void PGPopt::cubcover(bool compute_only_corr) {
        Attribute<vec3> lockU(m->vertices.attributes(), "lockU");
        Attribute<index_t> CCedgeid(m->edges.attributes(), "ccid");

        index_t nb_edge_var = number_edges_on_cut_graph(CCedgeid);
        index_t nb_vars = 3 * m->vertices.nb() + 3 * nb_edge_var;

        GEO::Logger::out("HexDom")  << "There are " << nb_edge_var << " edges on the cut graph" <<  std::endl;

        MatrixMixedConstrainedSolver intsolver(nb_vars);

        // define which variables are integers
        // 1/ period jumps
        for (index_t i = 3 * m->vertices.nb(); i < nb_vars; i++) intsolver.is_integer(i);
        // 2/ some coordinates on boundary
        if (!compute_only_corr) FOR(v, m->vertices.nb()) FOR(d, 3) if (fabs(lockU[v][d]) > .1) intsolver.is_integer(3 * v + d);


        intsolver.start_linear_contraints();

	{
            plop("set boundary equality constraints");

            vector<bool> border_edge(m->edges.nb(), false);
            FOR(c, m->cells.nb())FOR(cf, 4) if (m->cells.adjacent(c, cf) == NO_CELL) FOR(cfv, 3) {
                bool inv;
                border_edge[edge_from_vertices(m->cells.facet_vertex(c, cf, cfv), m->cells.facet_vertex(c, cf, next_mod(cfv, 3)), inv)] = true;
            }
            FOR(e, m->edges.nb()) {
                if (!border_edge[e]) continue;
                index_t v = m->edges.vertex(e, 0);
                index_t vopp = m->edges.vertex(e, 1);
                vec3 constr[2] = { lockU[v], lockU[vopp] };

                AxisPermutation ap = Rij(m,B,v,vopp);
                constr[1] = ap.inverse()*constr[1];
                FOR(d, 3) if (std::abs(constr[0][d]) > 0 && std::abs(constr[1][d]) > 0) {
                    vec3 edg = normalize(X(m)[m->edges.vertex(e, 1)] - X(m)[m->edges.vertex(e, 0)]);
                    if (std::abs(dot(normalize(col(B[v], d)), edg)) < .2) { //if the edge is more of less orthogonal to the coordinate to be locked
                        intsolver.begin_constraint();
                        intsolver.add_constraint_coeff(3 * v + d, -1);
                        FOR(d_opp, 3)
                            intsolver.add_constraint_coeff(3 * vopp + d_opp, ap.get_mat()(d_opp, d));
                        index_t eid = CCedgeid[e];
                        if (eid != index_t(-1)) intsolver.add_constraint_coeff(3 * m->vertices.nb() + 3 * eid + d, 1);// need to check that... it's not a new constraint ?
                        intsolver.end_constraint();
                    }
                }



            }
        }

        logt.add_step("SYMBOLIC set closed one form");

        TriangleEdgesInSameBasis com;
        FOR(c, m->cells.nb())FOR(cf, 4) {
            if ((c * 4 + cf) % 100000 == 0)  GEO::Logger::out("HexDom")  << " SYMBOLIC constraint for tet face  " << c * 4 + cf << " / " << 4 * m->cells.nb() <<  std::endl;

            if (!com.try_new_triangle(this, c, cf)) continue;
            FOR(dim, 3) {
                intsolver.begin_constraint();
                FOR(e, 3) {
                    index_t eid = CCedgeid[com.edge[e]];
                    if (eid != index_t(-1)) {
                        if (!com.inv[e])
                            FOR(dd,3)
                            intsolver.add_constraint_coeff(3 * m->vertices.nb() + 3 * eid + dd, com.edge_ap[e].get_mat()(dd, dim));
                        else
                            FOR(dd, 3)
                            intsolver.add_constraint_coeff(3 * m->vertices.nb() + 3 * eid + dd, -com.edge_ap[e].get_mat()(dd, dim));
                    }
                }
                intsolver.end_constraint();
            }
        }
        intsolver.end_linear_constraints();


        FOR(iter, 2) {
            intsolver.start_iter(iter);
            FOR(e, m->edges.nb()) {
                //AxisPermutation ap = edge_axis_permutation(e, false);
                AxisPermutation ap = Rij(m, B, m->edges.vertex(e, 0), m->edges.vertex(e, 1));
                vec3 theta = (1. / (2.*M_PI)) * wish_angle_edge_geom(e, false);
                FOR(d, 3) {
                    index_t off0 = 3 * m->edges.vertex(e, 0) + d;
                    intsolver.begin_energy();
                    intsolver.add_energy_coeff(off0, -1.);
                    FOR(dd,3) 
                        intsolver.add_energy_coeff(3 * m->edges.vertex(e, 1) + dd, ap.get_mat()(dd,d));
                    index_t eid = CCedgeid[e];
                    if (eid != index_t(-1))  intsolver.add_energy_coeff(3 * m->vertices.nb() + 3 * eid + d, 1);
                    intsolver.add_energy_rhs(-theta[d]);
                    intsolver.end_energy();
                }
            }
            intsolver.end_iter();
        }

        intsolver.check_constraints();



        //export debug info
        {
            Attribute<index_t> CCgrp(m->edges.attributes(), "ccgrp");
            FOR(e, m->edges.nb()) {
                CCgrp[e] = NOT_AN_ID;
                if (CCedgeid[e] != NOT_AN_ID) {
                    index_t mi = Numeric::Limits<Numeric::uint32>::max();
                    FOR(d, 3) mi = geo_min(mi, intsolver.G.line(3 * m->vertices.nb() + 3 * CCedgeid[e] + d)[0].index);
                    CCgrp[e] = mi;
                }
            }
        }




        FOR(v, m->vertices.nb())  FOR(d, 3) U[v][d] = intsolver.value(3 * v + d);

		snap_U_to_round();


        {
            logt.add_step("update correction 1-form for pts reconstruction");
            FOR(e, m->edges.nb()) {
                //AxisPermutation ap = edge_axis_permutation(e, false);
                AxisPermutation ap = Rij(m, B, m->edges.vertex(e, 0), m->edges.vertex(e, 1));
                vec3 dx(0, 0, 0);
                FOR(d, 3) {
                    index_t off0 = 3 * m->edges.vertex(e, 0) + d;
                    FOR(d1, 3) 
                        dx[d] += ap.get_mat()(d1, d) * intsolver.value(3 * m->edges.vertex(e, 1) + d1) ;
                    dx[d] -= intsolver.value(off0);                    
                    index_t eid = CCedgeid[e];
                    if (eid != index_t(-1))  dx[d] += intsolver.value(3 * m->vertices.nb() + 3 * eid + d);
                }
                FOR(d, 3) corr[e][d] = -2.*M_PI*dx[d] - wish_angle_edge_geom(e, false)[d];


                { // tij from corr
                    index_t i = m->edges.vertex(e, 0);
                    index_t j = m->edges.vertex(e, 1);
                    AxisPermutation rij = Rij(m, B, i, j);

                    vec3 gij = wish_angle(e, false) / (2.*M_PI);

                    FOR(d, 3) {
                        tij[e][d] = int(round(-(rij.inverse()*U[j])[d] - gij[d] + U[i][d]));
                    }
                }
            }
        }
		move_U_to_corner();


    }
#else
// [END_HEXMEN_ONLY]            
    void PGPopt::cubcover(bool compute_only_corr) {
	geo_argused(compute_only_corr);
	Logger::warn("PGP") << "Cubecover not implemented in the public version"
			    << std::endl;
    }
// [BEGIN_HEXMEN_ONLY]                
#endif
// [END_HEXMEN_ONLY]    
    
}
