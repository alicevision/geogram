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
#include <exploragram/hexdom/mixed_constrained_solver.h>
#include <exploragram/hexdom/time_log.h>

#include <geogram/points/colocate.h>
#include <geogram/NL/nl.h>

#include <algorithm>
#include <cmath>
#include <queue>

namespace GEO {

    PGPopt::PGPopt(Mesh* p_m) : m(p_m) {
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

    void PGPopt::make_compatible(index_t ref, index_t v) {
        // update U[v] and B[v] (rotation part)
        AxisPermutation change = Rij(m, B, ref, v);
        U[v] = change.inverse()  * U[v];
        B[v] = Frame(B[v]).apply_permutation(change);


        // just to check
        Attribute<vec3> lockU(m->vertices.attributes(), "lockU");
        lockU[v] = change.inverse().get_mat()*lockU[v];

        //update U[v]           (translation part)
        bool inv;
        index_t e = edge_from_vertices(ref, v, inv);
        geo_assert(e != NOT_AN_ID);
        vec3i t2 = !inv ? tij[e] : -(change.inverse()*tij[e]);
        U[v] += vec3(t2[0], t2[1], t2[2]);

        // update tij around v
        {
            index_t start = v2e[v];
            index_t end = m->edges.nb();
            if (v + 1 < m->vertices.nb()) end = v2e[v + 1];
            for (index_t cir = start; cir < end; cir++) {
                tij[cir] = change.inverse()*tij[cir] + t2;
            }

            FOR(i, v2eopp[v].size()) {
                index_t e2 = v2eopp[v][i];
                AxisPermutation change2 = Rij(m, B, m->edges.vertex(e2, 0), v);
                tij[e2] -= change2.inverse()*t2;
            }
        }
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
    }




    //    ___                   _
    //   | __|_ ___ __  ___ _ _| |_
    //   | _|\ \ / '_ \/ _ \ '_|  _|
    //   |___/_\_\ .__/\___/_|  \__|
    //           |_|

    vec3 PGPopt::change_tet_basis(vec3 in, vec3  P[4], vec3 P_img[4]) { // TODO sortir ça + trglgrad dans un fichier apart
        Matrix<12,double> M;
        M.load_zero();
        double RHS[12];
        double X[12];

        index_t cur_row = 0; // in fact cur_row can be replaced by v*3+dimx
        FOR(v, 4)  FOR(dimx, 3) {
            FOR(dimu, 3)
                M(4 * dimx + dimu, cur_row) = P[v][dimu];
            M(4 * dimx + 3, cur_row) = 1.0;
            RHS[cur_row] = P_img[v][dimx];
            cur_row++;

        }
        geo_assert(cur_row == 12);

        Matrix<12,double> inv = M.inverse();

        FOR(i, 12)            X[i] = 0;
        FOR(i, 12) FOR(j, 12) X[i] += inv(j, i)*RHS[j];

        vec3 res(0, 0, 0);
        FOR(dimx, 3)  FOR(dimu, 3)
            res[dimx] += in[dimu] * X[4 * dimx + dimu];

        FOR(dimx, 3) res[dimx] += X[4 * dimx + 3];
        return res;
    }

    bool PGPopt::in_tet(vec3 test, vec3 P[4], double eps) {
        int find[4][3] = { { 1, 3, 2 }, { 3, 0, 2 }, { 0, 3, 1 }, { 0, 1, 2 } };
        if (std::abs(dot(P[3] - P[0], cross(P[2] - P[0], P[1] - P[0]))) < 1e-15) return false;// check for flat tet
        FOR(f, 4) {
            vec3 vA = P[find[f][0]];
            vec3 vB = P[find[f][1]];
            vec3 vC = P[find[f][2]];
            vec3 n = cross(vB - vA, vC - vA);
            n = normalize(n);
            if (dot(test - vA, n) > eps) return false;
        }
        return true;
    }

    void PGPopt::get_grid_vertices(index_t t, std::vector<vec3>& psetX, std::vector<vec3>& psetU, bool dual) {
        //geo_assert(!tet_is_PGP_singular_fct(t));
        for (index_t i = 1; i < 4; i++) make_compatible(m->cells.vertex(t, 0), m->cells.vertex(t, i));
        // compute a bbox in the parametric domain
        int bbox[2][3] = { { 1000000, 1000000, 1000000 }, { -1000000, -1000000, -1000000 } };
        FOR(i, 4) FOR(dim, 3) {
            bbox[0][dim] = geo_min(bbox[0][dim], int(std::floor(U[m->cells.vertex(t, i)][dim])) - 1);
            bbox[1][dim] = geo_max(bbox[1][dim], int(std::ceil(U[m->cells.vertex(t, i)][dim])) + 1);
        }

        vec3 lX[4], lU[4];
        FOR(i, 4) {
            lX[i] = m->vertices.point(m->cells.vertex(t, i));
            lU[i] = U[m->cells.vertex(t, i)];
        }

        // raster the bbox and project to geometric space points included in the parametric space tet
        for (int i0 = bbox[0][0]; i0 <= bbox[1][0]; i0++) {
            for (int i1 = bbox[0][1]; i1 <= bbox[1][1]; i1++) {
                for (int i2 = bbox[0][2]; i2 <= bbox[1][2]; i2++) {
                    vec3 testpt(i0, i1, i2);
                    if (dual) testpt += vec3(.5, .5, .5);
                    if (in_tet(testpt, lU, .001)) {
                        psetX.push_back(change_tet_basis(testpt, lU, lX));
                        psetU.push_back(testpt);
                    }
                }
            }
        }
    }

    inline bool intersect_unit_box(vec3 p_boxcenter, vec3 tri[3]) {
        float boxcenter[3] = { float(p_boxcenter[0]), float(p_boxcenter[1]), float(p_boxcenter[2]) };
        float boxhalfsize[3] = { .5, .5, .5 };
        float triverts[3][3] = {
	    { float(tri[0][0]), float(tri[0][1]), float(tri[0][2]) },
	    { float(tri[1][0]), float(tri[1][1]), float(tri[1][2]) },
	    { float(tri[2][0]), float(tri[2][1]), float(tri[2][2]) }
        };
        return (triBoxOverlap(boxcenter, boxhalfsize, triverts) != 0);
    }

    inline void init_centered_unit_cube_face_bary(vec3 *cube_face_bary) {
        vec3 cubeU[8];
        FOR(k, 2)FOR(j, 2)FOR(i, 2)
            cubeU[4 * i + 2 * j + k] = vec3(i, j, k) - vec3(.5, .5, .5);
        CellDescriptor hexdescr = MeshCellsStore::cell_type_to_cell_descriptor(MESH_HEX);
        FOR(lf, 6) {
            cube_face_bary[lf] = vec3(0, 0, 0);
            for (int lv = 0; lv < 4; lv++) cube_face_bary[lf] += .5 * cubeU[hexdescr.facet_vertex[lf][lv]];
        }
    }


    void PGPopt::export_points(Mesh* pts) {
        vector<bool> PGP_singular(m->cells.nb(), false);
        FOR(c, m->cells.nb()) FOR(cf, 4) PGP_singular[c] = PGP_singular[c] || triangle_is_frame_singular(m, B, c, cf) || is_PGP_singular(c, cf);
        FOR(c, m->cells.nb()) {
            if (c % 100 == 0) GEO::Logger::out("HexDom")  << " EXPORT point set: tet" << c << " / " << m->cells.nb() <<  std::endl;
            if (PGP_singular[c]) continue;

            std::vector<vec3> seedU;
            std::vector<vec3> seedX;
            get_grid_vertices(c, seedX, seedU, false);
            index_t off = pts->vertices.create_vertices(index_t(seedX.size()));
            FOR(lv, seedX.size()) X(pts)[off + lv] = seedX[lv];
        }

        double eps = (1e-2)*get_cell_average_edge_size(m);
        vector<index_t> to_kill(pts->vertices.nb(), 0);
        vector<index_t> old2new(pts->vertices.nb());
        Geom::colocate(pts->vertices.point_ptr(0), 3, pts->vertices.nb(), old2new, eps);
        FOR(v, pts->vertices.nb()) if (old2new[v] != v) to_kill[v] = NOT_AN_ID;
        plop(pts->vertices.nb());
        pts->vertices.delete_elements(to_kill);
        plop(pts->vertices.nb());

    }


    void PGPopt::export_hexes(Mesh* hex) {
        Attribute<bool> touch_border(hex->cell_facets.attributes(), "tb");
        vec3 cube_face_bary[6];
        init_centered_unit_cube_face_bary(cube_face_bary);

        vector<bool> PGP_singular(m->cells.nb(), false);
        FOR(c, m->cells.nb()) FOR(cf, 4) PGP_singular[c] = PGP_singular[c] || triangle_is_frame_singular(m, B, c, cf) || is_PGP_singular(c, cf);


        FOR(c, m->cells.nb()) {

            if (c % 1000 == 0) GEO::Logger::out("HexDom")  << " EXPORT HEXES  tet = " << c << " / " << m->cells.nb() <<  std::endl;
            if (PGP_singular[c]) continue;
            std::vector<vec3> seedU; // contains centers of all cubes inside tet c
            {// init seedU
                std::vector<vec3> seedX;
                get_grid_vertices(c, seedX, seedU, true);
            }
            FOR(sid, seedU.size()) {
                vec3 ptsU[8]; // vertices of the cube centered in seedU[sid]
                vec3 ptsX[8];
                bool ptsdone[8] = {};
                bool have_boundary_face[6] = {};
                FOR(k, 2) FOR(j, 2) FOR(i, 2)
                    ptsU[4 * i + 2 * j + k] = seedU[sid] - vec3(.5, .5, .5) + vec3(i, j, k);


                std::vector<index_t> tet_stack;
                tet_stack.push_back(c);
                for (index_t i = 1; i < 4; i++) make_compatible(m->cells.vertex(c, 0), m->cells.vertex(c, i));


                std::vector<index_t> tet_done;
                bool has_sing_tet = false;

                while (!tet_stack.empty()) {
                    index_t curt = tet_stack.back(); tet_stack.pop_back();

                    // check if tet is singular
                    if (PGP_singular[curt]) { has_sing_tet = true; break; }


                    // create local geometry
                    vec3 lX[4], lU[4];
                    FOR(i, 4) {
                        lX[i] = m->vertices.point(m->cells.vertex(curt, i));
                        lU[i] = U[m->cells.vertex(curt, i)];
                    }

                    // find cubes corners located inside the current tet
                    FOR(i, 8) if (in_tet(ptsU[i], lU, 1e-5)) {
                        ptsX[i] = change_tet_basis(ptsU[i], lU, lX);
                        ptsdone[i] = true;
                    }

                    FOR(lf, 4) {
                        index_t verts[3] = {
                            m->cells.facet_vertex(curt, lf, 0),
                            m->cells.facet_vertex(curt, lf, 1),
                            m->cells.facet_vertex(curt, lf, 2)
                        };
                        vec3 tri[3] = { U[verts[0]], U[verts[1]], U[verts[2]] };

                        index_t oppt = m->cells.adjacent(curt, lf);
                        if (oppt == NO_CELL) {
                            // if the facet matches a unit cube facet in U coordinates, mark this unit cube facet
                            vec3 face_normal = normalize(cross(tri[1] - tri[0], tri[2] - tri[0]));
                            FOR(cubef, 6)
                                if ((face_normal - cube_face_bary[cubef]).length2() < 1e-10
                                    &&  round(.5 + dot(cube_face_bary[cubef], tri[0] - seedU[sid])) == 1) {
                                    have_boundary_face[cubef] = true;
                                    // need to check that it is not a concave hardedge
                                    FOR(cubef_in, 6) {
                                        if (cubef == cubef_in) continue;

                                        bool all_are_outside = true;
                                        FOR(d, 3) {
                                            all_are_outside = all_are_outside &&
                                                (.5 + dot(cube_face_bary[cubef_in], tri[d] - seedU[sid]) > 1 - 1e-5);
                                        }
                                        have_boundary_face[cubef] = have_boundary_face[cubef] && !all_are_outside;
                                    }
                                }
                            continue;
                        }

                        bool alreadydone = false;
                        FOR(p, tet_done.size()) alreadydone = alreadydone || (tet_done[p] == oppt);
                        if (alreadydone) continue;

                        if (intersect_unit_box(seedU[sid], tri)) {
                            FOR(i, 4) {
                                index_t v = m->cells.vertex(oppt, i);
                                if (v != verts[0] && v != verts[1] && v != verts[2]) {
                                    make_compatible(verts[0], v);
                                }
                            }
                            tet_stack.push_back(oppt);

                        }
                    }
                    tet_done.push_back(curt);
                }
                bool all_vertices_found = true;
                FOR(i, 8) all_vertices_found = all_vertices_found && ptsdone[i];
                if (all_vertices_found && !has_sing_tet) {
                    index_t off = hex->vertices.create_vertices(8);
                    index_t cid = hex->cells.create_hex(off, off + 1, off + 2, off + 3, off + 4, off + 5, off + 6, off + 7);
                    FOR(nb_f, 6)  touch_border[hex->cells.facet(cid, nb_f)] = have_boundary_face[nb_f]; //rand() % 2;
                    FOR(i, 8)      hex->vertices.point(off + i) = ptsX[i];
                }
            }
        }
    }



    void PGPopt::export_boundary_with_uv(Mesh* hex, Attribute<vec2>& uv, Attribute<index_t> &singtri) {
        // STEP 1: compute number of boundary vertices AND create a lookup table tet mesh vertex -> boundary surface vertex
        vector<index_t> tetV_to_facetV(m->vertices.nb(), NOT_AN_ID);
        index_t nb_boundary_V = 0;
        FOR(c, m->cells.nb()) FOR(cf, m->cells.nb_facets(c)) {
            if (NO_CELL != m->cells.adjacent(c, cf)) continue;
            FOR(cfv, m->cells.facet_nb_vertices(c, cf)) {
                index_t v = m->cells.facet_vertex(c, cf, cfv);
                if (NOT_AN_ID == tetV_to_facetV[v]) {
                    tetV_to_facetV[v] = nb_boundary_V++;
                }
            }
        }

        // STEP 2: copy boundary vertices to the new mesh
        hex->vertices.create_vertices(nb_boundary_V);
        FOR(v, m->vertices.nb())
            if (NOT_AN_ID != tetV_to_facetV[v])
                hex->vertices.point(tetV_to_facetV[v]) = m->vertices.point(v);

        // STEP 3: extract triangles and their parameterization (when possible)
        FOR(c, m->cells.nb()) FOR(cf, m->cells.nb_facets(c)) {
            if (m->cells.adjacent(c, cf) != NO_CELL) continue;
            index_t tet_verts[3];
            FOR(cfv, 3) {
                tet_verts[cfv] = m->cells.facet_vertex(c, cf, cfv);
            }

            for (index_t cfv = 1; cfv < 3; cfv++)
                make_compatible(tet_verts[0], tet_verts[cfv]);

            vec3 lX[3], lU[3];
            FOR(cfv, 3) {
                lX[cfv] = m->vertices.point(tet_verts[cfv]);
                lU[cfv] = U[tet_verts[cfv]];
            }

            // STEP 2.2: non singular boundary triangles are easy to extract
            index_t fid = hex->facets.create_triangle(
                tetV_to_facetV[tet_verts[0]],
                tetV_to_facetV[tet_verts[1]],
                tetV_to_facetV[tet_verts[2]]
            );

            bool has_valid_2d_param = false;

            if (!triangle_is_frame_singular(m, B, c, cf) && !is_PGP_singular(c, cf)) {
                FOR(dim, 3) { // we are looking for the param dimension that goes inside the volume
                    if (lU[0][dim] != lU[1][dim] || lU[0][dim] != lU[2][dim]) continue;
                    has_valid_2d_param = true;
                    FOR(lv, 3) {
                        uv[hex->facets.corner(fid, lv)] = vec2(lU[lv][(dim + 1) % 3], lU[lv][(dim + 2) % 3]);
                    }
                }
                // mirror the copy when necessary
                if (det(uv[hex->facets.corner(fid, 1)] - uv[hex->facets.corner(fid, 0)],
                    // TODO checker les inversions
                    uv[hex->facets.corner(fid, 2)] - uv[hex->facets.corner(fid, 0)]) < 0)
                    FOR(lv, 3) uv[hex->facets.corner(fid, lv)][0] *= -1.;
            }
            if (has_valid_2d_param) { // check param quality and invalidate it, if necessary
                TrglGradient grd(lX[0], lX[1], lX[2]);
                vec3 grduv[2];
                FOR(d, 2) grduv[d] = grd.gradient_3d(uv[hex->facets.corner(fid, 0)][d], uv[hex->facets.corner(fid, 1)][d], uv[hex->facets.corner(fid, 2)][d]);
                FOR(d, 2)FOR(dd, 3) if (GEO::Numeric::is_nan(grduv[d][dd])) has_valid_2d_param = false;
                if (grduv[0].length() > 10. * grduv[1].length()) has_valid_2d_param = false;
                if (grduv[1].length() > 10. * grduv[0].length()) has_valid_2d_param = false;
                if (std::abs(dot(normalize(grduv[0]), normalize(grduv[1]))) > cos(M_PI / 4.)) has_valid_2d_param = false;
            }

            if (!has_valid_2d_param) {
                FOR(lv, 3) {
                    uv[hex->facets.corner(fid, lv)] = vec2(.5, .5);
                }
            }

            singtri[fid] = !has_valid_2d_param;
        }
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

        //if (false) // disable for debug
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
    }
}
