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

#include <exploragram/hexdom/hex_cruncher.h>
#include <exploragram/hexdom/intersect_tools.h>
#include <exploragram/hexdom/polygon.h>
#include <exploragram/hexdom/mesh_inspector.h>
#include <geogram/points/colocate.h>
#include <geogram/mesh/triangle_intersection.h>
#include <geogram/mesh/mesh_repair.h> //used in bourrin subdivide hex
#define FPG_UNCERTAIN_VALUE 0
#include <geogram/numerics/predicates/orient3d.h>

namespace GEO {
    
    static bool intersect_quad(Mesh* m, HBoxes &hb, vector<vec3>& Q, double shrink = 0) {

        vec3 G = Poly3d(Q).barycenter();
        FOR(v, 4)Q[v] = (1. - shrink)*Q[v] + shrink*G;
        vec3 n = Poly3d(Q).normal();
        double decal = 0;
        FOR(cfv, 4) decal += .2*.25*(Q[cfv] - Q[next_mod(cfv, 4)]).length(); // .2*ave edge length
        Q.push_back(G + decal*n);
        Q.push_back(G - decal*n);

        BBox b;
        FOR(cfv, 6) b.add(Q[cfv]);

        vector<index_t> primitives;
        hb.intersect(b, primitives);

        FOR(i, primitives.size()) {
            index_t f = primitives[i];
            vector<vec3> P;
            FOR(fv, m->facets.nb_vertices(f))
                P.push_back(X(m)[m->facets.vertex(f, fv)]);
            geo_assert(P.size() < 5);

            vector<TriangleIsect> trash;
            FOR(diam, 12) {
                if (P.size() == 3) {
                    if (triangles_intersections(P[0], P[1], P[2],
                        Q[diamon_split[diam][0]], Q[diamon_split[diam][1]], Q[diamon_split[diam][2]], trash)) return true;
                }
                else {
                    FOR(qu, 4)
                        if (triangles_intersections(P[quad_split[qu][0]], P[quad_split[qu][1]], P[quad_split[qu][2]],
                            Q[diamon_split[diam][0]], Q[diamon_split[diam][1]], Q[diamon_split[diam][2]], trash)) return true;
                }
            }
        }
        return false;
    }


    inline double tetra_volume(vec3 A, vec3 B, vec3 C, vec3 D) {
        return dot(cross(B - A, C - A), D - A);
    }
    
    inline double tetra_volume_sign(vec3 A, vec3 B, vec3 C, vec3 D) {
        double res = tetra_volume(A, B, C, D);
        if (std::abs(res) > 1e-15) return res;
        return dot(normalize(cross(normalize(B - A), normalize(C - A))), normalize(D - A));
    }

    inline bool same_sign(double a, double b) { return (a>0) == (b > 0); }

    static bool canonical_order_set_A(vec3 A_ref, vec3 B_ref, vec3 C_ref, vec3& A, vec3& B, vec3& C) {
        // h* gives the halfspace delimited by triangle Aref,Bref,Cref
        double ha, hb, hc;
        vec3 n = cross(B_ref - A_ref, C_ref - A_ref);
        ha = dot(n, A - A_ref);
        hb = dot(n, B - A_ref);
        hc = dot(n, C - A_ref);
        if (std::abs(ha) < 1e-15 || std::abs(hb) < 1e-15 || std::abs(hc) < 1e-15) {
            ha = orient_3d_filter(A_ref.data(), B_ref.data(), C_ref.data(), A.data());
            hb = orient_3d_filter(A_ref.data(), B_ref.data(), C_ref.data(), B.data());
            hc = orient_3d_filter(A_ref.data(), B_ref.data(), C_ref.data(), C.data());
        }
        // avoid that B is alone on his side
        if (!same_sign(hb, hc) && !same_sign(hb, ha)) {
            std::swap(A, B);
            std::swap(ha, hb);
        }
        // if all points are on the same side there is no intersection
        if (same_sign(ha, hc)) return false;
        // set A to be the vertex that is alone on his side
        if (same_sign(hb, ha)) std::swap(A, C);
        return true;
    }


    static bool naive_tri_tri_intersect(vec3 v0, vec3 v1, vec3 v2, vec3 u0, vec3 u1, vec3 u2) {

        // the support planes of each triangles produces 2 halfspace.
        // vertices order is changed such that v0 and u0 are alone in the halfspace of the other triangle
        // if it is impossible, there is no intersection
        if (!canonical_order_set_A(v0, v1, v2, u0, u1, u2)) return false;
        if (!canonical_order_set_A(u0, u1, u2, v0, v1, v2))	return false;

        if (orient_3d_filter(v0.data(), v1.data(), v2.data(), u0.data()) > 0)   std::swap(v1, v2);
        if (orient_3d_filter(u0.data(), u1.data(), u2.data(), v0.data()) > 0)   std::swap(u1, u2);

        if (orient_3d_filter(u0.data(), v1.data(), v0.data(), u1.data())>=0)    return false;
        if (orient_3d_filter(u0.data(), v2.data(), v0.data(), u2.data())<= 0)  return false;

        return true;
    }

    /*
    static bool intersects(Mesh* quadtri, HBoxes &hb, vector<index_t> Qi) {
        vector<vec3> Q;
        FOR(i, 4) Q.push_back(X(quadtri)[Qi[i]]);

        vec3 G = Poly3d(Q).barycenter();
        vec3 n = Poly3d(Q).normal();
        double decal = 0;
        FOR(cfv, 4) decal += .2*.25*(Q[cfv] - Q[next_mod(cfv, 4)]).length(); // .2*ave edge length
        Q.push_back(G + decal*n);
        Q.push_back(G - decal*n);

        BBox b;
        FOR(cfv, 6) b.add(Q[cfv]);

        vector<index_t> primitives;
        hb.intersect(b, primitives);

        FOR(i, primitives.size()) {
            index_t f = primitives[i];
            vector<vec3> P;
            FOR(fv, quadtri->facets.nb_vertices(f))
                P.push_back(X(quadtri)[quadtri->facets.vertex(f, fv)]);
            geo_assert(P.size() < 5);

            if (4 == P.size()) { // check direct coincidence without diamonding
                bool have_same_vertices = true;
                FOR(fv1, 4) {
                    bool isin = false;
                    FOR(fv2, 4) {
                        isin = isin || quadtri->facets.vertex(primitives[i], fv1) == Qi[fv2];
                    }
                    have_same_vertices = have_same_vertices && isin;
                }
                if (have_same_vertices) continue;
            }

            vector<TriangleIsect> result;
            FOR(diam, 12) {
                if (P.size() == 3) {
                    if (triangles_intersections(P[0], P[1], P[2],
                        Q[diamon_split[diam][0]], Q[diamon_split[diam][1]], Q[diamon_split[diam][2]], result)) return true;
                }
                else {
                    FOR(qu, 4)
                        if (triangles_intersections(P[quad_split[qu][0]], P[quad_split[qu][1]], P[quad_split[qu][2]],
                            Q[diamon_split[diam][0]], Q[diamon_split[diam][1]], Q[diamon_split[diam][2]], result)) return true;
                }
            }
        }
        return false;
    }
    */

    struct QuadExtraConnectivity {

        void init(Mesh* p_m) {
            m = p_m;
            FOR(f, m->facets.nb()) geo_assert(m->facets.nb_vertices(f) == 4);// check that surface is quadragulated
            opp_h.clear();
            opp_h.resize(4 * m->facets.nb(), NOT_AN_ID); // NOT facet_corners !!!
            create_non_manifold_facet_adjacence(m);
            FOR(f, m->facets.nb()) FOR(lc, 4) {
                index_t opp_f = m->facets.adjacent(f, lc);
                if (opp_f == NOT_AN_ID) { plop("bad adjacency detected "); continue; }

                index_t opp_lc = NOT_AN_ID;
                FOR(opp_lc_it, 4) if (f == m->facets.adjacent(opp_f, opp_lc_it)
                    && m->facets.vertex(f, (lc + 1) % 4) == m->facets.vertex(opp_f, opp_lc_it)) opp_lc = opp_lc_it;

                if (opp_lc == NOT_AN_ID) plop("not a symetric opposite !");

                set_opp(4 * f + lc, 4 * opp_f + opp_lc);
                if (vertex(4 * f + lc) != vertex(next(opp(4 * f + lc)))) plop("wrong opposites");
                if (vertex(next(4 * f + lc)) != vertex(opp(4 * f + lc))) plop("wrong opposites");
            }
        }

        void debug_export_adjacence() {
            FOR(f, m->facets.nb()) FOR(lc, 4) {
                if (opp(4 * f + lc)!=NOT_AN_ID) 
                    m->facets.set_adjacent(f, lc, face(opp(4 * f + lc)));
                else m->facets.set_adjacent(f, lc, NOT_AN_ID);
            }
        }

        bool valid(index_t h) { return h < 4 * m->facets.nb(); }

        void set_opp(index_t i, index_t j) { geo_assert(valid(i) && valid(j));  opp_h[i] = j; opp_h[j] = i; }
        index_t face(index_t e) { geo_assert(valid(e)); return e / 4; }
        index_t local_id(index_t e) { geo_assert(valid(e)); return e % 4; }
        index_t next(index_t e, index_t nb = 1) { geo_assert(valid(e)); return 4 * face(e) + ((e + nb) % 4); }
        index_t opp(index_t e) { geo_assert(valid(e)); return opp_h[e]; }
        index_t vertex(index_t e) { geo_assert(valid(e)); return m->facets.vertex(face(e), local_id(e)); }
        index_t next_around_vertex(index_t e) { geo_assert(valid(e)); return opp(next(e, 3)); }

        void set_vertex(index_t e, index_t v) { geo_assert(valid(e) && v<m->vertices.nb()); m->facets.set_vertex(face(e), local_id(e), v); }


        bool has_valid_one_ring(index_t e) {
            index_t cir = e;
            int count = 0;
            do {
                if (cir == NOT_AN_ID)   return false;
                if (count++ == 1000) return false;
                cir = next_around_vertex(cir);
            } while (cir != e);
            return true;
        }
        int valence(index_t e) {
            geo_assert(has_valid_one_ring(e));
            index_t cir = e;
            int count = 0;
            do {
                count++;
                cir = next_around_vertex(cir);
            } while (cir != e);
            return count;
        }

        bool is_closed() {
            FOR(h, opp_h.size()) if (opp_h[h] == NOT_AN_ID) {
                Attribute<double> deb(m->vertices.attributes(), "debug");
                deb[vertex(h)] = 10;
                deb[vertex(next(h))] = 10;
                return false;
            }
            FOR(v, m->vertices.nb()) if (!has_valid_one_ring(v)) {
                Attribute<double> deb(m->vertices.attributes(), "debug");
                deb[v] = 10;
                return false;
            }
            return true;
        }


        Mesh* m;
        vector<index_t> opp_h;
    };



    struct CutSingularity {
        Mesh* m;                // ;)
        vector<bool> visited;   // prevents iterating more than once on the same cut
        vec3 N;                 // normal to the current cut
        QuadExtraConnectivity qem;
        vector<index_t> border;
        vector<vec3> pts;
        vector<index_t> quads;

        CutSingularity(Mesh* p_m) { 
            m = p_m; 
            qem.init(m);
            visited.resize(4*m->facets.nb(), false);
        }

        bool  create_edge_loop(index_t  h, vector<index_t>& test_border) {
            Attribute<bool> isquad(m->facets.attributes(), "isquad");
            double sigma_angle = 0;
            index_t cir = h;
            do {
                visited[cir] = true;
                vec3 Nup = facet_normal(m, qem.face(cir));
                if (qem.opp(cir) == NOT_AN_ID) { test_border.clear(); break; }
                vec3 Ndown = facet_normal(m, qem.face(qem.opp(cir)));
                if (dot(N, Nup) > .5 || dot(N, Ndown) < -.5) { test_border.clear(); break; }


                vec3 cir_dir = X(m)[qem.vertex(qem.next(cir))] - X(m)[qem.vertex(cir)];
                cir_dir = normalize(cir_dir);

                test_border.push_back(cir);
                index_t next = NOT_AN_ID;
                index_t in_cir = qem.next(cir);


                //double best_dot = -1e20;
                double best_angle = -M_PI;
                do {
                    if (!isquad[qem.face(in_cir)] || qem.opp(in_cir) == NOT_AN_ID) {
                        next = NOT_AN_ID;
                        break;
                    }

                    vec3 in_cir_dir = X(m)[qem.vertex(qem.next(in_cir))] - X(m)[qem.vertex(in_cir)];
                    in_cir_dir = normalize(in_cir_dir);

                    if (fabs(dot(in_cir_dir, N)) < sin(M_PI / 8.)	    // stay in the cut plane (orthogonal to N)
                        && in_cir != qem.opp(cir)				        // do not go back
                        ) {
                        double newangle = atan2(dot(N, cross(in_cir_dir, cir_dir)), dot(in_cir_dir, cir_dir));
                        if (best_angle < newangle) {
                            best_angle = newangle;
                            next = in_cir;
                        }
                    }
                    in_cir = qem.next(qem.opp(in_cir));
                } while (in_cir != qem.next(cir));
                sigma_angle += best_angle;

                if (next == NOT_AN_ID || test_border.size() > 30) {
                    test_border.clear();
                    break;
                }
                cir = next;
            } while (cir != h);
            if (sigma_angle > 0) return false;
            return true;
        }

        bool cut_separates_2_shorts_closed_quads_strips(index_t h) {
            // check that it separates two smalls and close quads strips
            index_t quad_strip_start[2] = { qem.next(h), qem.next(qem.opp(h)) };
            FOR(q, 2) {
                index_t it = quad_strip_start[q];
		int i=0;
		for(;;) {
                    if (i == 30 || qem.opp(it) == NOT_AN_ID)
                        return false;
                    it = qem.next(qem.opp(it), 2);
                    if (it == quad_strip_start[q]) break;
		    i++;
                }
            }
            return true;
        }


            bool can_easily_discard_current_edge_loop(vector<index_t>& test_border) {
                if (test_border.size() < 3) return true;
                if (!border.empty() && test_border.size() >= border.size()) return true;
                if (test_border.size() == 4
                    && (qem.face(qem.opp(test_border[0])) == qem.face(qem.opp(test_border[2]))
                        || qem.face(test_border[0]) == qem.face(test_border[2])
                        )
                    ) return true;
                if (test_border.size() % 2 != 0) return true;
                return false;
            }



        bool apply() {
            Attribute<bool> isquad(m->facets.attributes(), "isquad");
            vec3 best_N;
            double best_cut_area = 1e20;
            
            double ave_edge_length = 0;
            FOR(f, m->facets.nb()) FOR(v, 4) ave_edge_length += (X(m)[m->facets.vertex(f, v)] - X(m)[m->facets.vertex(f, (v + 1) % 4)]).length();
            ave_edge_length /= 4.*double(m->facets.nb());
            
            // STEP 1: determine a valid cut along halfedges "border", its quadrangulation "quads", with vertices "pts[i]" (starting with vertices of "border")
           

            int nb_border_tried = 0;
            int nb_border_success = 0;
            FOR(h, 4*m->facets.nb()) {
                if (visited[h]) continue;
                index_t f = qem.face(h);
                if (!isquad[f]) continue;
                if (qem.opp(h) == NOT_AN_ID) continue;
                if (!cut_separates_2_shorts_closed_quads_strips(h)) continue;
                N = X(m)[qem.vertex(qem.next(h,3))] - X(m)[qem.vertex(h)];
                N = normalize(N);

                vector<index_t> test_border;

                // STEP 1.1 create the edge loop starting from h
                if (!create_edge_loop(h, test_border)) continue;


                // STEP 1.2 check if the edge loop starting from h is valid and better than previous loop
                if (can_easily_discard_current_edge_loop(test_border)) continue;
                
                vector<vec3> test_pts;
                vector<index_t> test_quads;
                FOR(e, test_border.size()) test_pts.push_back(X(m)[qem.vertex(test_border[e])]);


                double cut_area = 0;
                FOR(e, test_pts.size()) cut_area -= dot(cross(test_pts[e], test_pts[(e + 1) % test_pts.size()]), N);
                if (best_cut_area < cut_area)continue;



                Poly3d p3(test_pts);
                plop("\n-------------------try_quadrangulate------------------------\n");
                nb_border_tried++;



                if (!p3.try_quadrangulate(test_quads)) {
                    plop("FAIL ");
                    continue;
                }
                else {
                    nb_border_success++;
                    bool cut_will_intersect = false;
                    vector<BBox> inboxes(m->facets.nb());
                    FOR(ff, m->facets.nb()) FOR(fv, m->facets.nb_vertices(ff)) inboxes[ff].add(X(m)[m->facets.vertex(ff, fv)]);
                    HBoxes hb(inboxes);
                    FOR(q, test_quads.size() / 4) {
                        vector<vec3> quad;
                        FOR(lc, 4) quad.push_back(test_pts[test_quads[4 * q + lc]]);
                        cut_will_intersect = cut_will_intersect || intersect_quad(m, hb, quad, 1e-5);
                    }
                    if (cut_will_intersect) {
                        plop("cannot split due to intersections");
                        continue;
                    }
                }

                // STEP 1.3 validate the new loop
                {
                    border.swap(test_border);
                    pts.swap(test_pts);
                    test_quads.swap(quads);
                    best_cut_area = cut_area;
                    best_N = N;
                    //plop("\n-------------------last border is current winner------------------------\n");
                }
            }

            plop(nb_border_tried);
            plop(nb_border_success);

            if (border.empty()) return false;


            vector<index_t> upper_v;
            vector<index_t> lower_v;
            {
                index_t off_v = m->vertices.create_vertices(border.size());
                FOR(e, border.size()) {
                    upper_v.push_back(qem.vertex(border[e]));
                    pts.push_back(X(m)[upper_v[e]]);
                    lower_v.push_back(off_v + e);
                    X(m)[off_v + e] = pts[e];
                };
            }
            vector<vector<index_t> > opp_fan(border.size());
            FOR(e, border.size()) {
                index_t cir = qem.opp(border[e]);
                do {
                    opp_fan[e].push_back(cir);
                    cir = qem.opp(qem.next(cir,3));
                } while (cir != border[next_mod(e, border.size())]);
            }
            FOR(e, border.size()) FOR(v, opp_fan[e].size()) m->facet_corners.set_vertex(opp_fan[e][v], lower_v[next_mod(e, border.size())]);


            if (pts.size() > border.size()) {
                index_t off_v = m->vertices.create_vertices(2 * (pts.size() - border.size()));
                FOR(i, pts.size() - border.size()) {
                    X(m)[off_v + 2 * i] = pts[border.size() + i];
                    X(m)[off_v + 2 * i + 1] = pts[border.size() + i];
                    upper_v.push_back(off_v + 2 * i);
                    lower_v.push_back(off_v + 2 * i + 1);
                }
            }
            FOR(v, lower_v.size()) X(m)[lower_v[v]] = X(m)[lower_v[v]] - .1*ave_edge_length *best_N;
            FOR(q, quads.size() / 4) {
                isquad[m->facets.create_quad(
                    upper_v[quads[4 * q + 0]],
                    upper_v[quads[4 * q + 3]],
                    upper_v[quads[4 * q + 2]],
                    upper_v[quads[4 * q + 1]]
                )] = true;
                isquad[m->facets.create_quad(
                    lower_v[quads[4 * q + 0]],
                    lower_v[quads[4 * q + 1]],
                    lower_v[quads[4 * q + 2]],
                    lower_v[quads[4 * q + 3]]
                )] = true;
            }
            // debug output

            Attribute<int> date(m->edges.attributes(), "date");
            if (!border.empty()) {
                //m->edges.clear();
                index_t off_e = m->edges.create_edges(border.size());
                FOR(e, border.size()) {
                    date[off_e + e] = int(off_e);
                    FOR(extr, 2)
                        m->edges.set_vertex(off_e + e, extr, qem.vertex(border[(e + extr) % border.size()]));
                }
            }
            return true;

        }
    };




    inline double cos_corner(vec3 B, vec3 A, vec3 C) {
        return dot(normalize(B - A), normalize(C - A));
    }

    struct VertexPuncher {
        Mesh* m;
        Mesh* newhex;
        QuadExtraConnectivity qem;
        DynamicHBoxes hb;                              //      -> a static BBox tree
        double ave_edge_length;
        Attribute<bool> dead_face;

        int nb_punchs;
        index_t punch_v;
        index_t nv_punch_v;
        index_t H[3][4];
        index_t oppH[3][4];
        vec3 old_vertex_position;
        vec3 new_vertex_position;
        vector<int> v2nb_facets;           //      -> facets that have moved
        index_t via_facet[3] ;


        Attribute<double> failt; ///DEBUG
        Attribute<bool> isquad;


        VertexPuncher(Mesh* p_m, Mesh* p_newhex)  {
            m = p_m;
            newhex = p_newhex;
            dead_face.bind(m->facets.attributes(), "dead_face");
        }



        void unglue_duplicates() {
            FOR(lf, 3) {
                index_t f = via_facet[lf];
                if (f == NOT_AN_ID) continue;
                index_t cir_h[4];
                index_t cir_opp[4];
                index_t opp_f= NOT_AN_ID;
                FOR(lh, 4) {
                    index_t h = 4 * f + lh;
                    index_t h_opp = qem.opp(h);
                    if (qem.vertex(qem.next(h_opp, 2)) == qem.vertex(qem.next(h, 3))
                        && qem.vertex(qem.next(h_opp, 3)) == qem.vertex(qem.next(h, 2))) {
                        FOR(i, 4) {
                            cir_h[i] = 4 * f + (lh + i) % 4;
                            cir_opp[i] = 4 * qem.face(h_opp) + (qem.local_id(h_opp) + 4 - i) % 4;
                        }
                        opp_f = qem.face(h_opp);
                        continue;
                    }
                }
                if (opp_f == NOT_AN_ID) continue;

                dead_face[f] = true;
                dead_face[opp_f] = true;
                FOR(i, 4) {
                    //std::cerr << qem.vertex(cir_h[i]) << "  " << qem.vertex(qem.next(cir_opp[i])) << "  \n";
                    if (qem.opp(cir_h[i]) == qem.opp(cir_opp[i])) continue;
                    qem.set_opp(qem.opp(cir_opp[i]), qem.opp(cir_h[i]));
                    qem.set_opp(cir_h[i], cir_opp[i]);
                };
            }
        }




        BBox facet_bbox(index_t f) {
            BBox res;
            FOR(fv, m->facets.nb_vertices(f)) res.add(X(m)[m->facets.vertex(f, fv)]);
            return res;
        }

        void init() {
            isquad.bind(m->facets.attributes(), "isquad");

            qem.init(m);
            nb_punchs = 0;
            v2nb_facets.clear();
            //moved_facets.clear();

            v2nb_facets.resize(m->vertices.nb(), 0);
            FOR(f, m->facets.nb()) FOR(v, 4) v2nb_facets[m->facets.vertex(f, v)]++;



            FOR(f, m->facets.nb())dead_face[f] = false;

            // mesh resolution
            ave_edge_length = 0;
            FOR(f, m->facets.nb()) FOR(v, 4) ave_edge_length += (X(m)[m->facets.vertex(f, v)] - X(m)[m->facets.vertex(f, (v + 1) % 4)]).length();
            ave_edge_length /= 4.*double(m->facets.nb());

            // structures to find facets

            vector<BBox> inboxes(m->facets.nb());
            FOR(f, m->facets.nb()) inboxes[f] = facet_bbox(f);//FOR(fv, m->facets.nb_vertices(f)) inboxes[f].add(X(m)[m->facets.vertex(f, fv)]);
            hb.init(inboxes);


        }



        bool init_one_ring(index_t h) {
            FOR(f, 3) {
                FOR(v, 4) {
                    H[f][v] = qem.next(h, v);
                    if (H[f][v] == NOT_AN_ID)   return false;// TO REMOVE if m is closed                    
                }
                if (!isquad[qem.face(H[f][0])])   return false;
                h = qem.next_around_vertex(h);
            }
            if (qem.next_around_vertex(H[2][0]) != H[0][0])
                return false;

            FOR(f, 3) FOR(e, 4) {
                oppH[f][e] = qem.opp(H[f][e]);
                if (oppH[f][e] == NOT_AN_ID)       return false;// TO REMOVE if m is closed
                
            }
            return true;
        }
        // create new hex
        void create_new_hex() {
            index_t off_v = newhex->vertices.create_vertices(8);
            Attribute<double> init(newhex->vertices.attributes(), "init");
            init[off_v] = -100; FOR(i, 7)init[off_v + i + 1] = nb_punchs;
            X(newhex)[off_v + 0] = old_vertex_position;
            X(newhex)[off_v + 1] = X(m)[qem.vertex(H[0][1])];
            X(newhex)[off_v + 2] = X(m)[qem.vertex(H[0][3])];
            X(newhex)[off_v + 3] = X(m)[qem.vertex(H[0][2])];
            X(newhex)[off_v + 4] = X(m)[qem.vertex(H[2][1])];
            X(newhex)[off_v + 5] = X(m)[qem.vertex(H[2][2])];
            X(newhex)[off_v + 6] = X(m)[qem.vertex(H[1][2])];
            X(newhex)[off_v + 7] = X(m)[nv_punch_v];
            newhex->cells.create_hex(off_v + 0, off_v + 1, off_v + 2, off_v + 3, off_v + 4, off_v + 5, off_v + 6, off_v + 7);
        }

        bool new_hex_geometry_is_crappy() {
           // return false;
            FOR(front, 2) {// front=0 for actual faces, front =1 for new faces 
                FOR(fid, 3) {
                    vector<vec3> v(4);
                    if (front ==0) {
                        v[0] = old_vertex_position;
                        v[1] = X(m)[qem.vertex(H[fid][1])];
                        v[2] = X(m)[qem.vertex(H[fid][2])];
                        v[3] = X(m)[qem.vertex(H[fid][3])];
                    } else {
                        v[0] = X(m)[nv_punch_v];
                        v[1] = X(m)[qem.vertex(H[next_mod(fid, 3)][2])];
                        v[2] = X(m)[qem.vertex(H[next_mod(fid,3)][1])];
                        v[3] = X(m)[qem.vertex(H[fid][2])];
                    }
                    vec3 n = Poly3d(v).normal();
                    //vec3 bary = Poly3d(v).barycenter();
                    FOR(lv,4) if (dot(n,cross(normalize(v[(lv + 2) % 4] - v[(lv + 1) % 4]), normalize(v[(lv ) % 4] - v[(lv + 1) % 4]))) < .1)return true;

                }
            }
            return false;

        }

        //topo_punch
        void topo_punch() {
            FOR(f, 3) qem.set_vertex(H[f][0], nv_punch_v);
            FOR(f, 3) qem.set_vertex(H[f][1], qem.vertex(oppH[(f + 1) % 3][1]));
            FOR(f, 3) qem.set_vertex(H[f][2], qem.vertex(oppH[(f + 1) % 3][2]));
            FOR(f, 3) qem.set_vertex(H[f][3], qem.vertex(oppH[(f + 2) % 3][1]));
            FOR(f, 3) qem.set_opp(H[f][1], oppH[(f + 1) % 3][2]);
            FOR(f, 3) qem.set_opp(H[f][2], oppH[(f + 2) % 3][1]);
        }
        bool topo_can_punch() {
            vector<index_t> neigh;
            // check if a vertex is duplicated
            FOR(ring, 3) FOR(lv, 2) neigh.push_back(H[ring][lv + 1]);
            FOR(n0, 6) if (qem.opp(neigh[n0]) == NOT_AN_ID) return false;
            FOR(n0, 6)for (index_t n1 = 0; n1 < n0; n1++)
                if (qem.vertex(neigh[n0]) == qem.vertex(neigh[n1])) return false;;


            // check is two boundary edges are opposite
            FOR(n0, 6)for (index_t n1 = 0; n1 < n0; n1++)
                if (qem.opp(neigh[n0]) == neigh[n1]) return false;
            return true;
        }


        void produce_diamon(vector<vec3>& P, vec3 A, vec3 B, vec3 C, vec3 D, double h) {
            P.reserve(6);
            P.resize(4);
            P[0] = A; P[1] = B; P[2] = C; P[3] = D;
            vec3 G = Poly3d(P).barycenter();
            vec3 n = Poly3d(P).normal();
            P.push_back(G + h*n);
            P.push_back(G - h*n);
        }
        

        bool punch_will_produce_intersection() {
            // check for geometric intersections
            FOR(fid, 3) {

                // construct a little structure including the new face
                vector<vec3> Q;

                index_t Qv[4] = { nv_punch_v, qem.vertex(H[fid][2]),qem.vertex(H[next_mod(fid, 3)][1]),qem.vertex(H[next_mod(fid, 3)][2]) };
                produce_diamon(Q, X(m)[Qv[0]], X(m)[Qv[1]], X(m)[Qv[2]], X(m)[Qv[3]], .1*ave_edge_length);


                // get candidate quads for intersecting current new face
                BBox b;
                FOR(cfv, 6) b.add(Q[cfv]);
                
                // scale B to include facets that have vertex too close, but do not necessary intersect
                b.max = b.max + ave_edge_length * vec3(1, 1, 1);
                b.min = b.min - ave_edge_length * vec3(1, 1, 1);

                vector<index_t> primitives;
                hb.intersect(b, primitives); 
                //FOR(i, m->facets.nb()) primitives.push_back(i);

                FOR(primid, primitives.size()) {
                    index_t f = primitives[primid];
                    if (f == qem.face(H[0][0]) || f == qem.face(H[1][0]) || f == qem.face(H[2][0])) continue;
                    if (dead_face[f]) continue;

                    vector<vec3> P;

                    index_t Pv[4] = { m->facets.vertex(f, 0),m->facets.vertex(f, 1),m->facets.vertex(f, 2),m->facets.vertex(f, 3) };
                    produce_diamon(P, X(m)[Pv[0]], X(m)[Pv[1]], X(m)[Pv[2]], X(m)[Pv[3]], .1*ave_edge_length);




                    // check that the new vertex is far enough to other pts
                    FOR(fv, 4) if (nv_punch_v == punch_v
                        && m->facets.vertex(f, fv) != nv_punch_v
                        && (P[fv] - X(m)[nv_punch_v]).length() < .25*ave_edge_length)
                    {
                        failt[punch_v] = geo_max(failt[punch_v], 70.);
                        return true;
                    }
                    
                    // shrink shared elements, grow others
                    bool shrink_Q[4] = { false, false, false, false };
                    bool shrink_P[4] = { false, false, false, false };
                    vec3 GQ = .25*(Q[0] + Q[1] + Q[2] + Q[3]);
                    vec3 GP = .25*(P[0] + P[1] + P[2] + P[3]);
                    FOR(q, 4) FOR(p, 4)  if (Qv[q] == Pv[p]) { shrink_Q[q] = true; shrink_P[p] = true; }
                    FOR(i, 4) {
                        if (shrink_Q[i]) Q[i] = .99 * Q[i] + .01 * GQ;
                        if (shrink_P[i]) P[i] = .99 * P[i] + .01 * GP;
                    }

                    // precise bbox filter
                    BBox ba; FOR(pv, 6) ba.add(P[pv]);
                    BBox bb; FOR(pv, 6) bb.add(Q[pv]);
                    if (!ba.intersect(bb)) continue;

                    // check for triangle / triangle intersections
                    FOR(diam, 16) FOR(qu, 16)
                        if (naive_tri_tri_intersect(P[diamon_quad_split[qu][0]], P[diamon_quad_split[qu][1]], P[diamon_quad_split[qu][2]],
                                Q[diamon_quad_split[diam][0]], Q[diamon_quad_split[diam][1]], Q[diamon_quad_split[diam][2]])){
                            index_t off = newhex->vertices.create_vertices(2);
                            X(newhex)[off] = Poly3d(Q).barycenter();
                            X(newhex)[off+1] = Poly3d(P).barycenter();
                            newhex->edges.create_edge(off, off + 1);
                            {
                                failt[punch_v] = geo_max(failt[punch_v], 80. );
                                return true;
                            }
                       }

                    //vector<TriangleIsect> result;
                    //FOR(diam, 16) FOR(qu, 4)
                    //    if (triangles_intersections(P[quad_split[qu][0]], P[quad_split[qu][1]], P[quad_split[qu][2]],
                    //            Q[diamon_quad_split[diam][0]], Q[diamon_quad_split[diam][1]], Q[diamon_quad_split[diam][2]], result))
                    //        return true;
                }
            }
            return false;
        }


        bool apply(int& nbmaxpunch ) {
            if (m->facets.nb() == 0) return false;
            init();
            plop("init done");
            bool finished = false;
            failt.bind(m->vertices.attributes(), "failt");
            FOR(v, m->vertices.nb()) failt[v] = 0;

            
            while (!finished) {
                //plop("new loop"); plop(nb_punchs);
                bool found_vertex_to_punch = false;
                FOR(seed, 4 * m->facets.nb()) {
                    //if (!qem.is_closed()) { plop(seed); plop("!qem.is_closed(), look at debug attrib"); return false; }
                    if (qem.next_around_vertex(seed) < seed || qem.next_around_vertex(qem.next_around_vertex(seed)) < seed) continue;
                    punch_v = qem.vertex(seed);
                    if (!init_one_ring(seed))   { failt[punch_v] = geo_max(failt[punch_v],10.); continue; }
                    if (!topo_can_punch())      { failt[punch_v] = geo_max(failt[punch_v], 20.); continue; }

                    old_vertex_position = X(m)[punch_v];

                    nv_punch_v = punch_v;
                    // do we already have 4+ faces of the hex ?
                    bool bad_config = false;
                    index_t punch_cand [3]= { NOT_AN_ID, NOT_AN_ID, NOT_AN_ID };
                    FOR(i, 3) via_facet[i] =NOT_AN_ID;
                    index_t punch_cand_ref= NOT_AN_ID;
                    FOR(i, 3) if (qem.valence(oppH[i][2]) == 3) {
                        punch_cand[i] = qem.vertex(qem.next(oppH[i][2], 2));
                        punch_cand_ref = punch_cand[i];
                        via_facet[i] = qem.face(oppH[i][2]);
                    }

                    FOR(i, 3) {
                        if (punch_cand[i] != NOT_AN_ID && punch_cand[i] != punch_cand_ref)  
                            bad_config = true;
                        if (qem.valence(H[i][2]) == 3)
                            if ((via_facet[i] == NOT_AN_ID) != (via_facet[(i+2)%3] == NOT_AN_ID))
                                bad_config = true;
                        if (via_facet[i] != NOT_AN_ID && !isquad[via_facet[i]])    
                            bad_config = true;
                    }
                    if (bad_config) { failt[punch_v] = geo_max(failt[punch_v], 20.); continue; }

                    if (punch_cand_ref != NOT_AN_ID){
                        nv_punch_v = punch_cand_ref;
                        FOR(i, 3) if (via_facet[i] != NOT_AN_ID) dead_face[via_facet[i]] = true;             // the face will be cancelled by its opposite
                    }
                    


                    if (nv_punch_v == punch_v) {
                        //if (search_for_existing_vertex_to_punch) continue;
                        // check geometry of existing faces
                        bool have_bad_angle = false;
                        FOR(ring, 3) FOR(lv, 4)
                            have_bad_angle = have_bad_angle || std::abs(cos_corner(
                                X(m)[qem.vertex(H[ring][lv])],
                                X(m)[qem.vertex(H[ring][(lv + 1) % 4])],
                                X(m)[qem.vertex(H[ring][(lv + 2) % 4])]
                            )) > .8;
                        if (have_bad_angle) { failt[punch_v] = geo_max(failt[punch_v], 30.); continue; }

                        //-------------------------
                        {
                            vec3 cubebary(0, 0, 0);
                            FOR(i, 3) FOR(e, 2) cubebary = cubebary + X(m)[qem.vertex(H[i][1 + e])];
                            cubebary = (1. / 6.) *cubebary;
                            vec3 decal = cubebary - old_vertex_position;
                            vec3 n[3];
                            FOR(i, 3) n[i] = facet_normal(m, qem.face(H[i][0]));
                            if (dot(decal, n[0] + n[1] + n[2]) > 0) { 
                                failt[punch_v] = geo_max(failt[punch_v], 40.);
                                continue; 
                            }
                            new_vertex_position = 2. * cubebary - old_vertex_position;
                        }
                        // new way to compute new position
                        vec3 n[3];
                        mat3 mat;
                        index_t tri[3][3];// 3 triplet of vertices that miss the last point
                        FOR(f, 3) {
                            index_t next_f = next_mod(f, 3);
                            tri[f][0] = qem.vertex(H[f][2]);
                            tri[f][1] = qem.vertex(H[next_f][1]);
                            tri[f][2] = qem.vertex(H[next_f][2]);
                        }

                        FOR(f, 3) {
                            have_bad_angle = have_bad_angle
                                || std::abs(cos_corner(X(m)[tri[f][0]], X(m)[tri[f][1]], X(m)[tri[f][2]])) > cos(M_PI/8.);
                        }
                        if (have_bad_angle) { failt[punch_v] = geo_max(failt[punch_v], 50.); continue; }
                        
                        FOR(f, 3) {
                            n[f] = normalize(cross(X(m)[tri[f][1]] - X(m)[tri[f][0]], X(m)[tri[f][2]] - X(m)[tri[f][0]]));
                            FOR(j, 3) mat(f, j) = n[f][j];
                        }
                        mat3 inv = mat.inverse();
                        vec3 b;
                        FOR(i, 3) b[i] = dot(n[i], X(m)[tri[i][0]]);
                        mult(inv, b.data(), new_vertex_position.data());
                        X(m)[nv_punch_v] = new_vertex_position;
                    }

                    if (new_hex_geometry_is_crappy()) {
                        X(m)[punch_v] = old_vertex_position;
                        { failt[punch_v] = geo_max(failt[punch_v], 60.); continue; }
                    }



                    if (punch_will_produce_intersection()) {
                        X(m)[punch_v] = old_vertex_position;
                        continue; 
                    }

                    // split vertex if needed (non manifold)
                    if (punch_v == nv_punch_v && v2nb_facets[punch_v] != 3) {
                        index_t nvv = m->vertices.create_vertex();
                        v2nb_facets[punch_v] -= 3;
                        v2nb_facets.push_back(3);
                        X(m)[punch_v] = old_vertex_position;
                        X(m)[nvv] = new_vertex_position;
                        punch_v = nvv;
                        nv_punch_v = nvv;
                    }

                    FOR(f, 3) { v2nb_facets[qem.vertex(H[f][1])]--; v2nb_facets[qem.vertex(H[f][2])]++; }
                    create_new_hex();
                    topo_punch();
                    FOR(lf, 3) hb.update_bbox(qem.face(H[lf][0]), facet_bbox(qem.face(H[lf][0])));
                    found_vertex_to_punch = true;

                    unglue_duplicates();

                    nb_punchs++;
                    if (!(nb_punchs%100)) plop(nb_punchs);
                    if (-1 == --nbmaxpunch) { 
                        plop(nbmaxpunch); 
                        goto cleanup;// return true;
                    }// debug only... to be removed
                }
                finished = !found_vertex_to_punch;
            }

        cleanup:

            qem.debug_export_adjacence();
            vector<index_t> to_kill(m->facets.nb(), true);
            FOR(f, m->facets.nb()) {
                index_t opp = qem.face(qem.opp(4 * f));
                for (index_t h = 4 * f; h < 4 * (f + 1); h++)
                    if (qem.face(qem.opp(h)) != opp) {
                        to_kill[f] = false;
                        to_kill[opp] = false;
                    }
            }         
            m->facets.delete_elements(to_kill);
            plop(nb_punchs);
            return nb_punchs > 0;
        }


    };

    void bourrin_subdivide_hexes(Mesh* m) {
        if (m->cells.nb() == 0) return;
        index_t nb_cells = m->cells.nb();
        index_t off_v = m->vertices.create_vertices(64 * nb_cells);
        index_t off_c = m->cells.create_hexes(8 * nb_cells);
        FOR(c, nb_cells) FOR(nc, 8)FOR(nv, 8)
            m->cells.set_vertex(off_c + 8 * c + nc, nv, off_v + 64 * c + 8 * nc + nv);
        FOR(c, nb_cells) {
            vec3 pts[8];
            FOR(cv, 8) pts[cv] = X(m)[m->cells.vertex(c, cv)];

            FOR(nci, 2)FOR(ncj, 2)FOR(nck, 2) {				// new cell position in cell c

                FOR(nvi, 2)FOR(nvj, 2)FOR(nvk, 2) {			// new vertex position in new cell
                    index_t nc = nci + 2 * ncj + 4 * nck;
                    index_t nv = nvi + 2 * nvj + 4 * nvk;
                    vec3 coeff[2];
                    coeff[0] = vec3(double(nci + nvi) / 2., double(ncj + nvj) / 2., double(nck + nvk) / 2.);
                    coeff[1] = vec3(1, 1, 1) - coeff[0];
                    vec3 pos = vec3(0, 0, 0);
                    FOR(i, 2)FOR(j, 2)FOR(k, 2) {
                        pos += coeff[i][0] * coeff[j][1] * coeff[k][2] * pts[i + 2 * j + 4 * k];
                    }
                    X(m)[off_v + 64 * c + 8 * nc + nv] = pos;
                }

            }
        }


        vector<index_t> to_kill(m->cells.nb(), false);
        FOR(c, nb_cells) to_kill[c] = true;
        m->cells.delete_elements(to_kill);

        mesh_repair(*m, MESH_REPAIR_COLOCATE, 1e-15);

    }

    void bourrin_quadrangulate_facets(Mesh* m) {
        if (m->facets.nb() == 0) return;
        Attribute<bool> isquad(m->facets.attributes(), "isquad");
        index_t nb_facets = m->facets.nb();
        FOR(f, nb_facets) {
            index_t nbv = m->facets.nb_vertices(f);
            vec3 bary(0, 0, 0);
            FOR(fv, nbv) bary = bary + (1. / double(nbv))*X(m)[m->facets.vertex(f, fv)];
            index_t nb_v = m->facets.nb_corners(f);
            FOR(lc, nb_v) {
                vec3 v[3] = {
                    X(m)[m->facets.vertex(f, prev_mod(lc, nb_v))],
                    X(m)[m->facets.vertex(f, lc)],
                    X(m)[m->facets.vertex(f, next_mod(lc, nb_v))]
                };
                index_t off_v = m->vertices.create_vertices(4);
                X(m)[off_v] = bary;
                X(m)[off_v + 1] = 0.5*(v[0] + v[1]);
                X(m)[off_v + 2] = v[1];
                X(m)[off_v + 3] = 0.5*(v[1] + v[2]);
                isquad[m->facets.create_quad(off_v + 0, off_v + 1, off_v + 2, off_v + 3)] = (nb_v == 4);
            }
        }
        {
            vector<index_t> to_kill(m->facets.nb(), false);
            FOR(f, nb_facets) to_kill[f] = true;
            m->facets.delete_elements(to_kill);
        }

        // merge vertices
        {
            vector<index_t> to_kill(m->vertices.nb(), 0);
            vector<index_t> old2new(m->vertices.nb());
            Geom::colocate(m->vertices.point_ptr(0), 3, m->vertices.nb(), old2new, 1e-15);
            FOR(f, m->facets.nb()) FOR(fv, 4) m->facets.set_vertex(f, fv, old2new[m->facets.vertex(f, fv)]);
            FOR(v, m->vertices.nb()) if (old2new[v] != v) to_kill[v] = NOT_AN_ID;
            m->vertices.delete_elements(to_kill);
        }


    }




    static void remove_scabs(Mesh* m,Mesh* newhex) {
        GEO::Logger::out("HexDom")  << "try to remove_scabs" <<  std::endl;

        Attribute<int> ft(m->facets.attributes(), "ft");            // facet type 

        Attribute<bool> isquad(m->facets.attributes(), "isquad"); 
        vector<index_t> to_kill(m->facets.nb(), false);             // ;( faces may be hurt and even killed in this fonction

        vector<index_t> local_id(m->vertices.nb(), NOT_AN_ID);  // ids in the array of new vertices (projected onto the boundary)


        QuadExtraConnectivity qem;
        qem.init(m);
        FOR(h, 4 * m->facets.nb()) {
            if (to_kill[qem.face(h)]) continue;
            if (qem.opp(h) == NOT_AN_ID) continue;
            if (ft[qem.face(qem.opp(h)) ]!=0) continue;

            vector<index_t> contour;
            bool fail = false;
            {index_t cir = h;
            do {
                contour.push_back(cir);
                if (qem.opp(qem.next(cir, 3)) == NOT_AN_ID) { fail = true; break; }
                if (isquad[qem.face(qem.opp(qem.next(cir, 3)))]) { fail = true; break; }
                if (!isquad[qem.face(cir)]) { fail = true; break; }
                if (qem.opp(cir) == NOT_AN_ID) { fail = true; break; }
                cir = qem.next(qem.opp(cir), 2);
            } while (cir != h);
            }
            if (fail) continue;
            
            //{// debug output
            //    plop(contour.size());
            //    index_t off_e = newhex->edges.create_edges(contour.size());
            //    index_t off_v = newhex->vertices.create_vertices(contour.size());
            //    FOR(c, contour.size()) X(newhex)[off_v + c] = X(m)[qem.vertex(contour[c])];
            //    FOR(c, contour.size()) newhex->edges.set_vertex(off_e + c, 0, off_v + c);
            //    FOR(c, contour.size()) newhex->edges.set_vertex(off_e + c, 1, off_v + ((c + 1) % contour.size()));
            //}


            FOR(c, contour.size()) ft[qem.face(contour[c])] = 1;

            vector<index_t> in_facets;
            vector<index_t> out_facets;
            in_facets.push_back(qem.face(qem.opp(qem.next(h))));
            out_facets.push_back(qem.face(qem.opp(qem.next(h,3))));

            ft[in_facets[0]] = 2;
            ft[out_facets[0]] = 3;

            vector<index_t> stack;
            stack.push_back(in_facets[0]);
            stack.push_back(out_facets[0]);
            while (!stack.empty() && !fail) {
                index_t f = stack.back(); 
                stack.pop_back();
                for (index_t cir = 4 * f; cir < 4 * f + 4; cir++) {
                    // need to have an opposite
                    if (qem.opp(cir) == NOT_AN_ID) {
                        fail = true;
                        break;
                    }
                    // do not link to another isquad type OR with the contour
                    index_t oppf = qem.face(qem.opp(cir));
                    if (ft[oppf] != 1 && isquad[oppf] != isquad[f] ) {
                        fail = true;
                        break;
                    }

                  


                    if (ft[oppf] == 0) {
                        // linked by a reasonably flat edge---only for inside until outside have better mesh quality            
                        if (dot(facet_normal(m, f), facet_normal(m, oppf))<.8 && ft[f] == 2) {
                            fail = true;
                            break;
                        }
                        ft[oppf] = ft[f];
                        if (ft[f] == 2) in_facets.push_back(oppf);
                        if (ft[f] == 3) out_facets.push_back(oppf);
                        stack.push_back(oppf);
                    }
                }                
            }
            if (fail) continue;

            // check that in_facets and outfacets are topo disks 
            index_t in_facets_nb_neig = 0;
            FOR(i, in_facets.size()) {
                geo_assert(isquad[in_facets[i]]);
                for (index_t cir = 4 * in_facets[i]; cir < 4 * in_facets[i] + 4; cir++) {
                    if (ft[qem.face(qem.opp(cir))] != 2)
                        in_facets_nb_neig++;
                }
            }
            if(in_facets_nb_neig != contour.size()) continue;// assert is only for debug, to be replaced by if 

            index_t out_facets_nb_neig = 0;
            FOR(i, out_facets.size()) {
                geo_assert(!isquad[out_facets[i]]);
                for (index_t cir = 4 * out_facets[i]; cir < 4 * out_facets[i] + 4; cir++) {
                    if (ft[qem.face(qem.opp(cir))] != 3)
                        out_facets_nb_neig++;
                }
            }
            if (out_facets_nb_neig != contour.size()) continue;// assert is only for debug, to be replaced by if


            //generate new vertices
            index_t nbv = 0;
            vector<vec3> packed_v; // for each vertex of infacets, stores the 3-uplet (pos, normal, projected pos)
            FOR(i, in_facets.size()) {
                for (index_t cir = 4 * in_facets[i]; cir < 4 * in_facets[i] + 4; cir++) {
                    index_t v = qem.vertex(cir);
                    if (local_id[v] == NOT_AN_ID) {
                        local_id[v] = nbv;
                        nbv++;
                        packed_v.push_back(X(m)[v]);
                        packed_v.push_back(vec3(0, 0, 0));
                        packed_v.push_back(vec3(0, 0, 0));
                    }
                    // accumulate facet normals into vertex normal
                    packed_v[3*local_id[v]+1] = packed_v[3 * local_id[v] + 1]  + facet_normal(m, in_facets[i]);
                }
            }




            // place nex vertices... TODO
            //FOR(lv, nbv) packed_v[3 * lv + 2] = packed_v[3 * lv + 0] - .01 * normalize(packed_v[3 * lv + 1]);

            //double ave_edge_length = 0;
            //FOR(f, m->facets.nb()) FOR(v, 4) ave_edge_length += (X(m)[m->facets.vertex(f, v)] - X(m)[m->facets.vertex(f, (v + 1) % 4)]).length();
            //ave_edge_length /= 4.*double(m->facets.nb());

            double ave_decal_length = 0;
            int nb_intersections = 0;
            FOR(v, nbv) {
                vec3 O = packed_v[3 * v];
                vec3 O2 = packed_v[3 * v] + normalize(packed_v[3 * v + 1]);
                FOR(i, out_facets.size()) {
                    FOR(tr, 2) {
                        index_t cir = 4 * out_facets[i] + 2 * tr; // each quad is decomposed into 2 triangles 
                        vec3 P[3];
                        FOR(p, 3)  P[p] = X(m)[qem.vertex(qem.next(cir, p))];
                        double sign[3];
                        FOR(p, 3) sign[p] = tetra_volume_sign(P[p], P[(p + 1) % 3], O, O2);
                        if (!same_sign(sign[0], sign[1]) || !same_sign(sign[0], sign[2])) continue;
                        double c0 = tetra_volume(P[0], P[1], P[2], O);
                        double c1 = tetra_volume(P[0], P[1], P[2], O2);
                        plop("intersection found");
                        packed_v[3 * v + 2] = O+(c0 / (c0 - c1))*(O2 - O);
                        ave_decal_length += (packed_v[3 * v + 2] - O).length();
                        nb_intersections++;
                        break;
                    }

                    if (packed_v[3 * v + 2].length2() != 0) break;
                }
            }
            ave_decal_length /=double(nb_intersections );

            FOR(v, nbv) {
                if (packed_v[3 * v + 2].length2() == 0) {
                    plop("panic mode no intersection found");
                    packed_v[3 * v + 2] = packed_v[3 * v + 0] - ave_decal_length * normalize(packed_v[3 * v + 1]);
                }  
            }



            // vertices on border must match the contour
            FOR(c, contour.size()) packed_v[3 * local_id[qem.vertex(qem.next(contour[c]))] + 2] = X(m)[qem.vertex(contour[c])];


            // generate hexes
            index_t off_c = newhex->cells.create_hexes(in_facets.size());
            index_t off_v = newhex->vertices.create_vertices(2* nbv);
            FOR(lv, nbv) {
                X(newhex)[off_v + 2 * lv] = packed_v[3 * lv ];
                X(newhex)[off_v + 2 * lv + 1] = packed_v[3 * lv+2];
            }
            FOR(i, in_facets.size()) {
                index_t f = in_facets[i];
                newhex->cells.set_vertex(off_c + i, 0, off_v+ 2 * local_id[qem.vertex(4 * f)]);
                newhex->cells.set_vertex(off_c + i, 1, off_v + 2 * local_id[qem.vertex(4 * f + 1)]);
                newhex->cells.set_vertex(off_c + i, 2, off_v + 2 * local_id[qem.vertex(4 * f + 3)]);
                newhex->cells.set_vertex(off_c + i, 3, off_v + 2 * local_id[qem.vertex(4 * f + 2)]);
                newhex->cells.set_vertex(off_c + i, 4, off_v + 1 + 2 * local_id[qem.vertex(4 * f)]);
                newhex->cells.set_vertex(off_c + i, 5, off_v + 1 + 2 * local_id[qem.vertex(4 * f + 1)]);
                newhex->cells.set_vertex(off_c + i, 6, off_v + 1 + 2 * local_id[qem.vertex(4 * f + 3)]);
                newhex->cells.set_vertex(off_c + i, 7, off_v + 1 + 2 * local_id[qem.vertex(4 * f + 2)]);
            }

            FOR(c, contour.size())      to_kill[qem.face(contour[c])] = true;
            FOR(i, in_facets.size())    to_kill[in_facets[i]] = true;
            FOR(i, out_facets.size())   to_kill[out_facets[i]] = true;

            FOR(v, m->vertices.nb()) local_id[v] = NOT_AN_ID; //not clearly needed, but brainless
        }
        m->facets.delete_elements(to_kill);
    }


    static void evaluate_edges_valence(Mesh* m) {
        Attribute<int> edge_angle(m->facet_corners.attributes(), "edgeangle");
        FOR(f, m->facets.nb()) FOR(h, 4) edge_angle[m->facets.corner(f, h)] = rand() % 3;

    }

    static bool next_crunch(Mesh* m, Mesh* newhex,int nb_max_punch) {
        GEO::Logger::out("HexDom")  << "next_crunch" <<  std::endl;
        //newhex->clear();

        if (m->facets.nb() == 0) return false;
        //static index_t nb_splits = 0;

        static index_t iter = index_t(-1);

	if(iter == index_t(-1)) {
	    iter = 0;
	} else {
	    iter++;
	}

        if (iter==0) evaluate_edges_valence(m);


        plop(iter);

        bool did_something = false;
        while (nb_max_punch > 0) {
            
            plop(nb_max_punch);
            VertexPuncher punch(m, newhex);
            if (punch.apply(nb_max_punch)) {
                did_something = true;
                if (nb_max_punch <= 0) return true;
            } else {
                //plop(did_something);
                if (did_something) return true;
                break;
            }
        }
        GEO::Logger::out("HexDom")  << "try to cut" <<  std::endl;
        static int cutit = 0;
        CutSingularity cut(m);
        if (cut.apply()) {
            GEO::Logger::out("HexDom")  << "------------------------" <<  std::endl;
	    GEO::Logger::out("HexDom")  << "CUT success... doing a small laplacian smoothing to separate collocated vertices" <<  std::endl;
            plop(cutit);
            cutit++;
            return true;
        }
        GEO::Logger::out("HexDom")  << "------------------------" <<  std::endl;
	GEO::Logger::out("HexDom")  << "CUT FAILED" <<  std::endl;

        remove_scabs(m,newhex);
        return false;
    }

    void prepare_crunch(Mesh* m,bool subdivide, bool revert) {
        if (subdivide) bourrin_quadrangulate_facets(m);
        if (revert) FOR(f, m->facets.nb()) {
            index_t save_v = m->facets.vertex(f, 0);
            m->facets.set_vertex(f, 0, m->facets.vertex(f, 2));
            m->facets.set_vertex(f, 2, save_v);
        }
    }

    void punch_and_cut(Mesh* m, Mesh* newhex, int nb_iter, int nb_punch_per_iter, bool check_validity) {
        FOR(i, nb_iter) {
            m->edges.clear();
            GEO::Logger::out("HexDom")  << "iteration " << i << " / " << nb_iter << "  with " << nb_punch_per_iter << "  punch per iter" <<  std::endl;
            if (!next_crunch(m, newhex, nb_punch_per_iter))break;
            if (check_validity) {
                m->edges.clear();
            }
        }
    }

    //bool in_volume(Mesh* surface, vec3 request) {
    //    int accum = 0;
    //    vec2 R(request[0], request[1]);
    //    FOR(f, surface->facets.nb()) {
    //        FOR(fan, surface->facets.nb_vertices(f) - 2) {
    //            index_t v[3] = {
    //                surface->facets.vertex(f,0),
    //                surface->facets.vertex(f,1 + fan),
    //                surface->facets.vertex(f,2 + fan)
    //            };
    //            vec2 P[3];
    //            FOR(vid, 3) P[vid] = vec2(X(surface)[v[vid]][0], X(surface)[v[vid]][1]) - R;
    //            bool in_triangle = true;
    //            double tr_orient = det(P[1] - P[0], P[2] - P[0]);
    //            if (tr_orient > 0)tr_orient = 1; else tr_orient = -1;
    //            FOR(vid, 3) in_triangle = in_triangle && ((det(P[vid], P[(vid + 1) % 3]) > 0) == (tr_orient >0)); ;
    //            //plop(in_triangle);
    //            if (!in_triangle) continue;

    //            int sign = orient_3d_filter(X(surface)[v[0]].data(), X(surface)[v[1]].data(), X(surface)[v[2]].data(), request.data());
    //            if (sign == 0) return false; // I don't want to deal with degenerate case
    //            accum += sign;
    //        }
    //    }
    //    //plop(accum);
    //    return accum != 0;

    //}

    //struct GrowPt {
    //    GrowPt(vec3 p_pos, AxisAngleRot p_r) { pos = p_pos; r = p_r; }
    //    vec3 pos;
    //    AxisAngleRot r;
    //};

 //   void Baudoin_mesher(Mesh* m) {
	//geo_argused(m);

        //// init rot
        //Attribute<bool> isquad(m->facets.attributes(), "isquad");
        //Attribute<mat3> B(m->vertices.attributes(), "B");

        //{
        //    vector<vector<vec3> >  dir(m->vertices.nb());
        //    FOR(f, m->facets.nb()) FOR(e, 4)
        //        dir[m->facets.vertex(f, e)].push_back(normalize(X(m)[m->facets.vertex(f, e)] - X(m)[m->facets.vertex(f, (e + 1) % 4)]));

        //    FOR(f, m->facets.nb()) {
        //        if (isquad[f]) continue;
        //        FOR(e, 4) dir[m->facets.vertex(f, e)].clear();
        //    }
        //    FOR(v, m->vertices.nb())
        //        if (dir[v].empty()) B[v].load_identity();
        //        else {
        //            B[v] = Frame::representative_frame(dir[v]);
        //            //if (rot[v].v.length2() < 1e-10) rot[v].v = vec3(1e-10, 0, 0);
        //        }
        //}

        //// compute ave QUAD edge length
        //double ave = 0;
        //double nb = 0;
        //FOR(f, m->facets.nb()) {
        //    if (!isquad[f]) continue;
        //    FOR(e, 4) {
        //        ave += (X(m)[m->facets.vertex(f, e)] - X(m)[m->facets.vertex(f, (e + 1) % 4)]).length();
        //        nb += 1.;
        //    }

        //}
        //ave /= nb;



        //vector<vec3> nvvertices;

        //vector<GrowPt> g;
        //int cur = 0;
        //// add useless boundary vertices (to prevent intersections)
        //FOR(v, m->vertices.nb()) if (B[v].is_identity()) g.push_back(GrowPt(X(m)[v], AxisAngleRot()));
        //cur = g.size();
        //FOR(v, m->vertices.nb()) if (!B[v].is_identity()) g.push_back(GrowPt(X(m)[v], B[v]));


        //// original paper

        //int nb_max_pts = 300;
        //while (cur < g.size() && nb_max_pts-->0) {
        //    FOR(d, 3) FOR(s, 2) {
        //        vec3 cand = ave * g[cur].r.rot_axis(d);
        //        if (s > 0) cand *= -1;
        //        cand = cand + g[cur].pos;
        //        bool fail = false;
        //        FOR(i, g.size())
        //            if ((g[i].pos - cand).length2() < pow(.5*ave, 2.))
        //                fail = true;
        //        if (!fail && in_volume(m, cand)) {
        //            g.push_back(GrowPt(cand, g[cur].r));
        //            nvvertices.push_back(cand);
        //        }
        //    }
        //    cur++;
        //}



        //index_t off_v = m->vertices.create_vertices(nvvertices.size());
        //FOR(i, nvvertices.size())  X(m)[off_v + i] = nvvertices[i];

    //}

}
