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

#include <exploragram/hexdom/PGP_export.h>
#include <exploragram/hexdom/PGP.h>
#include <exploragram/hexdom/basic.h>
#include <exploragram/hexdom/extra_connectivity.h>
#include <geogram/numerics/matrix_util.h>
#include <geogram/basic/permutation.h>
#include <algorithm>
#include <geogram/mesh/triangle_intersection.h>

#include <geogram/mesh/mesh_tetrahedralize.h>
#include <geogram/delaunay/delaunay.h>

#include <geogram/points/nn_search.h>
#include <geogram/points/colocate.h>
#include <queue>

#include <exploragram/hexdom/mesh_inspector.h>
#include <exploragram/hexdom/intersect_tools.h>
#include <exploragram/hexdom/polygon.h>


namespace GEO {

    /**
    * INPUT: facets with xyz geometry and uv coordinates s.t. no edge is of size 0
             singular bool per triangle
    *         integer values of uv --- interpolated along edge 'e' ---  matches on both side of 'e'
    on presuppose que les U ont été pre-snappés sur la grille entiere
    * OUTPUT: facets with uv coordinates s.t.
    *               a new vertex is inserted (with interpolated xyz and uv) at every integer values of uv along edges
    *               nothing prevents some vertices around a facet to share the same uv's
     ATTENTION: it modifies the parameterization when snaps grid corners on edges!
    */

    void split_edges_by_iso_uvs(Mesh* m, Attribute<vec2>& uv, Attribute<index_t> &singular) {
        index_t nb_init_facets = m->facets.nb();

        typedef std::pair<index_t /*v_id*/, vec2 /*uv*/> NewCorner;
        vector<vector<NewCorner> > new_corners(m->facet_corners.nb());

        reach("create vertices and init arrays of NewCorner to insert on edges");
        {
            FacetsExtraConnectivity fec(m);
            FOR(h, m->facet_corners.nb()) {
                index_t opp = fec.opposite(h);

                geo_assert(NOT_AN_ID != opp); // on suppose que la surface est 2-varieté

                if (singular[fec.facet(h)])             continue; // the face is responsible for the edge if opp<h or if its opposite is singular
                if (opp > h && !singular[fec.facet(opp)]) continue;

                vector<double> coeff;
                vec2 lU[2] = { uv[h], uv[fec.next(h)] };
                FOR(coord, 2) { // find barycentric coordinates of both u and v integer isos
                    double v[2] = { lU[0][coord], lU[1][coord] };
                    for (double iso = floor(std::min(v[0], v[1])) + 1; iso < std::max(v[0], v[1]); iso += 1.) {
                        geo_assert(std::abs(v[0] - v[1]) > 1e-3);
                        double c = (iso - v[0]) / (v[1] - v[0]); // v[0] is far from v[1] (U was pre-snapped to integers with .05 tolerance)
                        if (c > 0 && c < 1) coeff.push_back(c);
                    }
                }
                std::sort(coeff.begin(), coeff.end(), std::less<double>()); // we need to sort in order to merge close values

                vector<vec3> pts;
                vector<vec2> lu;
                vector<vec2> luopp;

                FOR(i, coeff.size()) {
                    // a + c*(b-a) returns exactly a when a==b and no NaNs involved
                    // no need to worry about the cases when c==0 and c==1, because of the presnapping of the parameterization
                    vec3 pt = X(m)[fec.org(h)] + coeff[i] * (X(m)[fec.dest(h)] - X(m)[fec.org(h)]);
                    vec2 u = lU[0] + coeff[i] * (lU[1] - lU[0]);
                    vec2 uopp = uv[fec.next(opp)] + coeff[i] * (uv[opp] - uv[fec.next(opp)]);

                    if (i) { // if a quad corner is close to the edge, it generates two points, let us merge them
                        vec2 &u_p = lu.back();
                        vec2 &uopp_p = luopp.back();
                        if (std::abs(u_p[0] - u[0]) < 1e-4 && std::abs(u_p[1] - u[1]) < 1e-4) { // ATTENTION! New mesh wont match perfectly the original hex mesh
                            geo_assert(round(u_p[0]) == round(u[0]) && round(u_p[1]) == round(u[1]));

                            FOR(d, 2) { // in addition to skipping this vertex we must not also round the previous point parameterization
                                u_p[d] = round(u_p[d]);
                                if (!singular[fec.facet(opp)]) {
                                    uopp_p[d] = round(uopp_p[d]);
                                }
                            }
                            continue;
                        }
                    }

                    FOR(d, 2) { // we must guarantee that new vertices have (at least) one integer component, do not forget pre-snapped U (.05 tolerance), thus it is safe to snap here at 1e-10
                        if (std::abs(u[d] - round(u[d])) < 1e-10) {
                            u[d] = round(u[d]);
                        }
                        if (std::abs(uopp[d] - round(uopp[d])) < 1e-10) {
                            uopp[d] = round(uopp[d]);
                        }
                    }

                    pts.push_back(pt);
                    lu.push_back(u);
                    luopp.push_back(uopp);
                }

                // create vertices
                index_t off = m->vertices.create_vertices(pts.size());
                FOR(i, pts.size()) {
                    m->vertices.point(off + i) = pts[i];
                    new_corners[h].push_back(NewCorner(off + i, lu[i]));
                    new_corners[opp].push_back(NewCorner(off + i, luopp[i]));
                }
                std::reverse(new_corners[opp].begin(), new_corners[opp].end());
            }
        }

        // now we have all the corners to insert, we create new facets for all the surface, old facets are to delete

        reach("split edges");
        FOR(f, nb_init_facets) {
            vector<index_t> polyV;
            vector<vec2> poly_uv;
            FOR(fc, m->facets.nb_corners(f)) {
                polyV.push_back(m->facets.vertex(f, fc));
                index_t c = m->facets.corner(f, fc);
                poly_uv.push_back(uv[c]);
                FOR(i, new_corners[c].size()) {
                    polyV.push_back(new_corners[c][i].first);
                    poly_uv.push_back(new_corners[c][i].second);
                }
            }
            index_t nf = m->facets.create_polygon(polyV);
            singular[nf] = singular[f];
            FOR(fc, m->facets.nb_corners(nf))
                uv[m->facets.corner(nf, fc)] = poly_uv[fc];
        }

        reach("kill facets");
        vector<index_t> to_kill(nb_init_facets, true); // kill old (pre-split) facets
        to_kill.resize(m->facets.nb(), false);
        m->facets.delete_elements(to_kill);
    }


    /**
    * INPUT:        facets with uv coordinates, where many vertices have integer values in uv
    * OUTPUT:       facets with uv coordinates, inclusing edges that have one coordinate of uv that is constant and integer valued
    *
    * remark:  some polylines of these new edges are likely to be pre-images of edges of the regular grid...
    */


    void facets_split_ca_va(Mesh* m, Attribute<vec2>& uv, Attribute<index_t> &singular) {
        vector<index_t> to_kill(m->facets.nb(), 0);
        FOR(f, m->facets.nb()) {
            if (singular[f]) continue;
            index_t nbc = m->facets.nb_corners(f);
            vector<index_t> hid(nbc);
            FOR(fc, nbc) hid[fc] = m->facets.corner(f, fc);

            // STEP 1: find a couple (coordinate, iso value) to split the facet
            index_t coord = index_t(-1);
            double iso = 0;
            FOR(test_coord, 2) {
                double min_v = 1e20;
                double max_v = -1e20;
                FOR(fc, nbc) {
                    min_v = std::min(min_v, uv[hid[fc]][test_coord]);
                    max_v = std::max(max_v, uv[hid[fc]][test_coord]);
                }
                if (floor(min_v) + 1 < max_v) {  // floor(min_v)+1 is the first integer value strictly superior to min_v
                    coord = test_coord;
                    iso = floor(min_v) + 1;
                    break;
                }
            }
            if (coord == index_t(-1)) continue;

            // STEP 2: if STEP 1 succedeed, compute the extremities of the new edge
            index_t cut[2] = { NOT_AN_ID, NOT_AN_ID };

            FOR(fc, nbc) {
                if (uv[hid[fc]][coord] != iso) continue;
                if (cut[0] == NOT_AN_ID) {
                    cut[0] = fc;
                }
                else if (fc != next_mod(cut[0], nbc) && next_mod(fc, nbc) != cut[0]) {
                    cut[1] = fc;
                }
            }

            geo_assert(cut[1] != NOT_AN_ID); // in fact, when two adjacent triangles are extracted with different coordinates e.g. (u,v) and (v,w), there might be conflicts on the incident edge, leading to the lack of the opposite on the non-responsible triangle
            geo_assert(next_mod(cut[1], nbc) != cut[0] && cut[1] != next_mod(cut[0], nbc));

            to_kill[f] = true;

            // STEP 3: compute new vertices (iso-integer value of uv) on the new edge
            vector<vec3> nv_pts;
            vector<vec2> nv_uv;
            {
                vec3 lX[2];
                FOR(i, 2) lX[i] = m->vertices.point(m->facets.vertex(f, cut[i]));
                vec2 lU[2];
                FOR(i, 2) lU[i] = uv[hid[cut[i]]];
                vector<double> coeff;
                double v[2] = { lU[0][(coord + 1) % 2], lU[1][(coord + 1) % 2] }; // recall that coord is the cutting dimension

                for (double cur_iso = ceil(std::min(v[0], v[1])); cur_iso < std::max(v[0], v[1]); cur_iso += 1.0) {
                    geo_assert(std::abs(v[0] - v[1]) > 1e-8); // TODO create a proper error handling; if we crash here, it means we have a VERY distorted triangle in the mesh (if I am not mistaken, 10^5 anisotropy)
                    double c = (cur_iso - v[0]) / (v[1] - v[0]);    // v[0] is far from v[1] (U was pre-snapped to integers with .05 tolerance)
                    if (c > 0 && c < 1)
                        coeff.push_back(c);
                }

                std::sort(coeff.begin(), coeff.end(), std::less<double>());
                FOR(i, coeff.size()) {
                    vec3 x = lX[0] + coeff[i] * (lX[1] - lX[0]); // it guarantees x==lX[0] when lX[0]==lX[1]
                    vec2 u = lU[0] + coeff[i] * (lU[1] - lU[0]); // no need to worry about coeff[i]==0 and coeff[i]==1 because of the parameterization pre-snapping
                    nv_pts.push_back(x);
                    nv_uv.push_back(u);
                }
                // new vertices must have only integer values of uv --- remove possible numerical imprecision
                FOR(i, nv_pts.size()) {
                    FOR(d, 2) {
                        nv_uv[i][d] = round(nv_uv[i][d]);
                    }
                }
            }

            // STEP 4: create new vertices and new faces
            index_t off = m->vertices.create_vertices(nv_pts.size());
            FOR(i, nv_pts.size()) X(m)[off + i] = nv_pts[i];
            FOR(half, 2) {
                vector <index_t> lv;
                vector <vec2> luv;
                //                    vector <bool> iso;

                                    // add original vertices
                index_t cir = cut[half];
                do {
                    //                        iso.push_back(cir == cut[half] ? true : false);
                    lv.push_back(m->facets.vertex(f, cir));
                    luv.push_back(uv[hid[cir]]);
                    cir = next_mod(cir, nbc);
                } while (cir != cut[(half + 1) % 2]);
                //                    iso.push_back(false);
                lv.push_back(m->facets.vertex(f, cir));
                luv.push_back(uv[hid[cir]]);

                // add new vertices
                FOR(i, nv_pts.size()) {
                    index_t ind = i;
                    if (half == 0) ind = nv_pts.size() - 1 - i;
                    lv.push_back(off + ind);
                    luv.push_back(nv_uv[ind]);
                    //                        iso.push_back(true);
                }

                index_t fid = m->facets.create_polygon(lv);
                FOR(fc, m->facets.nb_corners(fid)) {
                    uv[m->facets.corner(fid, fc)] = luv[fc];
                    //                        isovalue[m->facets.corner(fid, fc)] = iso[fc];
                }

                to_kill.push_back(false);
            }
        }
        m->facets.delete_elements(to_kill);
    }

    static void triangulate_surface_preserve_attributes(Mesh* m) {
        Attribute<int>     chart(m->facets.attributes(), "chart");   // TODO put these in parameters of the function
        Attribute<index_t> singular(m->facets.attributes(), "singular");
        Attribute<bool> quadelement(m->facets.attributes(), "quadelement");

        plop("triangulating surface after embedding isoUV");
        vector<index_t> to_kill(m->facets.nb(), false);
        FOR(f, m->facets.nb()) {
            geo_assert(m->facets.nb_corners(f) >= 3);
            if (m->facets.nb_corners(f) == 3) continue;

            vector<vec3> pts;
            FOR(lc, m->facets.nb_corners(f)) {
                pts.push_back(X(m)[m->facets.vertex(f, lc)]);
            }

            vector<index_t> triangles;
            bool success = Poly3d(pts).try_triangulate_minweight(triangles);
            geo_assert(success); // Our input polygons are guaranteed to be flat, convex and without too much distortion. Are'nt you happy?
            geo_assert(0 == triangles.size() % 3);

            FOR(new_face, triangles.size() / 3) {
                index_t new_f = m->facets.create_triangle(
                    m->facets.vertex(f, triangles[3 * new_face + 0]),
                    m->facets.vertex(f, triangles[3 * new_face + 1]),
                    m->facets.vertex(f, triangles[3 * new_face + 2]));
                chart[new_f] = chart[f];
                quadelement[new_f] = quadelement[f];
                singular[new_f] = singular[f];
                to_kill.push_back(false);
            }
            to_kill[f] = true;
        }
        m->facets.delete_elements(to_kill);
    }


    /**
    * INPUT:        facets with uv coordinates
    * OUTPUT:       chart facet attribute s.t. the chart frontier is included in edges that are iso-integer value of 'uv'
    */


    inline index_t indir_root(index_t i, vector<index_t>& indir) {
        while (i != indir[i]) i = indir[i];
        return i;
    }

    // merges charts for adjacent facets under two conditions:
    // 1) both facets are not singular
    // 2) the shared edge is not integer iso in both facets
    void mark_charts(Mesh* m, Attribute<vec2>& uv, Attribute<int>& chart, Attribute<index_t> &singular) { // 2-manifold surface is supposed
        vector<bool> isovalue(m->facet_corners.nb(), false);
        Attribute<bool> quadelement(m->facets.attributes(), "quadelement");

        FacetsExtraConnectivity fec(m);
        vector<index_t> indir(m->facets.nb());
        FOR(f, m->facets.nb()) {
            indir[f] = f;
        }

        FOR(h, m->facet_corners.nb()) {
            index_t hopp = fec.opposite(h);
            index_t f = fec.facet(h);
            index_t fopp = fec.facet(hopp);

            if (singular[f] || hopp > h) continue;

            bool cut = false;
            index_t test[2] = { h, hopp };
            FOR(lh, 2) {
                FOR(coord, 2) {
                    cut = cut || (uv[test[lh]][coord] == uv[fec.next(test[lh])][coord] && uv[test[lh]][coord] == round(uv[test[lh]][coord]));
                }
            }
            if (cut) {
                isovalue[h] = isovalue[hopp] = true;
            }

            if (!cut && !singular[fopp]) {
                indir[indir_root(fopp, indir)] = indir_root(f, indir);
            }
        }

        FOR(f, m->facets.nb()) {
            chart[f] = int(indir_root(f, indir));
        }

        // TODO this code verifies only if chart boundaries are marked as isovalues.
        // normally it is also necessary to check if the "quad" is a 2d disk (one boundary + Euler characterisic (imagine a torus with a hole))

        vector<bool> seen(m->facets.nb(), false);
        FOR(fseed, m->facets.nb()) {
            if (fseed != index_t(chart[fseed])) continue;

            std::deque<index_t> Q;
            Q.push_back(fseed);
            seen[fseed] = true;
            bool iso_only_at_boundaries = true;
            std::vector<index_t> C;
            while (Q.size()) {
                index_t f = Q.front();
                C.push_back(f);
                Q.pop_front();
                FOR(fc, m->facets.nb_corners(f)) {
                    index_t h = m->facets.corner(f, fc);
                    index_t hopp = fec.opposite(h);
                    geo_assert(hopp != NOT_AN_ID);
                    index_t fopp = fec.facet(hopp);
                    if (chart[fopp] != chart[fseed]) {
                        if (isovalue[h]) {
                            continue;
                        }
                        else {
                            Q.clear();
                            iso_only_at_boundaries = false;
                            break;
                        }
                    }
                    if (seen[fopp]) continue;
                    Q.push_back(fopp);
                    seen[fopp] = true;
                }
            }
            if (iso_only_at_boundaries) {
                FOR(i, C.size()) {
                    quadelement[C[i]] = true;
                }
            }
        }
        triangulate_surface_preserve_attributes(m);
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool try_remove_vertices_in_quad_charts(Mesh* m, vector<BBox>& locked_regions, bool on_boundary);
    bool lock_self_intersecting_regions(Mesh* facets_with_quads_and_tri_only, vector<BBox>& regions_to_lock);



    // Attention, sub-functions of this function need to access to attributes "chart" and "singular"
    void simplify_quad_charts(Mesh* m) {
        std::string msg;
        if (!surface_is_manifold(m, msg)) plop(msg);

        Mesh m_bak;
        m_bak.copy(*m);

        vector<BBox> locked_regions;
        int phase = 0;
        for(;;) {
            bool modified = false;
            if (0 == phase) {
                plop("try_remove_vertices_inside_quads(m, locked_regions)");
                modified = try_remove_vertices_in_quad_charts(m, locked_regions, false);
            }
            else {
                plop("try_remove_vertices_on_chart_frontier(m, locked_regions)");
                modified = try_remove_vertices_in_quad_charts(m, locked_regions, true);
            }

            plop("check for intersections");

            bool fail = lock_self_intersecting_regions(m, locked_regions);
            if (fail) {
                plop("conflict detected");
                m->copy(m_bak);
            }
            else {
                m_bak.copy(*m);
                if (!modified && ++phase > 1) break;
            }
        }
    }


    // test for self-intersection
    // input : a surface with only triangles and quads allowed
    // output : push_back into regions_to_lock BBoxes of self-intersections
    bool lock_self_intersecting_regions(Mesh* facets_with_quads_and_tri_only, vector<BBox>& regions_to_lock) {
        Mesh* m = facets_with_quads_and_tri_only;

        vector<BBox> inboxes(m->facets.nb());
        FOR(f, m->facets.nb()) {
            index_t nbv = m->facets.nb_vertices(f);
            geo_assert(4 == nbv || 3 == nbv);
            FOR(fv, nbv) {
                inboxes[f].add(X(m)[m->facets.vertex(f, fv)]);
            }
        }
        HBoxes hb(inboxes);

        bool conflict_detected = false;
        FOR(f, m->facets.nb()) {
            vector<index_t> primitives;
            hb.intersect(inboxes[f], primitives);
            FOR(i, primitives.size()) {
                index_t opp_f = primitives[i];
                vector<vec3> P;
                FOR(fv, m->facets.nb_vertices(f)) {
                    P.push_back(X(m)[m->facets.vertex(f, fv)]);
                }
                vector<vec3> Q;
                FOR(fv, m->facets.nb_vertices(opp_f)) {
                    Q.push_back(X(m)[m->facets.vertex(opp_f, fv)]);
                }

                bool conflict = false;
                FOR(trP, 4) {
                    FOR(trQ, 4) {
                        if (trP > 0 && P.size() == 3) continue;
                        if (trQ > 0 && Q.size() == 3) continue;
                        vector<TriangleIsect> trash;
                        conflict = conflict || triangles_intersections(
                            P[quad_split[trP][0]], P[quad_split[trP][1]], P[quad_split[trP][2]],
                            Q[quad_split[trQ][0]], Q[quad_split[trQ][1]], Q[quad_split[trQ][2]],
                            trash
                        );
                    }
                }

                if (conflict) {
                    regions_to_lock.push_back(inboxes[f]);
                    regions_to_lock.push_back(inboxes[opp_f]);
                }
                conflict_detected = conflict_detected || conflict;
            }
        }
        return  conflict_detected;
    }



    struct TrFan { // the idea is to call this guy in cases when we have either 1 chart (inside a quad) or exactly 2 charts in the fan (on a boundary of a quad)
        struct TrFanElt {
            TrFanElt(Mesh* m, index_t v, index_t p_f, int p_chart, int q) : f(p_f), chart(p_chart), quadelement(q != 0) {
                geo_assert(m->facets.nb_vertices(f) == 3);
                FOR(c, 3) {
                    if (m->facets.vertex(f, c) != v) continue;
                    org = m->facets.vertex(f, (c + 1) % 3);
                    dest = m->facets.vertex(f, (c + 2) % 3);
                    break;
                }
            }

            index_t f;	    
            int chart;
            bool quadelement;
            index_t org;
            index_t dest;
        };

        TrFan(index_t v, Mesh *m, vector<vector<index_t> > &v2f, Attribute<int> &chart_, Attribute<bool> &quadelement_) : m_(m), v_(v) {
            if (0 == v2f[v].size()) return;
            FOR(i, v2f[v].size()) { // collect triangles around v
                index_t f = v2f[v][i];
                fan.push_back(TrFanElt(m, v, f, chart_[f], quadelement_[f]));
            }

            FOR(f, fan.size() - 1) { // sort triangles in circular order, first triangle is not moved
                for (index_t i = f + 2; i < fan.size(); i++) {
                    if (fan[f].dest == fan[i].org) {
                        std::swap(fan[f + 1], fan[i]);
                        break;
                    }
                }
            }

            FOR(f, fan.size()) {
                if (fan[f].dest != fan[(f + 1) % fan.size()].org) {
                    GEO::Logger::out("HexDom")  << "Fan around vertex " << v_ << " is not valid" <<  std::endl;
                    FOR(ff, fan.size()) GEO::Logger::out("HexDom")  << "fan[ff].org = " << fan[ff].org << "\tfan[ff].dest = " << fan[ff].dest <<  std::endl;
                }
                geo_assert(fan[f].dest == fan[(f + 1) % fan.size()].org);
            }

            index_t rotate = 0; // find a boundary between charts (any boundary will do); if no boundary rotate will be zero
            FOR(f, fan.size()) {
                if (fan[rotate].chart != fan[(rotate - 1 + fan.size()) % fan.size()].chart) break;
                rotate++;
            }
            geo_assert(rotate <= fan.size());
            std::rotate(fan.begin(), fan.begin() + int(rotate), fan.end());

            chart_offset.push_back(0);
            FOR(f, fan.size() - 1) {
                if (fan[f].chart == fan[f + 1].chart) continue;
                chart_offset.push_back(int(f + 1));
            }
        }

        index_t ncharts() {
            return index_t(chart_offset.size());
        }

        bool quadelement(index_t ichart) {
            geo_assert(ichart < chart_offset.size());
            geo_assert(index_t(chart_offset[ichart]) < fan.size());
            return fan[chart_offset[ichart]].quadelement;
        }

        index_t chart(index_t ichart) {
            geo_assert(ichart < chart_offset.size());
            geo_assert(index_t(chart_offset[ichart]) < fan.size());
            return index_t(fan[chart_offset[ichart]].chart);
        }

        index_t nverts(index_t ichart) {
            geo_assert(ichart < chart_offset.size());
            if (1 == chart_offset.size()) { // if one chart only, no need to add vertices
                geo_assert(0 == ichart);
                return fan.size();
            }
            int off1 = chart_offset[ichart];
            int off2 = ichart + 1 < chart_offset.size() ? chart_offset[ichart + 1] : int(fan.size());
            return index_t(off2 - off1 + 1);
        }

        index_t vert(index_t ichart, index_t ivert) {
            geo_assert(ichart < chart_offset.size());
            int off1 = chart_offset[ichart];
            int off2 = ichart + 1 < chart_offset.size() ? chart_offset[ichart + 1] : int(fan.size());
            if (int(ivert) < off2 - off1) return fan[index_t(chart_offset[ichart]) + ivert].org;
            geo_assert(off2 >= 1);
            return fan[off2 - 1].dest;
        }

        bool triangulate() {
            bool result = true;
            FOR(ichart, ncharts()) {
                vector<vec3> pts3d;
                FOR(ivert, nverts(ichart)) {
                    pts3d.push_back(X(m_)[vert(ichart, ivert)]);
                }
                pts3d.push_back(X(m_)[v_]);
                Basis3d b(Poly3d(pts3d).normal());

                vector<vec2> pts2d;
                FOR(ivert, nverts(ichart)) { // attention it skips v_, it is useful only for computing normal
                    pts2d.push_back(b.project_xy(pts3d[ivert]));
                }

                triangles.push_back(vector<index_t>());
                if (nverts(ichart) > 2 && !Poly2d(pts2d).try_triangulate_minweight(triangles.back())) {
                    error("was not able to ear cut");
                    triangles.back().clear();
                    result = false;
                }

                geo_assert(0 == triangles.back().size() % 3);
            }
            return result;
        }

        bool belongs_to_a_quad() {
            FOR(ichart, ncharts()) {
                if (quadelement(ichart)) return true;
            }
            return false;
        }

        TrFanElt& operator[](int i) {
            geo_assert(i >= 0 && i < int(fan.size()));
            return fan[i];
        }

        TrFanElt& operator[](index_t i) {
            geo_assert(i < fan.size());
            return fan[i];
        }
	
        index_t nb_fan_triangles() {
            return index_t(fan.size());
        }

        vector<vector<index_t> > triangles;
        vector<TrFanElt> fan;
        vector<int> chart_offset;

        Mesh *m_;	
        index_t v_;	
    };
    
    // compute light connectivity (vertex to facets)
    static vector<vector<index_t> > generate_v2f(Mesh *m) {
        vector<vector<index_t> > v2f(m->vertices.nb());
        FOR(f, m->facets.nb()) {
            geo_assert(m->facets.nb_vertices(f) == 3);
            FOR(lc, m->facets.nb_corners(f)) {
                v2f[m->facets.vertex(f, lc)].push_back(f);
            }
        }
        return v2f;
    }

    static void fan_asserts(Mesh *m, vector<vector<index_t> > &v2f, Attribute<int> &chart, Attribute<bool> &quadelement) {
        //        std::cerr << "sanity check: triangulated surface with 3 distinct vertid for each triangle...";
        FOR(f, m->facets.nb()) {
            geo_assert(3 == m->facets.nb_corners(f));
            FOR(c1, m->facets.nb_corners(f)) {
                FOR(c2, m->facets.nb_corners(f)) {
                    if (c1 == c2) continue;
                    geo_assert(m->facets.vertex(f, c1) != m->facets.vertex(f, c2));
                }
            }
        }
        //        std::cerr << "ok\n";

        //        std::cerr << "sanity check: no valence 2 verts...";
        bool ok = true;
        FOR(v, v2f.size()) {
            if (v2f[v].size() > 0 && v2f[v].size() < 3) {
                ok = false;
            }
        }
        geo_assert(ok);

        //        std::cerr << "sanity check: all fans are correct...";
        FOR(v, m->vertices.nb()) {
            TrFan fan = TrFan(v, m, v2f, chart, quadelement);
        }
        //        std::cerr << "ok\n";
    }

    bool try_remove_vertices_in_quad_charts(Mesh* m, vector<BBox>& locked_regions, bool on_boundary) {
        Attribute<int>        chart(m->facets.attributes(), "chart");      // TODO document this
        Attribute<bool> quadelement(m->facets.attributes(), "quadelement");

        vector<index_t> to_kill(m->facets.nb(), false);
        vector<vector<index_t> > v2f = generate_v2f(m);
        bool mesh_is_modified = false;

        fan_asserts(m, v2f, chart, quadelement);

        FOR(v, m->vertices.nb()) {
            TrFan fan = TrFan(v, m, v2f, chart, quadelement);
            if (index_t(on_boundary ? 2 : 1) != fan.ncharts() || !fan.belongs_to_a_quad()) continue;

            {
                bool intersect_locked_region = false;
                BBox b;
                b.add(X(m)[v]);
                FOR(i, fan.nb_fan_triangles()) {
                    index_t vcur = fan[i].org;
                    b.add(X(m)[vcur]);
                }
                FOR(i, locked_regions.size()) {
                    intersect_locked_region = intersect_locked_region || locked_regions[i].intersect(b);
                }
                if (intersect_locked_region) continue;
            }

            if (!fan.triangulate()) {
                error("was not able to triangulate");
                continue;
            }

            { // block of geometrical and topological sanity checks, allows to reject the triangulation
                bool non_manifold_edge = false;
                FOR(ichart, fan.ncharts()) {
                    vector<index_t> &tri = fan.triangles[ichart];
                    FOR(t, tri.size() / 3) {
                        FOR(iv, 3) {
                            index_t a = fan.vert(ichart, tri[3 * t + iv]);
                            index_t b = fan.vert(ichart, tri[3 * t + (iv + 1) % 3]);

                            // (a,b) is an edge in new triangulation, let us check whether it is on the boundary of the chart
                            bool fan_boundary = false;
                            FOR(i, fan.nb_fan_triangles()) {
                                fan_boundary = fan_boundary || (a == fan[i].org && b == fan[i].dest);
                            }
                            if (fan_boundary) continue;

                            // here we know that (a,b) is an interior edge
                            FOR(f, v2f[a].size()) {
                                FOR(c, m->facets.nb_corners(v2f[a][f])) {
                                    non_manifold_edge = non_manifold_edge || (b == m->facets.vertex(v2f[a][f], c));
                                }
                            }
                        }
                    }
                }
                if (non_manifold_edge) {
                    error("this triangulation would introduce a non manifold edge, rejecting");
                    continue;
                }
            }


            mesh_is_modified = true;

            FOR(ichart, fan.ncharts()) {
                vector<index_t> &tri = fan.triangles[ichart];
                FOR(t, tri.size() / 3) {
                    index_t new_f = m->facets.create_triangle(
                        fan.vert(ichart, tri[3 * t + 0]),
                        fan.vert(ichart, tri[3 * t + 1]),
                        fan.vert(ichart, tri[3 * t + 2]));
                    quadelement[new_f] = fan.quadelement(ichart);
                    chart[new_f] = int(fan.chart(ichart));
                    to_kill.push_back(false);

                    FOR(lc, 3) {
                        v2f[m->facets.vertex(new_f, lc)].push_back(new_f);
                    }
                }
            }

            // update light connection info
            FOR(i, v2f[fan.v_].size()) {
                to_kill[v2f[fan.v_][i]] = true;
            }
            v2f[fan.v_].clear();

            FOR(i, fan.nb_fan_triangles()) {
                index_t f = fan[i].f;
                index_t org = fan[i].org;
                index_t dest = fan[i].dest;

                index_t vopp[2] = { org, dest };

                FOR(j, 2) {
                    FOR(fi, v2f[vopp[j]].size()) {
                        if (v2f[vopp[j]][fi] != f) continue;
                        v2f[vopp[j]][fi] = v2f[vopp[j]].back();
                        v2f[vopp[j]].pop_back();
                        break;
                    }
                }
            }
            //fan_asserts(m, v2f, chart, quadelement);
        }

        fan_asserts(m, v2f, chart, quadelement);

        m->facets.delete_elements(to_kill);
        kill_isolated_vertices(m);
        return mesh_is_modified;
    }


    static void try_export_quadtri_from_charts(Mesh* m, vector<BBox>& locked_regions) {
        Attribute<int>        chart(m->facets.attributes(), "chart");      // TODO document this
        Attribute<bool> quadelement(m->facets.attributes(), "quadelement");
        Attribute<bool> isquad(m->facets.attributes(), "isquad"); // TODO virer ça

        vector<index_t> to_kill(m->facets.nb(), false);

        index_t max_chart_no = 0;
        FOR(f, m->facets.nb()) {
            max_chart_no = std::max(max_chart_no, index_t(chart[f]) + 1);
        }
        vector<int> nbtri_in_chart(max_chart_no, 0);
        FOR(f, m->facets.nb()) {
            nbtri_in_chart[chart[f]]++;
        }

        FacetsExtraConnectivity fec(m);
        index_t nbf = m->facets.nb();
        FOR(f, nbf) {
            geo_assert(3 == m->facets.nb_corners(f));

            if (nbtri_in_chart[chart[f]] != 2 || !quadelement[f]) continue;
            FOR(ih, 3) {
                index_t h = m->facets.corner(f, ih);
                index_t hopp = fec.opposite(h);
                geo_assert(NOT_AN_ID != hopp);
                index_t fopp = fec.facet(hopp);
                geo_assert(NOT_AN_ID != fopp);
                if (chart[fopp] != chart[f]) continue;
                if (f < fopp) break;

                vector<index_t> pts;
                pts.push_back(fec.dest(h));
                pts.push_back(fec.dest(fec.next(h)));
                pts.push_back(fec.org(h));
                pts.push_back(fec.dest(fec.next(hopp)));

                bool intersect_locked_region = false;
                BBox b;
                FOR(v, 4) {
                    b.add(X(m)[pts[v]]);
                }
                FOR(i, locked_regions.size()) {
                    intersect_locked_region = intersect_locked_region || locked_regions[i].intersect(b);
                }

                if (!intersect_locked_region) {
                    index_t fid = m->facets.create_polygon(pts);
                    isquad[fid] = true;
                    to_kill.push_back(false);
                    to_kill[f] = to_kill[fopp] = true;
                }
            }
        }
        m->facets.delete_elements(to_kill);
    }

    void export_quadtri_from_charts(Mesh* m) {
        Mesh m_bak;
        m_bak.copy(*m);

        plop("export_quadtri_from_charts(Mesh* m, Mesh* quadmesh)");
        vector<BBox> locked_regions;
        for(;;) {
            try_export_quadtri_from_charts(m, locked_regions);
            bool fail = lock_self_intersecting_regions(m, locked_regions);
            if (fail) {
                plop("conflict dectected");
                m->copy(m_bak);
            }
            else {
                break;
            }
        }
    }

    struct IntersectQuadTriangle {
        IntersectQuadTriangle(Mesh* p_hex, Mesh* p_border, index_t p_c, index_t p_cf) {
            border = p_border;
            hex = p_hex;
            c = p_c;
            cf = p_cf;
            intersected = false;
            //std::cerr << "new intersection with an hex quad\n";
        }
        void operator()(index_t trid) {
            vec3 T[3];
            FOR(p, 3) T[p] = X(border)[border->facets.vertex(trid, p)];
            vec3 Q[4];
            FOR(p, 4) Q[p] = X(hex)[hex->cells.facet_vertex(c, cf, p)];

            //vec3 Qbary(0, 0, 0);
            //FOR (p, 4) Qbary = Qbary + .25* Q[p];
            //FOR (p, 4) Q[p] = .999999 *Q[p] + .000001*Qbary;

            vector<TriangleIsect> result;
            intersected = intersected || triangles_intersections(T[0], T[1], T[2], Q[0], Q[1], Q[2], result);
            intersected = intersected || triangles_intersections(T[0], T[1], T[2], Q[0], Q[2], Q[3], result);
            intersected = intersected || triangles_intersections(T[0], T[1], T[2], Q[0], Q[1], Q[3], result);
            intersected = intersected || triangles_intersections(T[0], T[1], T[2], Q[1], Q[2], Q[3], result);
            //std::cerr << "avt tps " << intersected<<"\n";
            FOR(it, 3)FOR(iq, 4) intersected = intersected || (T[it] - Q[iq]).length() == 0;
            //std::cerr << "after tps " << intersected << "\n";

        }
        bool intersected;
        index_t c;
        index_t cf;
        Mesh* hex;
        Mesh* border;
    };


    static void snap_vertices_to_ref_vertices(Mesh* m, Mesh* ref, double eps = 1e-10) {
        if (ref->vertices.nb() == 0) return;
        NearestNeighborSearch_var NN = NearestNeighborSearch::create(3, "default");
        NN->set_points(ref->vertices.nb(), (double*)ref->vertices.point_ptr(0));
        Attribute<bool> snapped(m->vertices.attributes(), "snapped");

        FOR(v, m->vertices.nb()) {
            index_t v_ref = NN->get_nearest_neighbor((double*)&(X(m)[v]));
            snapped[v] = ((X(m)[v] - X(ref)[v_ref]).length() < eps);
            if ((X(m)[v] - X(ref)[v_ref]).length() < eps)
                X(m)[v] = X(ref)[v_ref];
        }
    }






    void kill_intersecting_hexes(Mesh* hex) {
        int nb_intersecting_hex = 0;
        vector<index_t> to_kill(hex->cells.nb(), false);
        vector<BBox> inboxes(6 * hex->cells.nb());
        FOR(f, inboxes.size())FOR(cfv, 4)
            inboxes[f].add(X(hex)[hex->cells.facet_vertex(f / 6, f % 6, cfv)]);
        HBoxes hb(inboxes);

        FOR(c, hex->cells.nb()) FOR(cf, 6) {
            BBox b;
            vector<vec3> Q(4);
            FOR(cfv, 4) {
                Q[cfv] = X(hex)[hex->cells.facet_vertex(c, cf, cfv)];
                b.add(Q[cfv]);
            }
            vector<index_t> primitives;
            hb.intersect(b, primitives);

            FOR(i, primitives.size()) {
                index_t other_c = primitives[i] / 6;
                if (other_c == c) continue;
                if (to_kill[other_c]) continue;
                index_t other_cf = primitives[i] % 6;
                vector<vec3> P(4);
                FOR(other_cfv, 4) P[other_cfv] = X(hex)[hex->cells.facet_vertex(other_c, other_cf, other_cfv)];

                // check if opposite
                bool is_opp = false;
                FOR(v, 4) {
                    index_t it = 0;
                    while (it < 4 && (Q[it] - P[(v + it) % 4]).length2() == 0) it++;
                    is_opp = is_opp || (it == 4);
                }

                if (is_opp) continue;

                // check intersection
                vector<TriangleIsect> result;
                FOR(c0, 2)FOR(c1, 2) {
                    if (triangles_intersections(
                        P[quad_rand_split[c0][0]], P[quad_rand_split[c0][1]], P[quad_rand_split[c0][2]],
                        Q[quad_rand_split[c1][0]], Q[quad_rand_split[c1][1]], Q[quad_rand_split[c1][2]],
                        result
                    )) {
                        plop("auto intersection found");
                        to_kill[c] = true;
                        nb_intersecting_hex++;
                    }
                }
            }
        }
        plop(nb_intersecting_hex);
        hex->cells.delete_elements(to_kill);
    }


    void hex_set_2_hex_mesh(Mesh* hex, Mesh* quadtri) {
        if (hex->cells.nb() == 0) return;

        // merge vertices
        double eps = (1e-3)*get_cell_average_edge_size(hex);
        {
            vector<index_t> to_kill(hex->vertices.nb(), 0);
            vector<index_t> old2new(hex->vertices.nb());
            Geom::colocate(hex->vertices.point_ptr(0), 3, hex->vertices.nb(), old2new, eps);
            FOR(c, hex->cells.nb()) FOR(cv, 8) hex->cells.set_vertex(c, cv, old2new[hex->cells.vertex(c, cv)]);
            FOR(v, hex->vertices.nb()) if (old2new[v] != v) to_kill[v] = NOT_AN_ID;
            hex->vertices.delete_elements(to_kill);
        }
        snap_vertices_to_ref_vertices(hex, quadtri, eps);

        // remove duplicated hex
        {
            int nb_duplicated_hex = 0;
            vector<vec3> sumVpos(hex->cells.nb(), vec3(0, 0, 0));
            vector<index_t> to_kill(hex->cells.nb(), 0);
            FOR(c, hex->cells.nb()) {
                FOR(lv, 8) sumVpos[c] = sumVpos[c] + hex->vertices.point(hex->cells.vertex(c, lv));
            }
            vector<index_t> bary_old2new(hex->cells.nb());
            Geom::colocate((double*)(sumVpos.data()), 3, hex->cells.nb(), bary_old2new, eps);
            FOR(c, hex->cells.nb()) {
                if (bary_old2new[c] != c) {
                    to_kill[c] = NOT_AN_ID;
                    nb_duplicated_hex++;
                }
            }
            plop(nb_duplicated_hex);
            hex->cells.delete_elements(to_kill);
            plop(hex->cells.nb());

        }


        // remove bad shaped hex
        {
            int nb_bad_shaped_hex = 0;
            vector<index_t> to_kill(hex->cells.nb(), false);
            FOR(c, hex->cells.nb()) {
                FOR(cf, 6)FOR(cfv, 4) {
                    vec3 C = X(hex)[hex->cells.facet_vertex(c, cf, prev_mod(cfv, 4))];
                    vec3 A = X(hex)[hex->cells.facet_vertex(c, cf, cfv)];
                    vec3 B = X(hex)[hex->cells.facet_vertex(c, cf, next_mod(cfv, 4))];
                    if (std::abs(dot(normalize(B - A), normalize(C - A))) > .6) {
                        to_kill[c] = true;
                        nb_bad_shaped_hex++;
                    }
                }
            }
            plop(nb_bad_shaped_hex);
            hex->cells.delete_elements(to_kill);
            plop(hex->cells.nb());

        }

        // remove hexes that are linked by 3 vertices
        {
            index_t nb_hex_linked_by_3_vertices = 0;
            vector<index_t> to_kill(hex->cells.nb(), false);
            vector<vector<index_t> > v2hexface(hex->vertices.nb());
            FOR(c, hex->cells.nb()) FOR(cf, 6) FOR(cfv, 4) v2hexface[hex->cells.facet_vertex(c, cf, cfv)].push_back(6 * c + cf);
            FOR(v, hex->vertices.nb()) FOR(h0, v2hexface[v].size())FOR(h1, h0) {
                index_t hexf0 = v2hexface[v][h0];
                index_t hexf1 = v2hexface[v][h1];
                if (to_kill[hexf0 / 6] || to_kill[hexf1 / 6]) continue;
                index_t f[2][4];
                FOR(lc, 4) f[0][lc] = hex->cells.facet_vertex(hexf0 / 6, hexf0 % 6, lc);
                FOR(lc, 4) f[1][lc] = hex->cells.facet_vertex(hexf1 / 6, hexf1 % 6, lc);
                index_t i0 = index_t(-1), i1 = index_t(-1);
                FOR(lc, 4) if (f[0][lc] == v) i0 = lc;
                FOR(lc, 4) if (f[1][lc] == v) i1 = lc;
                geo_assert(i0 != index_t(-1) && i1 != index_t(-1));
                int nb_shared = 0;
                FOR(lc, 4) if (f[0][(i0 + lc) % 4] == f[1][(i1 + 4 - lc) % 4])nb_shared++;
                if (nb_shared == 3) {
                    to_kill[hexf0 / 6] = true;
                    nb_hex_linked_by_3_vertices++;
                }
            }
            plop(nb_hex_linked_by_3_vertices);
            hex->cells.delete_elements(to_kill);
            plop(hex->cells.nb());

        }
        kill_intersecting_hexes(hex);

        // remove hex that are incompatible with quadtri
        if (quadtri != NULL) {
            int nb_hex_incompatible_with_quadtri = 0;

            // remove all hex that touch border
            vector<index_t> to_kill(hex->cells.nb(), true);

            vector<BBox> inboxes(quadtri->facets.nb());
            FOR(f, quadtri->facets.nb())
                FOR(fv, quadtri->facets.nb_vertices(f))
                inboxes[f].add(X(quadtri)[quadtri->facets.vertex(f, fv)]);
            HBoxes hb(inboxes);


            // check for interesection
            FOR(c, hex->cells.nb()) {
                to_kill[c] = false;
                FOR(cf, 6) {
                    // Q contains the face of the hex PLUS two extra vertices to make a diamon shape
                    vector<vec3> Q;
                    Q.reserve(6);
                    FOR(cfv, 4) Q.push_back(X(hex)[hex->cells.facet_vertex(c, cf, cfv)]);
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

                    bool have_intersection = false;
                    bool have_tri_quad_intersect = false;
                    FOR(i, primitives.size()) {
                        bool is_quatri_face = false;
                        index_t f = primitives[i];
                        vector<vec3> P;
                        FOR(fv, quadtri->facets.nb_vertices(f))
                            P.push_back(X(quadtri)[quadtri->facets.vertex(f, fv)]);
                        geo_assert(P.size() < 5);

                        // check "is_quatri_face"
                        if (P.size() == 4) {
                            FOR(v, 4) {
                                index_t it = 0;
                                while (it < 4 && (Q[it] - P[(v + it) % 4]).length2() == 0) it++;
                                is_quatri_face = is_quatri_face || (it == 4);
                            }
                        }
                        else {
                            FOR(v, 3) {
                                index_t it = 0;
                                while (it < 3 && (Q[it] - P[(v + it) % 3]).length2() == 0) it++;
                                have_tri_quad_intersect = have_tri_quad_intersect || (it == 3);
                            }
                        }


                        // check for non degenerated intersections
                        if (!is_quatri_face) {
                            vector<TriangleIsect> result;
                            FOR(diam, 12) {
                                if (P.size() == 3)
                                    have_intersection = have_intersection || triangles_intersections(
                                        P[0], P[1], P[2],
                                        Q[diamon_split[diam][0]], Q[diamon_split[diam][1]], Q[diamon_split[diam][2]], result
                                    );
                                else {
                                    FOR(qu, 4)
                                        have_intersection = have_intersection || triangles_intersections(
                                            P[quad_split[qu][0]], P[quad_split[qu][1]], P[quad_split[qu][2]],
                                            Q[diamon_split[diam][0]], Q[diamon_split[diam][1]], Q[diamon_split[diam][2]], result
                                        );
                                }

                            }
                        }

                    }

                    if (have_intersection) to_kill[c] = true;
                    if (have_tri_quad_intersect) to_kill[c] = true;
                    if (have_intersection || have_tri_quad_intersect)
                        nb_hex_incompatible_with_quadtri++;
                }
            }
            plop(nb_hex_incompatible_with_quadtri);
            hex->cells.delete_elements(to_kill);
            plop(hex->cells.nb());

        }

        hex->cells.connect();
        hex->cells.compute_borders();

    }










    void merge_hex_boundary_and_quadtri(Mesh* hex, Mesh* quad) {
        if (hex->cells.nb() == 0) return;
        double eps = 1e-3*get_cell_average_edge_size(hex);

        Attribute<index_t> hexvid(quad->vertices.attributes(), "hexvid");
        FOR(v, quad->vertices.nb()) hexvid[v] = NOT_AN_ID;

        Attribute<bool> isquad(quad->facets.attributes(), "isquad");


        // add hex surface to quadt
        index_t off_v = quad->vertices.create_vertices(hex->vertices.nb());
        FOR(v, hex->vertices.nb()) X(quad)[off_v + v] = X(hex)[v];
        FOR(v, hex->vertices.nb()) hexvid[off_v + v] = v;

        index_t off_f = quad->facets.create_facets(hex->facets.nb(), 4);
        FOR(f, hex->facets.nb()) FOR(fc, 4)
            quad->facets.set_vertex(off_f + f, (3 - fc), off_v + hex->facets.vertex(f, fc));
        FOR(f, hex->facets.nb()) isquad[off_f + f] = true;


        // merge vertices from hex faces to quatri vertices

        {
            vector<index_t> to_kill(quad->vertices.nb(), 0);
            vector<index_t> old2new(quad->vertices.nb());
            FOR(v, quad->vertices.nb()) old2new[v] = v;

            NearestNeighborSearch_var NN = NearestNeighborSearch::create(3);
            NN->set_points(off_v, quad->vertices.point_ptr(0));

            for (index_t v = off_v; v < quad->vertices.nb(); v++) {
                index_t nearest = NN->get_nearest_neighbor(X(quad)[v].data());
                if ((X(quad)[v] - X(quad)[nearest]).length2() == 0) {
                    old2new[v] = nearest;
                    hexvid[nearest] = hexvid[v];
                }
            }

            FOR(f, quad->facets.nb()) FOR(fv, quad->facets.nb_vertices(f))
                quad->facets.set_vertex(f, fv, old2new[quad->facets.vertex(f, fv)]);
            FOR(v, quad->vertices.nb()) if (old2new[v] != v) to_kill[v] = NOT_AN_ID;
            quad->vertices.delete_elements(to_kill);
        }




        // remove useless quads
        vector<vec3> f_bary(quad->facets.nb(), vec3(0, 0, 0));

        FOR(f, quad->facets.nb()) FOR(fc, quad->facets.nb_vertices(f))
            f_bary[f] = f_bary[f] + (1. / quad->facets.nb_vertices(f))*quad->vertices.point(quad->facets.vertex(f, fc));

        vector<index_t> bary_old2new(quad->facets.nb());
        Geom::colocate((double*)(f_bary.data()), 3, quad->facets.nb(), bary_old2new, eps);

        vector<index_t> to_kill(quad->facets.nb(), 0);
        FOR(f, quad->facets.nb()) if (bary_old2new[f] != f) {
            if (!isquad[f]) continue;
            if (!isquad[bary_old2new[f]]) continue;

            bool have_same_vertices = true;
            FOR(fv1, quad->facets.nb_vertices(f)) {
                bool isin = false;
                FOR(fv2, quad->facets.nb_vertices(bary_old2new[f])) {
                    isin = isin || quad->facets.vertex(f, fv1) == quad->facets.vertex(bary_old2new[f], fv2);
                }
                have_same_vertices = have_same_vertices && isin;
            }
            if (!have_same_vertices) continue;
            to_kill[f] = true;
            to_kill[bary_old2new[f]] = true;
        }
        quad->facets.delete_elements(to_kill);
        create_non_manifold_facet_adjacence(quad);


        // debug output
        //Attribute<int> isquaddebug(quad->facets.attributes(), "id"); FOR (f, quad->facets.nb())isquaddebug[f] = isquad[f] ? 1 : 0;

    }



    void export_tb_faces_DEBUG(Mesh* hex, Mesh* tbfaces) {
        Attribute<bool> touch_border(hex->cell_facets.attributes(), "tb");
        tbfaces->vertices.create_vertices(hex->vertices.nb());
        FOR(v, tbfaces->vertices.nb()) X(tbfaces)[v] = X(hex)[v];
        FOR(c, hex->cells.nb()) FOR(cf, 6) {
            if (touch_border[hex->cells.facet(c, cf)]) {
                tbfaces->facets.create_quad(
                    hex->cells.facet_vertex(c, cf, 0),
                    hex->cells.facet_vertex(c, cf, 1),
                    hex->cells.facet_vertex(c, cf, 2),
                    hex->cells.facet_vertex(c, cf, 3)
                );
            }
        }
    }



    static void fill_quad_tri_surface(Mesh* m, bool with_pyramids) {
        if (m->vertices.nb() == 0)  return;
        m->edges.clear();

        if (with_pyramids) {

            vector<index_t> pyrindex;
            // Remark: spliting quads in 2 may be better on boundary... but it is not clear for constrained boundary... that would also require isquad attribute
            vector<index_t> to_kill(m->facets.nb(), 0);

            index_t init_nb_facets = m->facets.nb();
            FOR(f, init_nb_facets) {
                if (m->facets.nb_vertices(f) != 4) continue;

                index_t nvv = m->vertices.create_vertex();
                //vec3 n = facet_normal(m, f);
                double d = 0;
                FOR(e, 4) d += .25*(X(m)[m->facets.vertex(f, e)] - X(m)[m->facets.vertex(f, (e + 1) % 4)]).length();
                X(m)[nvv] = facet_bary(m, f);// +.2*d*n;
                to_kill[f] = 1;

                index_t off_f = m->facets.create_triangles(4);
                FOR(e, 4) {
                    to_kill.push_back(0);
                    m->facets.set_vertex(off_f + e, 0, m->facets.vertex(f, e));
                    m->facets.set_vertex(off_f + e, 1, m->facets.vertex(f, (e + 1) % 4));
                    m->facets.set_vertex(off_f + e, 2, nvv);
                }
                FOR(e, 4) pyrindex.push_back(m->facets.vertex(f, e));
                pyrindex.push_back(nvv);
            }
            m->facets.delete_elements(to_kill, false);

            m->facets.triangulate();
            create_non_manifold_facet_adjacence(m);
            try {
                mesh_tetrahedralize(*m, false, true, 1.);
                index_t off_c = m->cells.create_pyramids(pyrindex.size() / 5);
                FOR(p, pyrindex.size() / 5) FOR(lv, 5)
                    m->cells.set_vertex(off_c + p, lv, pyrindex[5 * p + lv]);
            } catch (const Delaunay::InvalidInput& error_report) {
                FOR(i, error_report.invalid_facets.size())
                    plop(error_report.invalid_facets[i]);
            }

        }
        else {// !with_pyramids
            m->facets.triangulate();
            create_non_manifold_facet_adjacence(m);
            try {
                mesh_tetrahedralize(*m, false, true, 1.);
            } catch (const Delaunay::InvalidInput& error_report) {
                FOR(i, error_report.invalid_facets.size())
                    plop(error_report.invalid_facets[i]);
            }

        }
    }

    void fill_cavity_with_tetgen(Mesh* input, Mesh* tri, bool propagate_hexvid, bool with_pyramid) {
        tri->copy(*input, false);
        fill_quad_tri_surface(tri, with_pyramid);


        if (propagate_hexvid) {
            Attribute<index_t> in_hexvid(input->vertices.attributes(), "hexvid");
            Attribute<index_t>  tri_hexvid(tri->vertices.attributes(), "hexvid");
            FOR(v, tri->vertices.nb()) tri_hexvid[v] = NOT_AN_ID;
            FOR(v, input->vertices.nb()) tri_hexvid[v] = in_hexvid[v];
            FOR(v, input->vertices.nb()) if ((X(input)[v] - X(tri)[v]).length2() > 0) {
                error("tetgen has jittered my vertices !! ");
                geo_assert_not_reached;
            }
        }

    }


    void add_hexes_to_tetmesh(Mesh* hex, Mesh* tet_mesh) {
        if (hex->cells.nb() == 0) return;
        Attribute<index_t> hexvid(tet_mesh->vertices.attributes(), "hexvid");
        vector<index_t> hex2tet_idmap(hex->vertices.nb(), NOT_AN_ID);

        index_t off_v = tet_mesh->vertices.create_vertices(hex->vertices.nb());
        FOR(v, hex->vertices.nb()) X(tet_mesh)[off_v + v] = X(hex)[v];

        FOR(v, hex->vertices.nb()) hexvid[off_v + v] = NOT_AN_ID;


        FOR(v, hex->vertices.nb()) hex2tet_idmap[v] = off_v + v;
        FOR(v, tet_mesh->vertices.nb())
            if (hexvid[v] != NOT_AN_ID) {
                geo_assert((X(hex)[hexvid[v]] - X(tet_mesh)[v]).length2() == 0); // we assume that tetgen does not change vertices ids on boundary
                hex2tet_idmap[hexvid[v]] = v;
            }
        index_t off_c = tet_mesh->cells.create_hexes(hex->cells.nb());
        FOR(c, hex->cells.nb())FOR(cv, 8) {
            index_t vid = hex2tet_idmap[hex->cells.vertex(c, cv)];
            tet_mesh->cells.set_vertex(off_c + c, cv, vid);
        }
    }


}





