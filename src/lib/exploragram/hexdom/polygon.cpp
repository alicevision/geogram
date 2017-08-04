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

#include <exploragram/hexdom/polygon.h>

namespace GEO {

    vec2 Poly2d::barycenter() {
	vec2 bary(0, 0);
	FOR(fv, index_t(pts.size())) {
	    bary = bary + (1. / double(pts.size()))*pts[fv];
	}
	return bary;
    }

    void Poly2d::dump_contour() {
	index_t nbv = pts.size();
	Mesh export_mesh;
	export_mesh.vertices.create_vertices(nbv);
	FOR(i, nbv) X(&export_mesh)[i] = vec3(pts[i][0], pts[i][1], 0);
	//            FOR(i, nbv) plop(X(&export_mesh)[i]);
	vector<index_t> num;
	FOR(i, nbv) num.push_back(i);
	export_mesh.facets.create_polygon(num);
	char filename[1024];
	sprintf(filename, "nimp%i.obj", rand());
	mesh_save(export_mesh, filename);
    }

    // returns 1024. if concave angle is encountered or if proposed triangle contains one of pts
    // otherwise returns max angle of the proposed triangle
    double Poly2d::cost(index_t i, index_t j, index_t k) {
	vec2 C[3] = { pts[i], pts[j], pts[k] };
	double m = 0;
	FOR(v, 3) {
	    // note that angle is not the angle inside the triangle, but its complement
	    // angle variable has the "direction" information, thus it is negative for concave angles (right turn) and positive for convex angles (left turn)
	    double angle = atan2(
  	        det(C[(v + 1) % 3] - C[(v + 0) % 3], C[(v + 2) % 3] - C[(v + 1) % 3]),
		dot(C[(v + 1) % 3] - C[(v + 0) % 3], C[(v + 2) % 3] - C[(v + 1) % 3])
	    );
	    if (angle <= 0) return 1024.;
	    m = std::max(m, M_PI - angle);
	}

	FOR(other, pts.size()) {
	    if (other == i || other == j || other == k) continue;
	    vec2 P = pts[other];
	    bool inside = true;
	    FOR(l, 3) {
		inside = inside && (det(normalize(C[(l + 1) % 3] - C[l]), normalize(P - C[l])) > 0);
	    }
	    if (inside) return 1024.;
	}
	return m;
    }

    bool Poly2d::try_triangulate_minweight(vector<index_t>& triangles) {
	triangles.clear();
	index_t n = pts.size();
	geo_assert(n >= 3);

	//if (n == 3) {
	//    FOR(v, 3) {
	//        triangles.push_back(v);
	//    }
	//    return true;
	//}

	// we store in this table results of subproblems
	// table[i*n + j] stores the triangulation cost for points from i to j
	// the entry table[0*n + n-1] has the final result.
	std::vector<double> table(n*n, 0.);

	// this table stores triangle indices: for each subproblem (i,j) we have table[i*n + j]==k, i.e. the triangle is (i,k,j)
	std::vector<index_t> tri(n*n, index_t(-1));

	// note that the table is filled in diagonals; elements below main diagonal are not used at all
	for (index_t pbsize = 2; pbsize<n; pbsize++) {
	    for (index_t i = 0, j = pbsize; j<n; i++, j++) {
		// recall that we are testing triangle (i,k,j) which splits the problem (i,j) into
		// two smaller subproblems (i,k) and (k,j)
		double minv = 1024.;
		index_t mink = index_t(-1);
		for (index_t k = i + 1; k<j; k++) {
		    double val = table[i*n + k] + table[k*n + j] + cost(i, k, j);
		    if (minv <= val) continue;
		    minv = val;
		    mink = k;
		}
		table[i*n + j] = minv;
		tri[i*n + j] = mink;
	    }
	}

	if (table[n - 1] >= 1024.) {
	    dump_contour();
	    return false;
	}

	vector<index_t> Q(1, n - 1);
	FOR(t, Q.size()) {
	    index_t idx = Q[t];

	    index_t i = idx / n;
	    index_t k = tri[idx];
	    index_t j = idx % n;

	    //geo_assert(i >= 0 && k >= 0 && j >= 0);
	    triangles.push_back(i);
	    triangles.push_back(k);
	    triangles.push_back(j);

	    if (k + 2 <= j) Q.push_back(k*n + j);
	    if (i + 2 <= k) Q.push_back(i*n + k);
	}
	return true;
    }


    // find parity of original points
    index_t Poly2d::parity_of_original_points() {
	index_t nbv = pts.size();
	index_t dec = 0;
	double dec_score[2] = { 0, 0 };
	FOR(q, nbv / 2) {
	    FOR(d, 2)
		dec_score[d] = geo_max(dec_score[d], std::abs(
							      det(pts[(q * 2 + 0 + d) % nbv] - pts[(q * 2 + 1 + d) % nbv],
								  pts[(q * 2 + 2 + d) % nbv] - pts[(q * 2 + 1 + d) % nbv])));
	}
	if (dec_score[1] < dec_score[0])dec = 1;
	return dec;
    }


    bool Poly2d::middle_point_quadrangulate(vector<index_t>& quads) {
	index_t nbv = pts.size();
	vec2 G = barycenter();
	index_t dec = parity_of_original_points();
	FOR(q, nbv / 2) {
	    quads.push_back(nbv);
	    FOR(v, 3) quads.push_back((q * 2 + 1 - dec + v) % nbv);
	}
	pts.push_back(G);
	return true;
    }


    bool Poly2d::quads_are_valid(vector<index_t>& quads) {
	FOR(q, quads.size() / 4) {
	    FOR(e, 4) {
		vec2 v0 = normalize(pts[quads[4 * q + next_mod(e, 4)]] - pts[quads[4 * q + e]]);
		vec2 v1 = normalize(pts[quads[4 * q + prev_mod(e, 4)]] - pts[quads[4 * q + e]]);
		if (det(v0, v1)<sin(M_PI / 8.)) return false;
	    }
	}
	return true;
    }


    bool Poly2d::try_quadrangulate(vector<index_t>& quads) {

	bool verbose = false;

	index_t nbv = pts.size();
	if (verbose) plop(nbv);
	if (nbv < 4) return false;
	if (nbv == 4) {
	    FOR(v, 4) quads.push_back(v);
	    if (!quads_are_valid(quads)) { GEO::Logger::out("HexDom")  << "FAIL" <<  std::endl; return false; }
	    return true;
	}
	if (nbv % 2 != 0) {
	    GEO::Logger::out("HexDom")  << "There is no way to quadrangulate a surface with an odd number of boundary edges" <<  std::endl;
	    return false;
	}


	// precompute a few things
	vector<double> angle(nbv);
	vector<double> length(nbv);
	double ave_length = 0;
	FOR(i, nbv) {
	    vec2 P[3];
	    FOR(p, 3) P[p] = aupp(i + p - 1, pts);
	    angle[i] = (180. / M_PI)*atan2(det(P[1] - P[0], P[2] - P[1]), dot(P[1] - P[0], P[2] - P[1]));
	    if (verbose)GEO::Logger::out("HexDom")  << "i= " << i << "angle = " << angle[i] <<  std::endl;
	    length[i] = (P[1] - P[0]).length() / double(nbv);
	    ave_length += length[i];
	}


	// define outputs of the search
	index_t start = index_t(-1);
	index_t end = index_t(-1);
	index_t nb_nv_pts = index_t(-1);
	double best_score = 0;



	index_t dec = parity_of_original_points();

	FOR(test_start, nbv) {
	    index_t test_end;
	    index_t test_nb_nv_pts;
	    double test_score;


	    FOR(d, nbv - 5) {
		test_end = test_start + d + 3;

		vec2 A[3];
		FOR(i, 3) A[i] = aupp(int(test_start + i) - 1, pts);
		vec2 B[3];
		FOR(i, 3) B[i] = aupp(int(test_end + i) - 1, pts);
		vec2 nA1A2 = normalize(A[2] - A[1]);
		vec2 nA1A0 = normalize(A[0] - A[1]);
		vec2 nB1B2 = normalize(B[2] - B[1]);
		vec2 nB1B0 = normalize(B[0] - B[1]);
		vec2 nAB = normalize(B[1] - A[1]);
		vec2 nBA = -nAB;

		double worst_det = 1;
		worst_det = geo_min(worst_det, det(nA1A2, nAB));
		worst_det = geo_min(worst_det, det(nAB, nA1A0));
		worst_det = geo_min(worst_det, det(nB1B2, nBA));
		worst_det = geo_min(worst_det, det(nBA, nB1B0));

		test_score = worst_det;


		double AB_relative_length = floor((B[1] - A[1]).length() / ave_length);
		test_nb_nv_pts = index_t(geo_max(0, int(AB_relative_length) - 1));

		if (test_nb_nv_pts % 2 != int(d % 2)) {
		    if (test_nb_nv_pts == 0) test_nb_nv_pts = 1; else test_nb_nv_pts--;
		}

		if (angle[test_start] < 1) test_score += 1;
		if (angle[test_end] < 1) test_score += 1;
		if (angle[test_start] < -45) test_score += 2;
		if (angle[test_end] < -45) test_score += 2;
		if ((test_start % 2) == 1 - dec) test_score -= 10;
		if ((test_end % 2) == 1 - dec) test_score -= 10;
		test_nb_nv_pts = 1;

		if (best_score < test_score) {
		    bool can_cut = true;
		    FOR(dd, nbv) {
			index_t ind = test_start + dd;
			if (ind > test_start && ind < test_end)
			    can_cut = can_cut && det(nAB, aupp(ind, pts) - A[1])<0;
			if (ind > test_end)
			    can_cut = can_cut && det(nAB, aupp(ind, pts) - A[1])>0;
		    }
		    if (verbose)
			std::cerr << "can_cut = " << can_cut << "  test_score = " << test_score
				  << "   test_start = " << test_start << "   test_end = " << test_end
				  << "   test_nb_nv_pts = " << test_nb_nv_pts << std::endl;
		    if (can_cut) {
			start = test_start;
			end = test_end;
			nb_nv_pts = test_nb_nv_pts;
			best_score = test_score;
		    }
		}
	    }


	}


	if (nbv>8)
	    if (nb_nv_pts != index_t(-1)) {
		if (verbose)GEO::Logger::out("HexDom")  << "remove quad strip from " << start << " with score = " << best_score << " with nbpts" << nb_nv_pts <<  std::endl;

		vector<index_t> global_vid[2]; // gives indices in "pts" from indices in "poly[i]"

		// fill both half with existing points
		//int end = start + nb_nv_pts + 3;
		FOR(d, end - start + 1)                                 global_vid[0].push_back((start + d) % nbv);
		FOR(d, nbv - (end - start) + 1)         global_vid[1].push_back((end + d) % nbv);


		// add new vertices along the cut
		FOR(i, nb_nv_pts)                                       global_vid[0].push_back(nbv + i);
		FOR(i, nb_nv_pts)                                       global_vid[1].push_back(nbv + (nb_nv_pts - 1 - i));
		FOR(i, nb_nv_pts) {
		    double c = 1.0 - double(i + 1) / double(nb_nv_pts + 1);
		    pts.push_back((1. - c)*pts[start] + c*pts[end% nbv]);
		}

		// solve on two halves
		vector<vec2> poly[2];
		FOR(i, 2) FOR(fv, global_vid[i].size()) poly[i].push_back(pts[global_vid[i][fv]]);

		vector<index_t> poly_quad[2];
		FOR(i, 2) if (!Poly2d(poly[i]).try_quadrangulate(poly_quad[i])) return false;
		// add new pts to global
		FOR(i, 2) for (index_t d = global_vid[i].size(); d < poly[i].size(); d++) {
		    global_vid[i].push_back(pts.size());
		    pts.push_back(poly[i][d]);
		}
		FOR(i, 2) FOR(qu, poly_quad[i].size()) quads.push_back(global_vid[i][poly_quad[i][qu]]);

		if (!quads_are_valid(quads)) { GEO::Logger::out("HexDom")  << "FAIL" <<  std::endl; return false; }
		return true;
	    }




	if (verbose) GEO::Logger::out("HexDom")  << "middle_point_quadrangulate(quads)" <<  std::endl;
	middle_point_quadrangulate(quads);
	if (!quads_are_valid(quads)) { GEO::Logger::out("HexDom")  << "FAIL" <<  std::endl; return false; }

	return true;
    }
    
    /*****************************************************************************************************/

    vec3 Poly3d::barycenter() {
	vec3 bary(0, 0, 0);
	FOR(fv, pts.size()) {
	    bary = bary + (1. / double(pts.size()))*pts[fv];
	}
	return bary;
    }

    vec3 Poly3d::normal() {
	vec3 n(0, 0, 0);
	vec3 bary = barycenter();
	FOR(fv, pts.size()) {
	    n = n + cross(pts[fv] - bary, pts[next_mod(fv, pts.size())] - bary);
	}
	n = normalize(n);
	return n;
    }


    bool Poly3d::try_triangulate_minweight(vector<index_t>& triangles) {
	index_t nbv = pts.size();
	if (nbv == 3) {
	    FOR(v, 3) {
		triangles.push_back(v);
	    }
	    return true;
	}
	geo_assert(nbv > 3);

	vector<vec2> pts2d;
	Basis3d b(normal());
	FOR(fv, nbv) {
	    pts2d.push_back(b.project_xy(pts[fv]));
	}

	return Poly2d(pts2d).try_triangulate_minweight(triangles);
    }

    /**
     * WARNING: it may introduce new vertices in pts
     */
    bool Poly3d::try_quadrangulate(vector<index_t>& quads) {
	index_t nbv = pts.size();
	geo_assert(nbv>0);
	vec3 G = barycenter();
	Basis3d b(normal());

	vector<vec2> pts2d;
	FOR(fv, nbv) pts2d.push_back(b.project_xy(pts[fv] - G));
	Poly2d p2d(pts2d);
	if (!p2d.try_quadrangulate(quads)) return false;
	for (index_t i = pts.size(); i < p2d.pts.size(); i++)
	    pts.push_back(G + b.un_project_xy(p2d.pts[i]));
	return true;
    }
    
} 
