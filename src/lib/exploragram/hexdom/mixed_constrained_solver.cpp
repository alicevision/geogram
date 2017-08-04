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

#include <exploragram/hexdom/mixed_constrained_solver.h>

namespace {
    using namespace GEO;

    void show_eq(std::vector<Coeff> & eq){
	FOR(i, eq.size()) {
	    std::cerr << eq[i].a << ".X[" << eq[i].index << "]   "; 
	    if (i + 1 < eq.size()) std::cerr << "+ ";
	}
	std::cerr << "\n";
    }

    inline bool coeffCmp(const Coeff& A, const Coeff& B){
	return (A.index> B.index);
    }

    inline bool eqCmp(const std::vector<Coeff> &A, const std::vector<Coeff> &B){
	if (A.size() != B.size()) return A.size() > B.size();
	for (size_t d = 0; d < A.size(); d++){
	    if (A[d].index < B[d].index) return true;
	    if (A[d].index > B[d].index) return false;
	}
	return false;
    }
    
    inline  bool eqeq(const std::vector<Coeff> &A, const std::vector<Coeff> &B){
	if (A.size() != B.size()) return false;
	for (size_t d = 0; d < A.size(); d++){
	    if (A[d].index != B[d].index) return false;
	}
	return true;
    }

    inline double nint(double x) {
	return floor(x + .5);
    }

    /*
    inline double d2nint(double x) {
	return std::fabs(x - nint(x));
    }
    */
    
}

namespace GEO {

    void NaiveSparseMatrix::show_dense(const char* header) {
	GEO::Logger::out("HexDom")  << "-------------------------------------------------------" <<  std::endl;
	std::cerr << header;
	FOR(l, M.size()) {
	    FOR(c, nb_col) {
		double val = 0;
		FOR(i, M[l].size()) if (M[l][i].index == c) val = M[l][i].a;
		std::cerr << "\t" << val;
	    }
	    GEO::Logger::out("HexDom")  <<  std::endl;
	}
    }
	    
    void NaiveSparseMatrix::show_sparse(const char* header, index_t nb_min_terms) {
	GEO::Logger::out("HexDom")  << "-------------------------------------------------------" <<  std::endl;
	std::cerr << header;
	FOR(l, M.size())  if (M[l].size() >= nb_min_terms) {
	    std::cerr << "Line "<< l<<":  " ;
	    show_eq(M[l]);
	}
    }


    MatrixMixedConstrainedSolver::MatrixMixedConstrainedSolver(index_t p_nb_vars) {
	nb_vars = p_nb_vars;
	intvar.resize(nb_vars, false);
	V.resize(nb_vars);
	C.init(nb_vars, 0);
	state = SET_INTEGER_CONSTRAINTS;
    }

    index_t MatrixMixedConstrainedSolver::first_term_of_in_equals_out(std::vector<Coeff> &in, std::vector<Coeff> &out) {
	geo_assert(in.size() > 0);
	geo_assert(std::fabs(in[0].a) > .001);
	index_t res = in[0].index;
	double scale = -(1. / in[0].a);
	for (index_t i = 1; i < in.size(); i++) {
	    out.push_back(in[i]);
	    out.back().a *= scale;
	}
	return res;
    }

    void MatrixMixedConstrainedSolver::sort_and_compress(std::vector<Coeff> &eq) {
	if (eq.empty()) return;
	std::sort(eq.begin(), eq.end(), coeffCmp);
	size_t nn = 0;
	for (index_t i = 1; i < eq.size(); i++) {
	    if (eq[i].index == eq[i - 1].index) eq[nn].a += eq[i].a;
	    else {
		if (std::fabs(eq[nn].a) > .0001) nn++; // otherwise it "erases" the coefficient
		eq[nn] = eq[i];
	    }
	}
	if (std::fabs(eq[nn].a) > .0001) eq.resize(nn + 1);
	else eq.resize(nn);
    }

    void MatrixMixedConstrainedSolver::check_G_validity() {
	GEO::Logger::out("HexDom")  << "check_G_validity ... " <<  std::endl;
	FOR(cl, C.nb_lines()) if (C.nb_coeffs_in_line(cl) == 2) {
	    index_t i = C.ith_coeff_in_line(cl, 0).index;
	    index_t j = C.ith_coeff_in_line(cl, 1).index;
	    if (G.line(i)[0].index != G.line(j)[0].index) std::cerr << "Groups are not correct !\n";
	    else if (
		     (G.line(i)[0].a*G.line(j)[0].a)
		     *
		     (C.ith_coeff_in_line(cl, 0).a*C.ith_coeff_in_line(cl, 1).a)
		     > 0) {
		std::cerr << "Bad sign !\n";
		nullgrp.push_back(int(G.line(i)[0].index));
	    }
	}
	GEO::Logger::out("HexDom")  << "OK" <<  std::endl;
    }
    
    void MatrixMixedConstrainedSolver::end_linear_constraints() {
	geo_assert(state == SET_LINEAR_CONSTRAINTS);
	// construct grps
	vector<index_t> group(nb_vars);
	vector<bool> sign(nb_vars);
	IdMapSigned::make_identity(group, sign);
	index_t nb_groups;
	FOR(cl, C.nb_lines()) if (C.nb_coeffs_in_line(cl) == 2) {
	    IdMapSigned::merge(group, sign, 
			       C.ith_coeff_in_line(cl,0).index, 
			       C.ith_coeff_in_line(cl, 1).index,
			       (C.ith_coeff_in_line(cl, 0).a*C.ith_coeff_in_line(cl, 1).a) < 0);
	}
	
	IdMapSigned::tree_to_map(group, sign);
	IdMapSigned::compress_id2grp(group, nb_groups);
	GEO::Logger::out("HexDom")  << "Create " << nb_groups << " from " << group.size() << " variables" <<  std::endl;
	G.init(nb_groups, nb_vars);
	FOR(i,nb_vars) G.add_coeff_to_line(group[i], i, (sign[i]) ? 1 : (-1.));
	check_G_validity();
	
	A.init_identity(nb_groups);
	At.init_identity(nb_groups);
	
	// create constraint matrix acting on group variables
	Cgrp.init(nb_groups);
	FOR(cl, C.nb_lines()) if (C.nb_coeffs_in_line(cl) != 2) {
	    Cgrp.add_line();
	    FOR(i, C.nb_coeffs_in_line(cl)) {
		Coeff c = C.ith_coeff_in_line(cl, i);
		Cgrp.add_coeff_to_last_line(G.ith_coeff_in_line(c.index, 0).index, G.ith_coeff_in_line(c.index, 0).a * c.a);
	    }
	}
	
	std::sort(Cgrp.M.begin(), Cgrp.M.end(), eqCmp);
	std::vector<std::vector<Coeff> >::iterator last = std::unique(Cgrp.M.begin(), Cgrp.M.end(), eqeq);
	Cgrp.M.erase(last, Cgrp.M.end());
	
	
	FOR(cl, Cgrp.nb_lines()) {
	    std::vector<Coeff> expended;
	    // remove variables that are already lin comb
	    FOR(t, Cgrp.line(cl).size()) {
		index_t lA = Cgrp.line(cl)[t].index;
		FOR(m, index_t(A.line(lA).size())) {
		    expended.push_back(A.line(lA)[m]);
		    expended.back().a *= Cgrp.line(cl)[t].a;
		}
	    }
	    sort_and_compress(expended);
            
	    // empty a column of A by replacing the corresponding terms by linear combination of others
	    std::vector<Coeff> lin;
	    if (expended.size()>0) {
		index_t rm = first_term_of_in_equals_out(expended, lin);
		while (!At.line(rm).empty()) {
		    index_t l = At.line(rm).back().index;
		    double coeff = value(l, rm);
		    FOR(m, lin.size()) {
			add(l, lin[m].index, coeff * lin[m].a);
		    }
		    remove(l, rm);
		}
	    }
	    
	}
	V.resize(nb_groups);
    }

    void MatrixMixedConstrainedSolver::start_full_real_iter() {
	geo_assert(state == SET_LINEAR_CONSTRAINTS);
	state = SET_FULL_REAL_PB;
	
	nlNewContext();
	nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
	nlSolverParameteri(NL_NB_VARIABLES, NLint(A.nb_lines()));
	nlBegin(NL_SYSTEM);
	nlBegin(NL_MATRIX);
    }
		
    void MatrixMixedConstrainedSolver::start_mixed_iter() {
	geo_assert(state == SET_FULL_REAL_PB);
	state = SET_MIXED_PB;
	nlNewContext();
	nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
	nlSolverParameteri(NL_NB_VARIABLES, NLint(A.nb_lines()));
	nlBegin(NL_SYSTEM);
	
	plop("enforce integer equality");
	FOR(i,intvar.size()) if (intvar[i]) {
	    nlSetVariable(G.line(i)[0].index, G.line(i)[0].a * nint(value(i)));
	    nlLockVariable(G.line(i)[0].index);
	}
	nlBegin(NL_MATRIX);
    }
    
    bool MatrixMixedConstrainedSolver::check_constraints() {
	GEO::Logger::out("HexDom")  << "Check that constraints are really respected..." <<  std::endl;
	FOR(cl, C.nb_lines()) {
	    double sum = 0;
	    FOR(i, C.nb_coeffs_in_line(cl)) {
		Coeff c = C.ith_coeff_in_line(cl, i);
		sum += c.a*value(c.index);
	    }
	    if (std::abs(sum) > .001) {
		std::cerr << "ERROR: condition violated\n";
		show_eq(C.line(cl));
		return false;
	    }
	}
	GEO::Logger::out("HexDom")  << "OK" <<  std::endl;
	return true;
    }
    
    
}
