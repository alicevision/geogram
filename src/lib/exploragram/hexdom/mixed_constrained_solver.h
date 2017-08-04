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

#ifndef H_HEXDOM_ALGO_MIXED_CONSTRAINED_SOLVER_H
#define H_HEXDOM_ALGO_MIXED_CONSTRAINED_SOLVER_H

#include <exploragram/basic/common.h>
#include <exploragram/hexdom/basic.h>
#include <exploragram/hexdom/id_map.h>
#include <geogram/NL/nl.h>
#include <algorithm>

namespace GEO{

    struct Coeff {
        Coeff(index_t p_index, double p_a) :
	    index(p_index), a(p_a) {
	}
	    
        Coeff() : index(index_t(-1)), a(-1.0) {
	}

	index_t index;	    
	double a;
    };
	
    struct EXPLORAGRAM_API NaiveSparseMatrix{
	NaiveSparseMatrix(index_t p_nb_col = index_t(-1), index_t p_nb_lines = 0){ init(p_nb_col, p_nb_lines); }
	void init(index_t p_nb_col, index_t p_nb_lines = 0)        { nb_col = p_nb_col; M.resize(p_nb_lines); }
	void init_identity(index_t size)                           { init(size, size); FOR(i, size) add_coeff_to_line(i, i, 1); }
	void add_line()                                            { M.resize(M.size() + 1); }
	void add_coeff_to_line(index_t j, index_t l, double coeff) {
	    geo_debug_assert(l < nb_lines()); geo_debug_assert(j<nb_col); if (coeff != 0.0) M[l].push_back(Coeff(j, coeff));
	}
	void add_coeff_to_last_line(index_t j, double coeff)       { geo_debug_assert(j < nb_col); if (coeff != 0.0) M.back().push_back(Coeff(j, coeff)); }
	Coeff& ith_coeff_in_line(index_t l, index_t i)             { return M[l][i]; }
	index_t nb_coeffs_in_line(index_t l)                       { return index_t(M[l].size()); }

	void show_dense(const char* header = " show matrix \n");
	void show_sparse(const char* header = " show matrix \n",index_t nb_min_terms = 0);

	std::vector<Coeff>& line(index_t l) { return M[l]; }
	index_t nb_lines() { return index_t(M.size()); }
	index_t nb_col;
	std::vector<std::vector<Coeff> > M;
    };


    struct EXPLORAGRAM_API MatrixMixedConstrainedSolver{
	// I do not check the construction of every single equation, but only the global organisation
	enum State{
	    SET_INTEGER_CONSTRAINTS,
	    SET_LINEAR_CONSTRAINTS,
	    SET_FULL_REAL_PB,
	    SET_MIXED_PB
	} state;

	MatrixMixedConstrainedSolver(index_t p_nb_vars);
		
	void is_integer(index_t id)             { geo_assert(state == SET_INTEGER_CONSTRAINTS); intvar[id] = true; }
	void start_linear_contraints()  { geo_assert(state == SET_INTEGER_CONSTRAINTS); state = SET_LINEAR_CONSTRAINTS; }
	void begin_constraint()                 { geo_assert(state == SET_LINEAR_CONSTRAINTS); C.add_line(); }
	void add_constraint_coeff(index_t id, double coeff){ geo_assert(state == SET_LINEAR_CONSTRAINTS); if (coeff!=0.0) C.add_coeff_to_last_line(id, coeff); }
	void end_constraint()                   { geo_assert(state == SET_LINEAR_CONSTRAINTS); if (C.M.back().empty()) C.M.pop_back(); }

	void remove_term(std::vector<Coeff> &eq, index_t i){ eq[i] = eq.back(); eq.pop_back(); }

	index_t first_term_of_in_equals_out(std::vector<Coeff> &in, std::vector<Coeff> &out);

	void sort_and_compress(std::vector<Coeff> &eq);

	void remove_term_sorted(std::vector<Coeff> &eq, index_t ind){
	    index_t i = 0;
	    while (eq[i].index != ind) { if (i >= eq.size()) abort(); i++; }
	    while (i + 1 < eq.size()){ eq[i] = eq[i + 1]; i++; }
	    geo_assert(!eq.empty());
	    eq.pop_back();
	}
	void create_term_sorted(std::vector<Coeff> &eq, index_t ind, double coeff){
	    Coeff c; c.index = ind; c.a = coeff;
	    eq.push_back(c);
	    int  i = int(eq.size()) - 1;
	    while (i > 0 && ind<eq[index_t(i - 1)].index) { std::swap(eq[index_t(i)], eq[index_t(i - 1)]); i--; }
	}

	void add_to_term(std::vector<Coeff> &eq, index_t ind, double coeff){
	    if (coeff == 0) return;
	    for (size_t i = 0; i < eq.size(); i++) if (eq[i].index == ind){
		    eq[i].a += coeff;
		    if (std::fabs(eq[i].a) < .0001)
			remove_term_sorted(eq, ind);
		    return;
		}
	    create_term_sorted(eq, ind, coeff);
	}


	void remove(index_t i, index_t j){
	    remove_term_sorted(A.line(i), j);
	    remove_term_sorted(At.line(j), i);
	}


	double value(index_t i, index_t j){
	    FOR(m,A.line(i).size()) if (A.line(i)[m].index == j) return A.line(i)[m].a;
	    plop("access to the value of a zero coefficient");
	    return 0;
	}

	void add(index_t i, index_t j, double coeff){
	    if (coeff == 0) return;
	    add_to_term(A.line(i), j, coeff);
	    add_to_term(At.line(j), i, coeff);
	}

	void check_G_validity();

	void end_linear_constraints();

	void start_full_real_iter();
		
	void start_mixed_iter();

	void start_iter(index_t it){ if (it == 0) start_full_real_iter(); else start_mixed_iter(); }

	void begin_energy(){ nlBegin(NL_ROW); }
	void add_energy_coeff(index_t id, double coeff){ 
	    if (coeff == 0) return;
	    index_t grp = G.line(id)[0].index;
	    double sign = G.line(id)[0].a;
	    FOR(m, index_t(A.line(grp).size())){
		Coeff t = A.ith_coeff_in_line(grp, m);
		nlCoefficient(t.index, sign * t.a*coeff);
	    }
	}
	void add_energy_rhs(double rhs){ nlRightHandSide(rhs); }
	void end_energy(){ nlEnd(NL_ROW); }


	void end_iter(){
	    nlEnd(NL_MATRIX);
	    nlEnd(NL_SYSTEM);
	    nlSolve();
	    FOR(i, A.nb_lines()) V[i] = nlGetVariable(i);
	    nlDeleteContext(nlGetCurrent());
	}
	void end_full_real_iter(){ end_iter(); }
	void end_mixed_iter(){ end_iter(); }


	double value(index_t i){
	    index_t grp = G.line(i)[0].index;
	    double res = 0;
	    FOR(m,A.line(grp).size()){
		Coeff t = A.ith_coeff_in_line(grp,m);
		res += V[t.index] * t.a;
	    }
	    return G.line(i)[0].a * res;
	}

	bool check_constraints();
                
	index_t nb_vars;

	NaiveSparseMatrix G;                    // matrix that gives the group and sign from each user variable
	NaiveSparseMatrix A;                    // matrix that gives group variables from the inner solver variables such that Cx=0
	NaiveSparseMatrix At;                   // transposed of A
	NaiveSparseMatrix C;                    // constraints matrix solution X of the problem must respect CX = 0;
                
	NaiveSparseMatrix Cgrp;
	std::vector<bool> intvar;               // which user variables are constrainted to be integer
	std::vector<double> V;                  // inner variables


	std::vector<int> nullgrp;               // groups that must by set to null do to opposite constraints
    };


}

#endif
