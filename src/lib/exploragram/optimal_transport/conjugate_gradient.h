
/*
 *  OGF/Graphite: Geometry and Graphics Programming Library + Utilities
 *  Copyright (C) 2000-2009 INRIA - Project ALICE
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
 *  Contact: Bruno Levy - levy@loria.fr
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
 * As an exception to the GPL, Graphite can be linked with 
 *     the following (non-GPL) libraries:
 *     Qt, SuperLU, WildMagic and CGAL
 */

#ifndef H_EXPLORAGRAM_OPTIMAL_TRANSPORT_CONJUGATE_GRADIENT_H
#define H_EXPLORAGRAM_OPTIMAL_TRANSPORT_CONJUGATE_GRADIENT_H

#include <exploragram/basic/common.h>
#include <geogram/basic/numeric.h>

/**
 * \file exploragram/optimal_transport/conjugate_gradient.h
 * \brief a simple in-place conjugate gradient solver for 
 *  implementing Newton iterations for optimal transport.
 */

namespace GEO {

    /**
     * \brief Solves a linear system using Jacobi-preconditioned
     *   conjugate gradient.
     * \param[in] n dimension of the system
     * \param[in] rowptr , colind , val the matrix of the system in the 
     *   compressed row storage format
     * \param[in] rhs the right-hand side of the system
     * \param[out] x the solution of the linear system
     * \param[in] max_iter maximum number of iterations
     * \param[in] threshold the maximum value of 
     *    \f$ \| Ax - b \| / \| b \| \f$ before the iterative solver is stopped.
     * \param[in] symmetric_storage true if the matrix is stored with symmetric
     *    storage (i.e. only lower triangular part is stored), false otherwise
     *    (i.e. then all the coefficients of the matrix are stored).
     */
    void EXPLORAGRAM_API solve_conjugate_gradient(
        index_t  n,
        index_t* rowptr,
        index_t* colind,
        double*  val,
        double*  rhs,
        double* x,
        index_t  max_iter,
        double   threshold,
        bool     symmetric_storage=false
    );
    
}

#endif

