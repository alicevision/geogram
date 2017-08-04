/*
 *  Copyright (c) 2004-2010, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     levy@loria.fr
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 */


#ifndef OPENNL_ITERATIVE_SOLVERS_H
#define OPENNL_ITERATIVE_SOLVERS_H

#include "nl_private.h"

/**
 * \file geogram/NL/nl_iterative_solvers.h
 * \brief Internal OpenNL functions that implement iterative solvers.
 */

/**
 * \brief Solves the system in nlCurrentContext 
 *  using the Conjugate Gradient solver
 * \details The implementation is inspired by 
 * the lsolver library, by Christian Badura, available from:
 * http://www.mathematik.uni-freiburg.de
 * /IAM/Research/projectskr/lin_solver/
 * \return the used number of iterations
 */
NLuint nlSolve_CG(void);

/**
 * \brief Solves the system in nlCurrentContext 
 *  using the preconditioned Conjugate Gradient solver
 * \details The implementation is inspired by 
 * the lsolver library, by Christian Badura, available from:
 * http://www.mathematik.uni-freiburg.de
 * /IAM/Research/projectskr/lin_solver/
 * \return the used number of iterations
 */
NLuint nlSolve_CG_precond(void);

/**
 * \brief Solves the system in nlCurrentContext 
 *  using the stabilized bi conjugate gradient solver.
 * \details The implementation is inspired by 
 * the lsolver library, by Christian Badura, available from:
 * http://www.mathematik.uni-freiburg.de
 * /IAM/Research/projectskr/lin_solver/
 * \return the used number of iterations
 */
NLuint nlSolve_BICGSTAB(void);

/**
 * \brief Solves the system in nlCurrentContext 
 *  using the preconditioned stabilized 
 *  bi conjugate gradient solver.
 * \details The implementation is inspired by 
 * the lsolver library, by Christian Badura, available from:
 * http://www.mathematik.uni-freiburg.de
 * /IAM/Research/projectskr/lin_solver/
 * \return the used number of iterations
 */
NLuint nlSolve_BICGSTAB_precond(void);

/**
 * \brief Solves the system in nlCurrentContext 
 *  using the GMRES solver.
 * \details The implementation is inspired by 
 * the lsolver library, by Christian Badura, available from:
 * http://www.mathematik.uni-freiburg.de
 * /IAM/Research/projectskr/lin_solver/
 * \return the used number of iterations
 */
NLuint nlSolve_GMRES(void);

#endif

