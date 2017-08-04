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
 *
 */

#ifndef OPENNL_PRECONDITIONERS_H
#define OPENNL_PRECONDITIONERS_H

#include "nl_private.h"

/**
 * \file geogram/NL/nl_preconditioners.h
 * \brief Internal OpenNL functions that implement preconditioners.
 */

/******************************************************************************/
/* preconditioners */

/**
 * \brief Computes the product between the Jacobi 
 *  preconditioner and a vector.
 * \details The Jacobi preconditioner corresponds to
 *  the inverse of the diagonal of the matrix stored
 *  in the current OpenNL context.
 * \param[in] x the vector to be multiplied, dimension = nlCurrentContext->n,
 *  remains unchanged
 * \param[out] y where to store the result, dimension = nlCurrentContext->n
 */
void nlPreconditioner_Jacobi(const NLdouble* x, NLdouble* y) ;

/**
 * \brief Computes the product between the SSOR
 *  preconditioner and a vector.
 * \details The SSOR preconditioner is computed from the
 *  matrix stored the current OpenNL context and the omega
 *  parameter, also stored in the current OpenNL context.
 * \param[in] x the vector to be multiplied, dimension = nlCurrentContext->n,
 *  remains unchanged
 * \param[out] y where to store the result, dimension = nlCurrentContext->n
 */
void nlPreconditioner_SSOR(const NLdouble* x, NLdouble* y) ;


/**
 * \brief Multiplies a vector by the diagonal of the matrix
 *  stored in the current OpenNL context.
 * \details \$ x \leftarrow 1/\omega \mbox{diag}(M) x \$,
 *   used to implement the SSOR preconditioner.
 * \param[in,out] x the vector to be multiplied, 
 *  size = nlCurrentContext->n
 * \param[in] omega all components are divided by omega
 */
void nlMultDiagonal(NLdouble* x, NLdouble omega) ;

/**
 * \brief Multiplies a vector by the inverse of the
 *  diagonal of the matrix stored in the current OpenNL context.
 * \details \$ x \leftarrow \omega \mbox{diag}(M)^{-1} x \$,
 *   used to implement the SSOR preconditioner.
 * \param[in,out] x the vector to be multiplied,
 *  size = nlCurrentContext->n
 * \param[in] omega all components are multiplied by omega
 */
void nlMultDiagonalInverse(NLdouble* x, NLdouble omega) ;

/**
 * \brief Multiplies a vector by the inverse of the
 *  lower triangular block of the matrix stored in 
 *  the current OpenNL context.
 * \details \$ x \leftarrow \omega \mbox{trilow}(M)^{-1} x \$,
 *   used to implement the SSOR preconditioner.
 * \param[in] x the vector to be multiplied,
 *  size = nlCurrentContext->n
 * \param[out] y where to store the result,
 *  size = nlCurrentContext->n
 * \param[in] omega all components are multiplied by omega
 */
void nlMultLowerInverse(const NLdouble* x, NLdouble* y, NLdouble omega) ;

/**
 * \brief Multiplies a vector by the inverse of the
 *  upper triangular block of the matrix stored in 
 *  the current OpenNL context.
 * \details \$ x \leftarrow \omega \mbox{triup}(M)^{-1} x \$,
 *   used to implement the SSOR preconditioner.
 * \param[in] x the vector to be multiplied,
 *  size = nlCurrentContext->n
 * \param[out] y where to store the result,
 *  size = nlCurrentContext->n
 * \param[in] omega all components are multiplied by omega
 */
void nlMultUpperInverse(const NLdouble* x, NLdouble* y, NLdouble omega) ;

#endif
