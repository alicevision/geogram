/*
 *  Copyright (c) 2012-2014, Bruno Levy
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
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#ifndef GEOGRAM_NUMERICS_SPARSE_MATRIX
#define GEOGRAM_NUMERICS_SPARSE_MATRIX

#include <geogram/basic/common.h>
#include <geogram/NL/nl_matrix.h>

/**
 * \file geogram/numerics/sparse_matrix.h
 * \brief C++ API around OpenNL sparse matrices
 */

namespace GEO {

    /**
     * \brief A C++ shell around OpenNL sparse matrix
     *  datastructure.
     */
    class GEOGRAM_API SparseMatrix {
    public:

        /**
         * \brief SparseMatrix default constructor.
         * \details Creates an uninitialized SparseMatrix.
         */
        SparseMatrix() {
            impl_.m = 0;
            impl_.n = 0;
        }
        
        /**
         * \brief SparseMatrix constructor.
         * \param[in] m number of rows
         * \param[in] n number of columns
         * \param[in] storage an '|' combination of 
         *  NL_MATRIX_STORE_ROWS, NL_MATRIX_STORE_COLUMNS, 
         *  NL_MATRIX_STORE_SYMMETRIC 
         */
        SparseMatrix(
            index_t m, index_t n, NLenum storage = NL_MATRIX_STORE_ROWS
        ) {
            nlSparseMatrixConstruct(&impl_,m,n,storage);
        }

        /**
         * \brief SparseMatrix destructor.
         */
        ~SparseMatrix() {
            clear();
        }

        /**
         * \brief Tests whether this SparseMatrix is initialized.
         * \retval true if this SparseMatrix is initialized
         * \retval false otherwise
         */
        bool initialized() const {
            geo_debug_assert((impl_.m == 0)==(impl_.n == 0));
            return (impl_.m != 0);
        }

        /**
         * \brief Clears this SparseMatrix.
         * \details On exit, this SparseMatrix is in the uninitialized
         *  state. If this SparseMatrix was uninitialized, then this function
         *  does nothing.
         */
        void clear() {
            if(initialized()) {
                sparseMatrixDestroy(&impl_);
                impl_.m = 0;
                impl_.n = 0;
            }
        }

        /**
         * \brief Initializes this SparseMatrix.
         * \details It is allowed to call this function on an already
         *  initialized SparseMatrix (then it will be cleared).
         */
        void initialize(
            index_t m, index_t n, NLenum storage = NL_MATRIX_STORE_ROWS
        ) {
            clear();
            sparseMatrixConstruct(&impl_, m, n, storage);
        }


        /**
         * \brief Tests whether a SparseMatrix is locked.
         * \details A locked SparseMatrix can no longer be modified.
         */  
        bool locked() const {
            return (impl_.storage & NL_MATRIX_STORE_COMPRESSED != 0);
        }

        /**
         * \brief Compresses this SparseMatrix.
         * \details When a SparseMatrix is compressed, matrix x vector 
         *  product can be computed more quickly, but it is no longer
         *  possible to modify the matrix. The matrix becomes locked.
         */
        void compress() {
            nlSparseMatrixCompress(&impl_);
        }
        
        /**
         * \brief Gets the number of rows.
         * \return the number of rows of this SparseMatrix, or 0 if it
         *  is not initialized.
         */
        index_t m() const {
            return impl_.m;
        }

        /**
         * \brief Gets the number of columns.
         * \return the number of columns of this SparseMatrix, or 0 if it
         *  is not initialized.
         */
        index_t n() const {
            return impl_.n;            
        }

        /**
         * \brief Adds a scalar to a coefficient of the matrix.
         * \details \$ m_{i,j} \leftarrow m_{i,j} + val \$
         * \param[in] i , j row and column indices
         * \param[in] val the value to be added
         */
        void add(index_t i, index_t j, double val) {
            geo_debug_assert(i < m());
            geo_debug_assert(j < n());
            geo_debug_assert(!locked());
            nlSparseMatrixAdd(&impl_, i, j, val);
        }

        /**
         * \brief Computes the product between this SparseMatrix and a vector.
         * \param[in] x a pointer to an array of n() doubles with the vector 
         *  to be multiplied
         * \param[out] y a pointer to an array of m() doubles where to store
         *  the result vector
         */
        void mult(const double* x, double* y) const {
            nlSparseMatrixMult(&impl_, x, y);
        }
        
    private:
        /**
         * \brief forbids copy.
         */
        SparseMatrix(const SparseMatrix& rhs);
        
        /**
         * \brief forbids copy.
         */
        SparseMatrix& operator=(const SparseMatrix& rhs);
    private:
        NLSparseMatrix impl_;
    };
    
}

#endif

