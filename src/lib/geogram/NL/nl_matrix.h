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

#include "nl_private.h"

#ifndef OPENNL_MATRIX_H
#define OPENNL_MATRIX_H

/**
 * \file geogram/NL/nl_matrix.h
 * \brief Internal OpenNL functions to manipulate sparse matrices.
 */

#ifdef __cplusplus
extern "C" {
#endif

/************************************************************************************/
/* Dynamic arrays for sparse row/columns */

/**
 * \brief Represents a coefficient in a sparse matrix
 * \relates NLSparseMatrix
 */
typedef struct  {
    /**
     * \brief index of the coefficient.
     */    
    NLuint index ;

    /**
     * \brief value of the coefficient. 
     */    
    NLdouble value ; 
} NLCoeff ;

/**
 * \brief Represents a row or a column of a sparse matrix
 * \relates NLSparseMatrix
 */
typedef struct {
    /**
     * \brief number of coefficients. 
     */    
    NLuint size ;
    
    /** 
     * \brief number of coefficients that can be 
     * stored without reallocating memory.
     */    
    NLuint capacity ;

    /**
     * \brief the array of coefficients, with enough
     * space to store capacity coefficients.
     */
    NLCoeff* coeff ;  
} NLRowColumn ;

/**
 * \brief Constructs a new NLRowColumn
 * \param[in,out] c a pointer to an 
 *  uninitialized NLRowColumn
 * \relates NLRowColumn
 */
void nlRowColumnConstruct(NLRowColumn* c) ;

/**
 * \brief Destroys a NLRowColumn
 * \details Only the memory allocated by the 
 *  NLRowColumn is freed. The NLRowColumn structure
 *  is not freed.
 * \param[in,out] c a pointer to an NLRowColumn
 * \relates NLRowColumn
 */
void nlRowColumnDestroy(NLRowColumn* c) ;

/**
 * \brief Allocates additional storage for
 *  the coefficients of an NLRowColumn
 * \details Operates like the class vector of the C++ standard library, 
 *  by doubling the capacity each time it is needed. This amortizes
 *  the cost of the growing operations as compared to re-allocating
 *  at each element insertion
 * \param[in,out] c a pointer to an NLRowColumn
 * \relates NLRowColumn
 */
void nlRowColumnGrow(NLRowColumn* c) ;

/**
 * \brief Adds a coefficient to an NLRowColumn.    
 * \details Performs the following operation:
 *  \f$ a_i \leftarrow a_i + value \f$. If the NLRowColumn
 *  already has a coefficient with index \p index, then 
 *  the value is added to that coefficient, else a new
 *  coefficient is created. Additional storage is allocated
 *  as need be.
 * \param[in,out] c a pointer to an NLRowColumn
 * \param[in] index index of the coefficient
 * \param[in] value value of the coefficient
 * \relates NLRowColumn
 */
void nlRowColumnAdd(NLRowColumn* c, NLuint index, NLdouble value) ;

/**
 * \brief Appends a coefficient to an NLRowColumn    .
 * \details In contrast with nlRowColumnAdd(), this function does
 *  not tests whether a coefficient with index \p index already exists
 *  in the NLRowColumn. A new coefficient is always created.
 * \param[in,out] c a pointer to an NLRowColumn
 * \param[in] index index of the coefficient
 * \param[in] value value of the coefficient
 * \relates NLRowColumn
 */
void nlRowColumnAppend(NLRowColumn* c, NLuint index, NLdouble value) ;

/**
 * \brief Zeroes an NLRowColumn.
 * \details No memory is deallocated, the capacity remains
 *  the same.
 * \param[in,out] c a pointer to an NLRowColumn
 * \relates NLRowColumn
 */
void nlRowColumnZero(NLRowColumn* c) ;

/**
 * \brief Zeroes an NLRowColumn and deallocates the memory
 *  used by the NLRowColumn.
 * \details On exit, capacity is zeroed
 * \param[in,out] c a pointer to an NLRowColumn
 * \relates NLRowColumn
 */
void nlRowColumnClear(NLRowColumn* c) ;

/**
 * \brief Sorts the coefficients of an NLRowColumn
 *  by increasing index
 * \param[in,out] c a pointer to an NLRowColumn
 * \relates NLRowColumn
 */
void nlRowColumnSort(NLRowColumn* c) ;

/************************************************************************************/
/* Compressed Row Storage */

/**
 * \brief A compact self-contained storage for 
 *  sparse matrices.
 * \details Unlike with NLSparseMatrix, it is not possible
 *  to add new coefficients in an NLCRSMatrix.
 */
typedef struct {
    /**
     * \brief number of rows. 
     */    
    NLuint m;
    
    /**
     * \brief number of columns. 
     */    
    NLuint n;

    /**
     * \brief array of coefficient values, 
     * size = NNZ (number of non-zero coefficients)
     */    
    NLdouble* val;    

    /**
     * \brief row pointers, size = m+1 
     */    
    NLuint* rowptr;

    /**
     * \brief column indices, size = NNZ 
     */    
    NLuint* colind;

    /**
     * \brief number of slices, used by parallel spMv
     */    
    NLuint nslices;

    /** 
     * \brief slice pointers, size = nslices + 1, 
     * used by parallel spMv 
     */    
    NLuint* sliceptr; 
} NLCRSMatrix;

/**
 * \brief Constructs a new NLCRSMatrix
 * \param[in,out] M pointer to an uninitialized NLCRSMatrix
 * \param[in] m number of rows
 * \param[in] n number of columns
 * \param[in] nnz number of non-zero coefficientsz
 * \param[in] nslices number of slices, used by parallel spMv
 *  (typically, nslices = number of cores)
 * \relates NLCRSMatrix
 */
void nlCRSMatrixConstruct(
    NLCRSMatrix* M, NLuint m, NLuint n, NLuint nnz, NLuint nslices
);

/**
 * \brief Destroys a NLCRSMatrix
 * \details Only the memory allocated by the NLCRSMatrix is freed,
 *  The NLCRSMatrix structure is not freed.
 * \param[in,out] M pointer to an NLCRSMatrix
 * \relates NLCRSMatrix
 */
void nlCRSMatrixDestroy(NLCRSMatrix* M);

/**
 * \brief Loads a NLCRSMatrix from a file
 * \param[out] M a pointer to an uninitialized NLCRSMatriix 
 * \param[in] filename the name of the file
 * \retval NL_TRUE on success
 * \retval NL_FALSE on error
 * \relates NLCRSMatrix
 */
NLboolean nlCRSMatrixLoad(NLCRSMatrix* M, const char* filename);

/**
 * \brief Saves a NLCRSMatrix into a file
 * \param[in] M a pointer to the NLCRSMatriix 
 * \param[in] filename the name of the file
 * \retval NL_TRUE on success
 * \retval NL_FALSE on error
 * \relates NLCRSMatrix
 */
NLboolean nlCRSMatrixSave(NLCRSMatrix* M, const char* filename);

/************************************************************************************/
/* SparseMatrix data structure */

/**
 * for NLSparseMatrix storage: indicates that rows are stored.
 * \relates NLSparseMatrix
 */
#define NL_MATRIX_STORE_ROWS          1

/**
 * for NLSparseMatrix storage: indicates that columns are stored.
 * \relates NLSparseMatrix
 */
#define NL_MATRIX_STORE_COLUMNS       2

/**
 * for NLSparseMatrix storage: indicates that symmetric storage
 * is used (only the lower triangular part is stored).
 * \relates NLSparseMatrix
 */
#define NL_MATRIX_STORE_SYMMETRIC     4

/**
 * for NLSparseMatrix storage: indicates that the matrix
 * is compressed into an NLCRSMatrix.
 * \relates NLSparseMatrix
 */
#define NL_MATRIX_STORE_COMPRESSED    8

/**
 * for NLSparseMatrix storage: indicates that the inverse
 * of the diagonal is stored.
 * \relates NLSparseMatrix
 */
#define NL_MATRIX_STORE_DIAG_INV      16
    
typedef struct {
    /**
     * \brief number of rows 
     */    
    NLuint m ;
    
    /**
     * \brief number of columns 
     */    
    NLuint n ;

    /**
     * \brief number of elements in the diagonal 
     */    
    NLuint diag_size ;

    /**
     * \brief indicates what is stored in this matrix 
     */    
    NLenum storage ;

    /**
     * \brief the rows if (storage & NL_MATRIX_STORE_ROWS), size = m,
     * NULL otherwise
     */     
    NLRowColumn* row ;

    /** 
     * \brief the columns if (storage & NL_MATRIX_STORE_COLUMNS), size = n,
     * NULL otherwise
     */     
    NLRowColumn* column ;

    /**
     * \brief the diagonal elements, size = diag_size 
     */     
    NLdouble*    diag ;

    /**
     * \brief the inverse of the diagonal if (storage & NL_MATRIX_STORE_DIAG_INV), 
     * size = diag_size, NULL otherwise 
     */    
    NLdouble*    diag_inv ;

    /**
     * \brief the compressed CRS representation 
     * if (storage & NL_MATRIX_STORE_COMPRESSED),
     * NULL otherwise
     */    
    NLCRSMatrix* compressed ;  
} NLSparseMatrix ;


/**
 * \brief Constructs a new NLSparseMatrix
 * \param[in,out] M a pointer to an uninitialized NLSparseMatrix
 * \param[in] m number of rows
 * \param[in] n number of columns
 * \param[in] storage a bitwise or combination of flags that
 *  indicate what needs to be stored in the matrix.
 * \relates NLSparseMatrix
 */
void nlSparseMatrixConstruct(
    NLSparseMatrix* M, NLuint m, NLuint n, NLenum storage
) ;

/**
 * \brief Destroys an NLSparseMatrix
 * \details Only the memory allocated by the NLSparseMatrix
 *  is freed. The NLSparseMatrix structure is not freed.
 * \param[in,out] M a pointer to an NLSparseMatrix
 * \relates NLSparseMatrix
 */
void nlSparseMatrixDestroy(NLSparseMatrix* M) ;

/**
 * \brief Adds a coefficient to an NLSparseMatrix
 * \details Performs the following operation:
 *  \$ a_{i,j} \leftarrow a_{i,j} + \mbox{value} \$
 * \param[in,out] M a pointer to an NLSparseMatrix
 * \param[in] i index of the row
 * \param[in] j index of the column
 * \param[in] value the coefficient to be added
 * \relates NLSparseMatrix
 */
void nlSparseMatrixAdd(
    NLSparseMatrix* M, NLuint i, NLuint j, NLdouble value
) ;

/**
 * \brief Zeroes an NLSparseMatrix
 * \details The memory is not freed.
 * \param[in,out] M a pointer to the NLSparseMatrix to zero
 * \relates NLSparseMatrix
 */
void nlSparseMatrixZero( NLSparseMatrix* M) ;

/**
 * \brief Clears an NLSparseMatrix
 * \details The memory is freed.
 * \param[in,out] M a pointer to the NLSparseMatrix to zero
 * \relates NLSparseMatrix
 */
void nlSparseMatrixClear( NLSparseMatrix* M) ;

/**
 * \brief Gets the number of non-zero coefficient
 *  in an NLSparseMatrix
 * \param[in] M a pointer to the NLSparseMatrix
 * \return the number of non-zero coefficients in \p M
 * \relates NLSparseMatrix
 */
NLuint nlSparseMatrixNNZ( NLSparseMatrix* M) ;

/**
 * \brief Sorts the coefficients in an NLSParseMatrix
 * \param[in,out] M a pointer to the NLSparseMatrix 
 * \relates NLSparseMatrix
 */
void nlSparseMatrixSort( NLSparseMatrix* M) ;

/**
 * \brief Computes the inverse of the diagonal and
 *  stores it into an NLSparseMatrix.
 * \details This is used for instance by the Jacobi
 *  pre-conditioner.
 * \param[in,out] M a pointer to the NLSparseMatrix 
 * \relates NLSparseMatrix
 */
void nlSparseMatrixComputeDiagInv( NLSparseMatrix* M);

/**
 * \brief Computes a Compressed Row Storage representation
 *  of a NLSparseMatrix.
 * \details Once this function is called, nlSparseMatrixMult()
 *  is faster by an order of magnitude
 * \param[in,out] M a pointer to the NLSparseMatrix 
 * \relates NLSparseMatrix
 */
void nlSparseMatrixCompress( NLSparseMatrix* M);
    
/************************************************************************************/
/* SparseMatrix x Vector routine */

/**
 * \brief Computes a matrix-vector product
 * \param[in] A a pointer to the matrix
 * \param[in] x the vector to be multiplied, size = A->n
 * \param[in] y where to store the result, size = A->m
 * \relates NLSparseMatrix
 */
void nlSparseMatrixMult(NLSparseMatrix* A, const NLdouble* x, NLdouble* y) ;

#ifdef __cplusplus
}
#endif

#endif
