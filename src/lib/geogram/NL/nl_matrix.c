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

#include "nl_matrix.h"

/*
 Some warnings about const cast in callback for
 qsort() function.
 */

#ifdef __clang__
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif

void nlRowColumnConstruct(NLRowColumn* c) {
    c->size     = 0 ;
    c->capacity = 0 ;
    c->coeff    = NULL ;
}

void nlRowColumnDestroy(NLRowColumn* c) {
    NL_DELETE_ARRAY(c->coeff) ;
#ifdef NL_PARANOID
    NL_CLEAR(NLRowColumn, c) ; 
#endif
}

void nlRowColumnGrow(NLRowColumn* c) {
    if(c->capacity != 0) {
        c->capacity = 2 * c->capacity ;
        c->coeff = NL_RENEW_ARRAY(NLCoeff, c->coeff, c->capacity) ;
    } else {
        c->capacity = 4 ;
        c->coeff = NL_NEW_ARRAY(NLCoeff, c->capacity) ;
    }
}

void nlRowColumnAdd(NLRowColumn* c, NLuint index, NLdouble value) {
    NLuint i ;
    for(i=0; i<c->size; i++) {
        if(c->coeff[i].index == index) {
            c->coeff[i].value += value ;
            return ;
        }
    }
    if(c->size == c->capacity) {
        nlRowColumnGrow(c) ;
    }
    c->coeff[c->size].index = index ;
    c->coeff[c->size].value = value ;
    c->size++ ;
}

/* Does not check whether the index already exists */
void nlRowColumnAppend(NLRowColumn* c, NLuint index, NLdouble value) {
    if(c->size == c->capacity) {
        nlRowColumnGrow(c) ;
    }
    c->coeff[c->size].index = index ;
    c->coeff[c->size].value = value ;
    c->size++ ;
}

void nlRowColumnZero(NLRowColumn* c) {
    c->size = 0 ;
}

void nlRowColumnClear(NLRowColumn* c) {
    c->size     = 0 ;
    c->capacity = 0 ;
    NL_DELETE_ARRAY(c->coeff) ;
}

static int nlCoeffCompare(const void* p1, const void* p2) {
    return (((NLCoeff*)(p2))->index < ((NLCoeff*)(p1))->index) ;
}

void nlRowColumnSort(NLRowColumn* c) {
    qsort(c->coeff, c->size, sizeof(NLCoeff), nlCoeffCompare) ;
}

/******************************************************************************/
/* CRSMatrix data structure */

void nlCRSMatrixConstruct(
    NLCRSMatrix* M, NLuint m, NLuint n, NLuint nnz, NLuint nslices
) {
    M->m = m;
    M->n = n;
    M->nslices = nslices;
    M->val = NL_NEW_ARRAY(double, nnz);
    M->rowptr = NL_NEW_ARRAY(NLuint, m+1);
    M->colind = NL_NEW_ARRAY(NLuint, nnz);
    M->sliceptr = NL_NEW_ARRAY(NLuint, nslices+1);
}

void nlCRSMatrixDestroy(NLCRSMatrix* M) {
    NL_DELETE_ARRAY(M->val);
    NL_DELETE_ARRAY(M->rowptr);
    NL_DELETE_ARRAY(M->colind);
    NL_DELETE_ARRAY(M->sliceptr);
    M->m = 0;
    M->n = 0;
    M->nslices = 0;
}


NLboolean nlCRSMatrixSave(NLCRSMatrix* M, const char* filename) {
    NLuint nnz = M->rowptr[M->m];
    FILE* f = fopen(filename, "rb");
    if(f == NULL) {
        nlError("nlCRSMatrixSave", "Could not open file");
        return NL_FALSE;
    }

    fwrite(&M->m, sizeof(NLuint), 1, f);
    fwrite(&M->n, sizeof(NLuint), 1, f);
    fwrite(&nnz, sizeof(NLuint), 1, f);

    fwrite(M->rowptr, sizeof(NLuint), M->m+1, f);
    fwrite(M->colind, sizeof(NLuint), nnz, f);
    fwrite(M->val, sizeof(double), nnz, f);
    
    return NL_TRUE;
}

NLboolean nlCRSMatrixLoad(NLCRSMatrix* M, const char* filename) {
    NLuint nnz = 0;
    FILE* f = fopen(filename, "rb");
    NLboolean truncated = NL_FALSE;
    
    if(f == NULL) {
        nlError("nlCRSMatrixLoad", "Could not open file");
        return NL_FALSE;
    }
    
    truncated = truncated || (
        fread(&M->m, sizeof(NLuint), 1, f) != 1 ||
        fread(&M->n, sizeof(NLuint), 1, f) != 1 ||
        fread(&nnz, sizeof(NLuint), 1, f) != 1
    );

    if(truncated) {
        M->rowptr = NULL;
        M->colind = NULL;
        M->val = NULL;
    } else {
        M->rowptr = NL_NEW_ARRAY(NLuint, M->m+1);
        M->colind = NL_NEW_ARRAY(NLuint, nnz);
        M->val = NL_NEW_ARRAY(double, nnz);
        truncated = truncated || (
            fread(M->rowptr, sizeof(NLuint), M->m+1, f) != M->m+1 ||
            fread(M->colind, sizeof(NLuint), nnz, f) != nnz ||
            fread(M->val, sizeof(double), nnz, f) != nnz
        );
    }

    if(truncated) {
        nlError("nlCRSMatrixSave", "File appears to be truncated");
        NL_DELETE_ARRAY(M->rowptr);
        NL_DELETE_ARRAY(M->colind);
        NL_DELETE_ARRAY(M->val);
        return NL_FALSE;
    } else {
        M->nslices = 1;    
        M->sliceptr = NL_NEW_ARRAY(NLuint, M->nslices+1);
        M->sliceptr[0] = 0;
        M->sliceptr[1] = M->m;
    }

    fclose(f);
    return NL_TRUE;
}

static NLuint nlCRSMatrixNNZ(NLCRSMatrix* M) {
    return M->rowptr[M->m];
}

static void nlCRSMatrixMultSlice(
    NLCRSMatrix* M, const double* x, double* y, NLuint Ibegin, NLuint Iend
) {
    NLuint i,j;
    for(i=Ibegin; i<Iend; ++i) {
        double sum=0.0;
        for(j=M->rowptr[i]; j<M->rowptr[i+1]; ++j) {
            sum += M->val[j] * x[M->colind[j]];
        }
        y[i] = sum; 
    }
}


static void nlCRSMatrixMult(
    NLCRSMatrix* M, const double* x, double* y
) {
    
    int slice;
    int nslices = (int)(M->nslices);
    
#if defined(_OPENMP)
#pragma omp parallel for private(slice)
#endif
    
    for(slice=0; slice<nslices; ++slice) {
        nlCRSMatrixMultSlice(
            M,x,y,M->sliceptr[slice],M->sliceptr[slice+1]
        );
    }
}

/******************************************************************************/
/* SparseMatrix data structure */

void nlSparseMatrixConstruct(
    NLSparseMatrix* M, NLuint m, NLuint n, NLenum storage
) {
    NLuint i ;
    M->m = m ;
    M->n = n ;
    M->storage = storage ;
    if(storage & NL_MATRIX_STORE_ROWS) {
        M->row = NL_NEW_ARRAY(NLRowColumn, m) ;
        for(i=0; i<n; i++) {
            nlRowColumnConstruct(&(M->row[i])) ;
        }
    } else {
        M->row = NULL ;
    }

    if(storage & NL_MATRIX_STORE_COLUMNS) {
        M->column = NL_NEW_ARRAY(NLRowColumn, n) ;
        for(i=0; i<n; i++) {
            nlRowColumnConstruct(&(M->column[i])) ;
        }
    } else {
        M->column = NULL ;
    }

    M->diag_size = MIN(m,n) ;
    M->diag = NL_NEW_ARRAY(NLdouble, M->diag_size) ;
}

static void nlSparseMatrixDestroyRowColumns(NLSparseMatrix* M) {
    NLuint i ;
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i=0; i<M->m; i++) {
            nlRowColumnDestroy(&(M->row[i])) ;
        }
        NL_DELETE_ARRAY(M->row) ;
    }
    M->storage = (NLenum)((int)(M->storage) & ~NL_MATRIX_STORE_ROWS);
    
    if(M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i=0; i<M->n; i++) {
            nlRowColumnDestroy(&(M->column[i])) ;
        }
        NL_DELETE_ARRAY(M->column) ;
    }
    M->storage = (NLenum)((int)(M->storage) & ~NL_MATRIX_STORE_COLUMNS);    
}

void nlSparseMatrixDestroy(NLSparseMatrix* M) {
    nlSparseMatrixDestroyRowColumns(M);
    NL_DELETE_ARRAY(M->diag) ;
    if(M->storage & NL_MATRIX_STORE_DIAG_INV) {
        NL_DELETE_ARRAY(M->diag_inv) ;
        M->storage = (NLenum)((int)(M->storage) & ~NL_MATRIX_STORE_DIAG_INV);
    }    
    if(M->storage & NL_MATRIX_STORE_COMPRESSED) {
        nlCRSMatrixDestroy(M->compressed);
        NL_DELETE(M->compressed);
        M->storage = (NLenum)((int)(M->storage) & ~NL_MATRIX_STORE_COMPRESSED);
    }
#ifdef NL_PARANOID
    NL_CLEAR(NLSparseMatrix,M) ;
#endif
}

void nlSparseMatrixAdd(NLSparseMatrix* M, NLuint i, NLuint j, NLdouble value) {
    nl_parano_range_assert(i, 0, M->m - 1) ;
    nl_parano_range_assert(j, 0, M->n - 1) ;
    nl_debug_assert(!(M->storage & NL_MATRIX_STORE_COMPRESSED));
    if((M->storage & NL_MATRIX_STORE_SYMMETRIC) && (j > i)) {
        return ;
    }
    if(i == j) {
        M->diag[i] += value ;
    }
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        nlRowColumnAdd(&(M->row[i]), j, value) ;
    }
    if(M->storage & NL_MATRIX_STORE_COLUMNS) {
        nlRowColumnAdd(&(M->column[j]), i, value) ;
    }
}

void nlSparseMatrixZero( NLSparseMatrix* M) {
    NLuint i ;
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i=0; i<M->m; i++) {
            nlRowColumnZero(&(M->row[i])) ;
        }
    }
    if(M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i=0; i<M->n; i++) {
            nlRowColumnZero(&(M->column[i])) ;
        }
    }
    NL_CLEAR_ARRAY(NLdouble, M->diag, M->diag_size) ;
    if(M->storage & NL_MATRIX_STORE_DIAG_INV) {
        NL_CLEAR_ARRAY(NLdouble, M->diag_inv, M->diag_size) ;        
    }
}

void nlSparseMatrixClear( NLSparseMatrix* M) {
    NLuint i ;
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i=0; i<M->m; i++) {
            nlRowColumnClear(&(M->row[i])) ;
        }
    }
    if(M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i=0; i<M->n; i++) {
            nlRowColumnClear(&(M->column[i])) ;
        }
    }
    NL_CLEAR_ARRAY(NLdouble, M->diag, M->diag_size) ;
    if(M->storage & NL_MATRIX_STORE_DIAG_INV) {
        NL_CLEAR_ARRAY(NLdouble, M->diag_inv, M->diag_size) ;        
    }
}

/* Returns the number of non-zero coefficients */
NLuint nlSparseMatrixNNZ( NLSparseMatrix* M) {
    NLuint nnz = 0 ;
    NLuint i ;
    if (M->storage & NL_MATRIX_STORE_COMPRESSED) {
        nnz = nlCRSMatrixNNZ(M->compressed);
    } else if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i = 0; i<M->m; i++) {
            nnz += M->row[i].size ;
        }
    } else if (M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i = 0; i<M->n; i++) {
            nnz += M->column[i].size ;
        }
    } else {
        nl_assert_not_reached ;
    }
    return nnz ;
}

void nlSparseMatrixSort( NLSparseMatrix* M) {
    NLuint i ;
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i = 0; i<M->m; i++) {
            nlRowColumnSort(&(M->row[i])) ;                
        }
    } 
    if (M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i = 0; i<M->n; i++) {
            nlRowColumnSort(&(M->column[i])) ;
        }
    } 
}

void nlSparseMatrixComputeDiagInv( NLSparseMatrix* M) {
    NLuint i;
    NLdouble s;
    if(!(M->storage & NL_MATRIX_STORE_DIAG_INV)) {
        M->diag_inv = NL_NEW_ARRAY(double, M->diag_size);
        M->storage |= NL_MATRIX_STORE_DIAG_INV;
    }
    for(i=0; i<M->diag_size; ++i) {
        s = M->diag[i];
        if(s != 0.0) {
            s = 1.0 / s;
        }
        M->diag_inv[i] = s;
    }
}

void nlSparseMatrixCompress( NLSparseMatrix* M) {
    NLuint nnz = nlSparseMatrixNNZ(M);
    NLuint nslices = 8; /* TODO: get number of cores */
    NLuint slice, cur_bound, cur_NNZ, cur_row;
    NLuint i,ij,k;
    NLuint slice_size = nnz / nslices;
    NLCRSMatrix* CRS = NL_NEW(NLCRSMatrix);
    
    nlCRSMatrixConstruct(CRS, M->m, M->n, nnz, nslices);
    
    nl_assert(M->storage & NL_MATRIX_STORE_ROWS);
    nl_assert(!(M->storage & NL_MATRIX_STORE_SYMMETRIC));
    nl_assert(!(M->storage & NL_MATRIX_STORE_COMPRESSED));

    nlSparseMatrixSort(M);

    /* Copy dynamic sparse matrix into compressed row storage */
    k=0;
    for(i=0; i<M->m; ++i) {
        NLRowColumn* Ri = &(M->row[i]) ;
        CRS->rowptr[i] = k;
        for(ij=0; ij<Ri->size; ij++) {
            NLCoeff* c = &(Ri->coeff[ij]) ;
            CRS->val[k] = c->value;
            CRS->colind[k] = c->index;
            ++k;
        }
    }
    CRS->rowptr[M->m] = nnz;

    /* Create "slices" to be used by parallel sparse matrix vector product */
    cur_bound = slice_size;
    cur_NNZ = 0;
    cur_row = 0;
    CRS->sliceptr[0]=0;
    for(slice=1; slice<nslices; ++slice) {
        while(cur_NNZ < cur_bound && cur_row < M->m) {
            ++cur_row;
            cur_NNZ += CRS->rowptr[cur_row+1] - CRS->rowptr[cur_row];
        }
        CRS->sliceptr[slice] = cur_row;
        cur_bound += slice_size;
    }
    CRS->sliceptr[nslices]=M->m; 
    
    M->compressed = CRS;
    M->storage = (M->storage | NL_MATRIX_STORE_COMPRESSED);
    nlSparseMatrixDestroyRowColumns(M);
}



/*****************************************************************************/
/* SparseMatrix x Vector routines, internal helper routines */

static void nlSparseMatrix_mult_rows_symmetric(
        NLSparseMatrix* A,
        const NLdouble* x,
        NLdouble* y
) {
    NLuint m = A->m ;
    NLuint i,ij ;
    NLCoeff* c = NULL ;
    for(i=0; i<m; i++) {
        NLRowColumn* Ri = &(A->row[i]) ;
        y[i] = 0 ;
        for(ij=0; ij<Ri->size; ij++) {
            c = &(Ri->coeff[ij]) ;
            y[i] += c->value * x[c->index] ;
            if(i != c->index) {
                y[c->index] += c->value * x[i] ;
            }
        }
    }
}

static void nlSparseMatrix_mult_rows(
        NLSparseMatrix* A,
        const NLdouble* x,
        NLdouble* y
) {
    /* 
     * Note: OpenMP does not like unsigned ints
     * (causes some floating point exceptions),
     * therefore I use here signed ints for all
     * indices.
     */
    
    int m = (int)(A->m) ;
    int i,ij ;
    NLCoeff* c = NULL ;
    NLRowColumn* Ri = NULL;

#if defined(_OPENMP)    
#pragma omp parallel for private(i,ij,c,Ri)
#endif
    
    for(i=0; i<m; i++) {
        Ri = &(A->row[i]) ;       
        y[i] = 0 ;
        for(ij=0; ij<(int)(Ri->size); ij++) {
            c = &(Ri->coeff[ij]) ;
            y[i] += c->value * x[c->index] ;
        }
    }
}

static void nlSparseMatrix_mult_cols_symmetric(
        NLSparseMatrix* A,
        const NLdouble* x,
        NLdouble* y
) {
    NLuint n = A->n ;
    NLuint j,ii ;
    NLCoeff* c = NULL ;
    for(j=0; j<n; j++) {
        NLRowColumn* Cj = &(A->column[j]) ;       
        y[j] = 0 ;
        for(ii=0; ii<Cj->size; ii++) {
            c = &(Cj->coeff[ii]) ;
            y[c->index] += c->value * x[j] ;
            if(j != c->index) {
                y[j] += c->value * x[c->index] ;
            }
        }
    }
}

static void nlSparseMatrix_mult_cols(
        NLSparseMatrix* A,
        const NLdouble* x,
        NLdouble* y
) {
    NLuint n = A->n ;
    NLuint j,ii ; 
    NLCoeff* c = NULL ;
    NL_CLEAR_ARRAY(NLdouble, y, A->m) ;
    for(j=0; j<n; j++) {
        NLRowColumn* Cj = &(A->column[j]) ;
        for(ii=0; ii<Cj->size; ii++) {
            c = &(Cj->coeff[ii]) ;
            y[c->index] += c->value * x[j] ;
        }
    }
}

/*****************************************************************************/
/* SparseMatrix x Vector routines, main driver routine */

void nlSparseMatrixMult(NLSparseMatrix* A, const NLdouble* x, NLdouble* y) {
    if(A->storage & NL_MATRIX_STORE_COMPRESSED) {
        nlCRSMatrixMult(A->compressed,x,y);
    } else if(A->storage & NL_MATRIX_STORE_ROWS) {
        if(A->storage & NL_MATRIX_STORE_SYMMETRIC) {
            nlSparseMatrix_mult_rows_symmetric(A, x, y) ;
        } else {
            nlSparseMatrix_mult_rows(A, x, y) ;
        }
    } else {
        if(A->storage & NL_MATRIX_STORE_SYMMETRIC) {
            nlSparseMatrix_mult_cols_symmetric(A, x, y) ;
        } else {
            nlSparseMatrix_mult_cols(A, x, y) ;
        }
    }
}

