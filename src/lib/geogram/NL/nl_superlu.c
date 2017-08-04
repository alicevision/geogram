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

#include "nl_superlu.h"
#include "nl_context.h"

/**
 * \file nl_superlu.c
 * \brief Weak-coupling adapter to call SuperLU from OpenNL, 
 *  works with both SuperLU 3.x and SuperLU 4.x.
 */

#ifdef NL_OS_UNIX
#  ifdef NL_OS_APPLE
#      define SUPERLU_LIB_NAME "libsuperlu_4.dylib"
#  else
#      define SUPERLU_LIB_NAME "libsuperlu.so"
#  endif
#else
#  define SUPERLU_LIB_NAME "libsuperlu.xxx"
#endif


/**********************************/
/** Exerpt from SuperLU includes **/
/** <slu_cdefs.h>                **/
/** <supermatrix.h>              **/
/**                              **/
/**********************************/

/*
 * Important note: 
 * the order of some constants and the size
 * of some structures have changed between
 * SuperLU 3.x and SuperLU 4.x.
 * See documentation of the SuperLU_version() 
 * function in nl_superlu.c for more details.
 */

typedef enum {
    /* SuperLU 3.x */
    SLU3_NC  =0,   /* column-wise, no supernode */
    SLU3_NR  =1,   /* row-wize, no supernode */
    SLU3_SC  =2,   /* column-wise, supernode */
    SLU3_SR  =3,   /* row-wise, supernode */
    SLU3_NCP =4,   /* column-wise, column-permuted, no supernode*/ 
    SLU3_DN  =5,   /* Fortran style column-wise storage for dense matrix */

    /* SuperLU 4.x */    
    SLU4_NC     =0, /* column-wise, no supernode */
    SLU4_NCP    =1, /* column-wise, column-permuted, no supernode */ 
    SLU4_NR     =2, /* row-wize, no supernode */
    SLU4_SC     =3, /* column-wise, supernode */
    SLU4_SCP    =4, /* supernode, column-wise, permuted */
    SLU4_SR     =5, /* row-wise, supernode */
    SLU4_DN     =6, /* Fortran style column-wise storage for dense matrix */
    SLU4_NR_loc =7  /* distributed compressed row format  */
} Stype_t;

typedef enum {
    SLU_S,     /* single */
    SLU_D,     /* double */
    SLU_C,     /* single complex */
    SLU_Z      /* double complex */
} Dtype_t;

typedef enum {
    SLU_GE,    /* general */
    SLU_TRLU,  /* lower triangular, unit diagonal */
    SLU_TRUU,  /* upper triangular, unit diagonal */
    SLU_TRL,   /* lower triangular */
    SLU_TRU,   /* upper triangular */
    SLU_SYL,   /* symmetric, store lower half */
    SLU_SYU,   /* symmetric, store upper half */
    SLU_HEL,   /* Hermitian, store lower half */
    SLU_HEU    /* Hermitian, store upper half */
} Mtype_t;

typedef int int_t;

typedef struct {
        Stype_t Stype; /* Storage type: interprets the storage structure 
                          pointed to by *Store. */
        Dtype_t Dtype; /* Data type. */
        Mtype_t Mtype; /* Matrix type: describes the mathematical property of 
                          the matrix. */
        int_t  nrow;   /* number of rows */
        int_t  ncol;   /* number of columns */
        void *Store;   /* pointer to the actual storage of the matrix */
} SuperMatrix;

typedef struct {
    int_t lda;    /* leading dimension */
    void *nzval;  /* array of size lda*ncol to represent a dense matrix */
} DNformat;

typedef enum {NO, YES}                                          yes_no_t;
typedef enum {DOFACT, SamePattern, SamePattern_SameRowPerm, FACTORED} fact_t;
typedef enum {NOROWPERM, LargeDiag, MY_PERMR}                   rowperm_t;

typedef enum {
    /* SuperLU 3.x */
    SLU3_NATURAL       = 0,
    SLU3_MMD_ATA       = 1,
    SLU3_MMD_AT_PLUS_A = 2,
    SLU3_COLAMD        = 3,
    SLU3_MY_PERMC      = 4,

    /* SuperLU 4.x */
    SLU4_NATURAL         = 0,
    SLU4_MMD_ATA         = 1,
    SLU4_MMD_AT_PLUS_A   = 2,
    SLU4_COLAMD          = 3,
    SLU4_METIS_AT_PLUS_A = 4,
    SLU4_PARMETIS        = 5,
    SLU4_ZOLTAN          = 6,
    SLU4_MY_PERMC        = 7
} colperm_t;


typedef enum {NOTRANS, TRANS, CONJ}                             trans_t;
typedef enum {NOEQUIL, ROW, COL, BOTH}                          DiagScale_t;
typedef enum {NOREFINE, SLU_SINGLE=1, SLU_DOUBLE, SLU_EXTRA}    IterRefine_t;
typedef enum {LUSUP, UCOL, LSUB, USUB, SLU4_LLVL, SLU4_ULVL}    MemType;
typedef enum {HEAD, TAIL}                                       stack_end_t;
typedef enum {SYSTEM, USER}                                     LU_space_t;
typedef enum {SLU4_ONE_NORM, SLU4_TWO_NORM, SLU4_INF_NORM}      norm_t;
typedef enum {
    SLU4_SILU, SLU4_SMILU_1, SLU4_SMILU_2, SLU4_SMILU_3
} milu_t;

typedef struct {
    fact_t        Fact;
    yes_no_t      Equil;
    colperm_t     ColPerm;
    trans_t       Trans;
    IterRefine_t  IterRefine;
    yes_no_t      PrintStat;
    yes_no_t      SymmetricMode;
    double        DiagPivotThresh;
    yes_no_t      PivotGrowth;
    yes_no_t      ConditionNumber;
    rowperm_t     RowPerm;
    yes_no_t      ReplaceTinyPivot;
    yes_no_t      SolveInitialized;
    yes_no_t      RefineInitialized;
} superlu3_options_t;

typedef struct {
    fact_t        Fact;
    yes_no_t      Equil;
    colperm_t     ColPerm;
    trans_t       Trans;
    IterRefine_t  IterRefine;
    double        DiagPivotThresh;
    yes_no_t      SymmetricMode;
    yes_no_t      PivotGrowth;
    yes_no_t      ConditionNumber;
    rowperm_t     RowPerm;
    int           ILU_DropRule;
    double        ILU_DropTol;    /* threshold for dropping */
    double        ILU_FillFactor; /* gamma in the secondary dropping */
    norm_t        ILU_Norm;       /* infinity-norm, 1-norm, or 2-norm */
    double        ILU_FillTol;    /* threshold for zero pivot perturbation */
    milu_t        ILU_MILU;
    double        ILU_MILU_Dim;   /* Dimension of PDE (if available) */
    yes_no_t      ParSymbFact;
    yes_no_t      ReplaceTinyPivot; /* used in SuperLU_DIST */
    yes_no_t      SolveInitialized;
    yes_no_t      RefineInitialized;
    yes_no_t      PrintStat;
    int           nnzL, nnzU;      /* used to store nnzs for now       */
    int           num_lookaheads;  /* num of levels in look-ahead      */
    yes_no_t      lookahead_etree; /* use etree computed from the
                                      serial symbolic factorization */
    yes_no_t      SymPattern;      /* symmetric factorization          */
} superlu4_options_t;

typedef void* superlu_options_ptr;

typedef float    flops_t;
typedef unsigned char Logical;

typedef struct {
    int     *panel_histo;    /* histogram of panel size distribution */
    double  *utime;          /* running time at various phases */
    flops_t *ops;            /* operation count at various phases */
    int     TinyPivots;      /* number of tiny pivots */
    int     RefineSteps;     /* number of iterative refinement steps */
    int     slu4_expansions; /* number of memory expansions (SuperLU4) */
} SuperLUStat_t;

/*****************************************/
/** End of exerpt from SuperLU includes **/
/**                                     **/
/*****************************************/

/*****************************************/
/** Functions pointers to allow dynamic **/
/** linking of superlu libs.            **/
/**                                     **/
/*****************************************/

typedef void (*FUNPTR_set_default_options)(superlu_options_ptr options);

typedef void (*FUNPTR_StatInit)(SuperLUStat_t *);
typedef void (*FUNPTR_StatFree)(SuperLUStat_t *);

typedef void (*FUNPTR_dCreate_CompCol_Matrix)(
    SuperMatrix *, int, int, int, const double *,
    const int *, const int *, Stype_t, Dtype_t, Mtype_t);

typedef void (*FUNPTR_dCreate_Dense_Matrix)(
    SuperMatrix *, int, int, const double *, int,
    Stype_t, Dtype_t, Mtype_t);

typedef void (*FUNPTR_Destroy_SuperNode_Matrix)(SuperMatrix *);
typedef void (*FUNPTR_Destroy_CompCol_Matrix)(SuperMatrix *);
typedef void (*FUNPTR_Destroy_SuperMatrix_Store)(SuperMatrix *);

typedef void (*FUNPTR_dgssv)(
    superlu_options_ptr, SuperMatrix *, int *, int *, SuperMatrix *,
    SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int *
);

/**
 * \brief The structure that stores the handle to 
 *  the SuperLU shared object, the function pointers
 *  and the detected version.
 */
typedef struct {
    FUNPTR_set_default_options set_default_options;
    FUNPTR_StatInit StatInit;
    FUNPTR_StatFree StatFree;
    FUNPTR_dCreate_CompCol_Matrix dCreate_CompCol_Matrix;
    FUNPTR_dCreate_Dense_Matrix dCreate_Dense_Matrix;
    FUNPTR_Destroy_SuperNode_Matrix Destroy_SuperNode_Matrix;
    FUNPTR_Destroy_CompCol_Matrix Destroy_CompCol_Matrix;
    FUNPTR_Destroy_SuperMatrix_Store Destroy_SuperMatrix_Store;
    FUNPTR_dgssv dgssv;

    NLdll DLL_handle;

    double version;
} SuperLUContext;

/**
 * \brief Gets the SuperLU context.
 * \return a pointer to the SuperLU context
 */
static SuperLUContext* SuperLU() {
    static SuperLUContext context;
    static NLboolean init = NL_FALSE;
    if(!init) {
        init = NL_TRUE;
        memset(&context, 0, sizeof(context));
    }
    return &context;
}

/*
 * \brief Gets the version of SuperLU that was dynamically loaded.
 * \details It is important to know the version of SuperLU because
 *  the order of some constants and the size of some structures 
 *  have changed between SuperLU 3.x and SuperLU 4.x.
 *  When there is a mismatch between both versions, constants are 
 *  prefixed by SLU3_ or SLU4_ according to the version.
 * This concerns:
 * - enum constants in Stype_t
 * - enum constants in colperm_t
 * - enum constants in norm_t and milu_t (that only exist in version 4.x)
 * - struct superlu_options_t (use instead superlu3_options_t 
 *  or superlu4_options_t according to version)
 */
static double SuperLU_version() {
    return SuperLU()->version;
}

/**
 * \brief Tests whether SuperLU extension is
 *  initialized.
 * \details Tests whether SuperLU shared object
 *  was successfuly loaded and whether all the
 *  function pointers where found.
 * \retval NL_TRUE if SuperLU was successfully
 *  loaded and initialized
 * \retval NL_FALSE otherwise
 */
static NLboolean SuperLU_is_initialized() {
    return
        SuperLU()->DLL_handle != NULL &&
        SuperLU()->set_default_options != NULL &&
        SuperLU()->StatInit != NULL &&
        SuperLU()->StatFree != NULL &&
        SuperLU()->dCreate_CompCol_Matrix != NULL &&
        SuperLU()->dCreate_Dense_Matrix != NULL &&
        SuperLU()->Destroy_SuperNode_Matrix != NULL &&
        SuperLU()->Destroy_CompCol_Matrix != NULL &&
        SuperLU()->Destroy_SuperMatrix_Store != NULL &&
        SuperLU()->dgssv != NULL;
}

/**
 * \brief Finds and initializes a function pointer to
 *  one of the functions in SuperLU.
 * \details Function pointers are stored into the 
 *  SuperLUContext returned by the function SuperLU().
 *  If a symbol is not found, returns NL_FALSE from the
 *  calling function.
 */
#define find_superlu_func(name)                                   \
    if(                                                           \
        (                                                         \
            SuperLU()->name =                                     \
            (FUNPTR_##name)nlFindFunction(SuperLU()->DLL_handle,#name) \
        ) == NULL                                                 \
    ) {                                                           \
        nlError("nlInitExtension_SUPERLU","function not found");  \
        return NL_FALSE;                                          \
    }

/**
 * \brief Version-independent constant for SLU_NR 
 *  (row-wize, no supernode)
 * \see SuperLU_version()
 */
#define SLU_NR ((SuperLU()->version >= 4.0) ? SLU4_NR : SLU3_NR)

/**
 * \brief Version-independent constant for SLU_DN 
 *  (Fortran style column-wise storage for dense matrix)
 * \see SuperLU_version()
 */
#define SLU_DN ((SuperLU()->version >= 4.0) ? SLU4_DN : SLU3_DN)

/************************************************************************/

NLboolean nlSolve_system_with_SUPERLU(
    NLSparseMatrix* M, double* x, const double* b,
    NLenum solver, NLboolean clear_M
) {
    /* Compressed Row Storage matrix representation */
    NLuint    n      = M->n;
    NLuint    nnz    = nlSparseMatrixNNZ(M); /* Number of Non-Zero coeffs */
    NLint*    xa     = NL_NEW_ARRAY(NLint, n+1);
    NLdouble* rhs    = NL_NEW_ARRAY(NLdouble, n);
    NLdouble* a      = NL_NEW_ARRAY(NLdouble, nnz);
    NLint*    asub   = NL_NEW_ARRAY(NLint, nnz);

    /* Permutation vector */
    NLint*    perm_r  = NL_NEW_ARRAY(NLint, n);
    NLint*    perm    = NL_NEW_ARRAY(NLint, n);

    /* SuperLU variables */
    SuperMatrix A, B; /* System       */
    SuperMatrix L, U; /* Factorization of A */
    NLint info;       /* status code  */
    DNformat *vals = NULL; /* access to result */
    double *rvals  = NULL; /* access to result */

    /* SuperLU options and stats */
    superlu3_options_t options3;
    superlu4_options_t options4;
    SuperLUStat_t     stat;

    /* Temporary variables */
    NLRowColumn* Ri = NULL;
    NLuint         i,jj,count;
    
    /* Sanity checks */
    nl_assert(!(M->storage & NL_MATRIX_STORE_SYMMETRIC));
    nl_assert(M->storage & NL_MATRIX_STORE_ROWS);
    nl_assert(M->m == M->n);

    if(!SuperLU_is_initialized()) {
        nlError(
            "nlSolve_SUPERLU",
            "SuperLU extension not initialized (nlInitExtension(\"SUPERLU\") missing or failed)"
        );
        return NL_FALSE;
    }
    
    /*
     * Step 1: convert matrix M into SuperLU compressed column 
     *   representation.
     * -------------------------------------------------------
     */

    count = 0;
    for(i=0; i<n; i++) {
        Ri = &(M->row[i]);
        xa[i] = (NLint)(count);
        for(jj=0; jj<Ri->size; jj++) {
            a[count]    = Ri->coeff[jj].value;
            asub[count] = (NLint)(Ri->coeff[jj].index);
            count++;
        }
    }
    xa[n] = (NLint)(nnz);

    /* Save memory for SuperLU */
    if(clear_M) {
        nlSparseMatrixClear(M);
    }


    /*
     * Rem: SuperLU does not support symmetric storage.
     * In fact, for symmetric matrix, what we need 
     * is a SuperLLt algorithm (SuperNodal sparse Cholesky),
     * TODO: interface other solvers from suitesparse.
     * However, this is not a big problem (SuperLU is just
     * a superset of what we really need.
     */
    SuperLU()->dCreate_CompCol_Matrix(
        &A, (int)n, (int)n, (int)nnz,
        a, asub, xa, 
        SLU_NR,              /* Row_wise, no supernode */
        SLU_D,               /* doubles                */ 
        SLU_GE               /* general storage        */
    );

    /* Step 2: create vector */
    SuperLU()->dCreate_Dense_Matrix(
        &B, (int)n, 1, b, (int)n, 
        SLU_DN, /* Fortran-type column-wise storage */
        SLU_D,  /* doubles                          */
        SLU_GE  /* general                          */
    );
            

    /* Step 3: set SuperLU options 
     * ------------------------------
     */

    if(SuperLU_version() >= 4.0) {
        SuperLU()->set_default_options(&options4);
        switch(solver) {
        case NL_SUPERLU_EXT: {
            options4.ColPerm = SLU4_NATURAL;
        } break;
        case NL_PERM_SUPERLU_EXT: {
            options4.ColPerm = SLU4_COLAMD;
        } break;
        case NL_SYMMETRIC_SUPERLU_EXT: {
            options4.ColPerm = SLU4_MMD_AT_PLUS_A;
            options4.SymmetricMode = YES;
        } break;
        default: 
            nl_assert_not_reached;
        }
    } else {
        SuperLU()->set_default_options(&options3);
        switch(solver) {
        case NL_SUPERLU_EXT: {
            options3.ColPerm = SLU3_NATURAL;
        } break;
        case NL_PERM_SUPERLU_EXT: {
            options3.ColPerm = SLU3_COLAMD;
        } break;
        case NL_SYMMETRIC_SUPERLU_EXT: {
            options3.ColPerm = SLU3_MMD_AT_PLUS_A;
            options3.SymmetricMode = YES;
        } break;
        default: 
            nl_assert_not_reached;
        }
    }
    
    SuperLU()->StatInit(&stat);

    /* Step 4: call SuperLU main routine
     * ---------------------------------
     */
    
    if(SuperLU_version() >= 4.0) {
        SuperLU()->dgssv(
            &options4, &A, perm, perm_r, &L, &U, &B, &stat, &info
        );
    } else {
        SuperLU()->dgssv(
            &options3, &A, perm, perm_r, &L, &U, &B, &stat, &info
        );
    }

    /* Step 5: get the solution
     * ------------------------
     * Fortran-type column-wise storage
     */
    vals = (DNformat*)B.Store;
    rvals = (double*)(vals->nzval);
    if(info == 0) {
        for(i = 0; i <  n; i++){
            x[i] = rvals[i];
        }
    } else {
        nlError("nlSolve", "SuperLU failed");
    }

    /* Step 6: cleanup
     * ---------------
     */

    /*
     *  For these two ones, only the "store" structure
     * needs to be deallocated (the arrays have been allocated
     * by us).
     */
    SuperLU()->Destroy_SuperMatrix_Store(&A);
    SuperLU()->Destroy_SuperMatrix_Store(&B);

    /*
     *   These ones need to be fully deallocated (they have been
     * allocated by SuperLU).
     */
    SuperLU()->Destroy_SuperNode_Matrix(&L);
    SuperLU()->Destroy_CompCol_Matrix(&U);

    /* There are some dynamically allocated vectors in the stats */
    SuperLU()->StatFree(&stat);

    NL_DELETE_ARRAY(xa);
    NL_DELETE_ARRAY(rhs);
    NL_DELETE_ARRAY(a);
    NL_DELETE_ARRAY(asub);
    NL_DELETE_ARRAY(perm_r);
    NL_DELETE_ARRAY(perm);

    return (info == 0);
}

/************************************************************************/

NLboolean nlSolve_SUPERLU(void) {
    /* Get current linear system from context */
    NLSparseMatrix* M  = &(nlCurrentContext->M);
    NLdouble* b          = nlCurrentContext->b;
    NLdouble* x          = nlCurrentContext->x;

    /* Sanity checks */
    nl_assert(!(M->storage & NL_MATRIX_STORE_SYMMETRIC));
    nl_assert(M->storage & NL_MATRIX_STORE_ROWS);
    nl_assert(M->m == M->n);
    
    return nlSolve_system_with_SUPERLU(
        M, x, b, nlCurrentContext->solver, NL_TRUE
    );
}


static void nlTerminateExtension_SUPERLU(void) {
    if(SuperLU()->DLL_handle != NULL) {
        nlCloseDLL(SuperLU()->DLL_handle);
        SuperLU()->DLL_handle = NULL;
    }
}

NLboolean nlInitExtension_SUPERLU(void) {
    
    if(SuperLU()->DLL_handle != NULL) {
        return SuperLU_is_initialized();
    }

    SuperLU()->DLL_handle = nlOpenDLL(SUPERLU_LIB_NAME);
    if(SuperLU()->DLL_handle == NULL) {
        return NL_FALSE;
    }

    /* 
     * Check for SuperLU version:
     * Since ILU (incomplete Cholesky) is only available in 4.x, if
     * we find one of the ILU-related symbols in there, then we got
     * a 4.x.
     * TODO: there may be a finer way to detect version.
     */
    if(nlFindFunction(SuperLU()->DLL_handle,"ilu_set_default_options") != NULL) {
        SuperLU()->version = 4.0;
    } else {
        SuperLU()->version = 3.0;
    }
    
    find_superlu_func(set_default_options);
    find_superlu_func(StatInit);
    find_superlu_func(StatFree);
    find_superlu_func(dCreate_CompCol_Matrix);
    find_superlu_func(dCreate_Dense_Matrix);
    find_superlu_func(Destroy_SuperNode_Matrix);
    find_superlu_func(Destroy_CompCol_Matrix);    
    find_superlu_func(Destroy_SuperMatrix_Store);
    find_superlu_func(dgssv);

    atexit(nlTerminateExtension_SUPERLU);
    return NL_TRUE;
}



