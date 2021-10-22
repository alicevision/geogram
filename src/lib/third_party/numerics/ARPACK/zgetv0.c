/* ../FORTRAN/ARPACK/SRC/zgetv0.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    integer logfil, ndigit, mgetv0, msaupd, msaup2, msaitr, mseigt, msapps, 
	    msgets, mseupd, mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, 
	    mneupd, mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd;
} debug_;

#define debug_1 debug_

struct {
    integer nopx, nbx, nrorth, nitref, nrstrt;
    real tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv, tnaupd, 
	    tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv, tcaupd, tcaup2, 
	    tcaitr, tceigh, tcgets, tcapps, tcconv, tmvopx, tmvbx, tgetv0, 
	    titref, trvec;
} timing_;

#define timing_1 timing_

/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;

/* \BeginDoc */

/* \Name: zgetv0 */

/* \Description: */
/*  Generate a random initial residual vector for the Arnoldi process. */
/*  Force the residual vector to be in the range of the operator OP. */

/* \Usage: */
/*  call zgetv0 */
/*     ( IDO, BMAT, ITRY, INITV, N, J, V, LDV, RESID, RNORM, */
/*       IPNTR, WORKD, IERR ) */

/* \Arguments */
/*  IDO     Integer.  (INPUT/OUTPUT) */
/*          Reverse communication flag.  IDO must be zero on the first */
/*          call to zgetv0. */
/*          ------------------------------------------------------------- */
/*          IDO =  0: first call to the reverse communication interface */
/*          IDO = -1: compute  Y = OP * X  where */
/*                    IPNTR(1) is the pointer into WORKD for X, */
/*                    IPNTR(2) is the pointer into WORKD for Y. */
/*                    This is for the initialization phase to force the */
/*                    starting vector into the range of OP. */
/*          IDO =  2: compute  Y = B * X  where */
/*                    IPNTR(1) is the pointer into WORKD for X, */
/*                    IPNTR(2) is the pointer into WORKD for Y. */
/*          IDO = 99: done */
/*          ------------------------------------------------------------- */

/*  BMAT    Character*1.  (INPUT) */
/*          BMAT specifies the type of the matrix B in the (generalized) */
/*          eigenvalue problem A*x = lambda*B*x. */
/*          B = 'I' -> standard eigenvalue problem A*x = lambda*x */
/*          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x */

/*  ITRY    Integer.  (INPUT) */
/*          ITRY counts the number of times that zgetv0 is called. */
/*          It should be set to 1 on the initial call to zgetv0. */

/*  INITV   Logical variable.  (INPUT) */
/*          .TRUE.  => the initial residual vector is given in RESID. */
/*          .FALSE. => generate a random initial residual vector. */

/*  N       Integer.  (INPUT) */
/*          Dimension of the problem. */

/*  J       Integer.  (INPUT) */
/*          Index of the residual vector to be generated, with respect to */
/*          the Arnoldi process.  J > 1 in case of a "restart". */

/*  V       Complex*16 N by J array.  (INPUT) */
/*          The first J-1 columns of V contain the current Arnoldi basis */
/*          if this is a "restart". */

/*  LDV     Integer.  (INPUT) */
/*          Leading dimension of V exactly as declared in the calling */
/*          program. */

/*  RESID   Complex*16 array of length N.  (INPUT/OUTPUT) */
/*          Initial residual vector to be generated.  If RESID is */
/*          provided, force RESID into the range of the operator OP. */

/*  RNORM   Double precision scalar.  (OUTPUT) */
/*          B-norm of the generated residual. */

/*  IPNTR   Integer array of length 3.  (OUTPUT) */

/*  WORKD   Complex*16 work array of length 2*N.  (REVERSE COMMUNICATION). */
/*          On exit, WORK(1:N) = B*RESID to be used in SSAITR. */

/*  IERR    Integer.  (OUTPUT) */
/*          =  0: Normal exit. */
/*          = -1: Cannot generate a nontrivial restarted residual vector */
/*                in the range of the operator OP. */

/* \EndDoc */

/* ----------------------------------------------------------------------- */

/* \BeginLib */

/* \Local variables: */
/*     xxxxxx  Complex*16 */

/* \References: */
/*  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in */
/*     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992), */
/*     pp 357-385. */

/* \Routines called: */
/*     second  ARPACK utility routine for timing. */
/*     zvout   ARPACK utility routine that prints vectors. */
/*     zlarnv  LAPACK routine for generating a random vector. */
/*     zgemv   Level 2 BLAS routine for matrix vector multiplication. */
/*     zcopy   Level 1 BLAS that copies one vector to another. */
/*     zdotc   Level 1 BLAS that computes the scalar product of two vectors. */
/*     dznrm2  Level 1 BLAS that computes the norm of a vector. */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \SCCS Information: @(#) */
/* FILE: getv0.F   SID: 2.3   DATE OF SID: 08/27/96   RELEASE: 2 */

/* \EndLib */

/* ----------------------------------------------------------------------- */

/* Subroutine */ int zgetv0_(integer *ido, char *bmat, integer *itry, logical 
	*initv, integer *n, integer *j, doublecomplex *v, integer *ldv, 
	doublecomplex *resid, doublereal *rnorm, integer *ipntr, 
	doublecomplex *workd, integer *ierr, ftnlen bmat_len)
{
    /* Initialized data */

    static logical inits = TRUE_;

    /* System generated locals */
    integer v_dim1, v_offset, i__1, i__2;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static real t0, t1, t2, t3;
    static integer jj, iter;
    static logical orth;
    static integer iseed[4], idist;
    static doublecomplex cnorm;
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical first;
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    dvout_(integer *, integer *, doublereal *, integer *, char *, 
	    ftnlen), zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zvout_(integer *, integer *, 
	    doublecomplex *, integer *, char *, ftnlen);
    extern doublereal dlapy2_(doublereal *, doublereal *), dznrm2_(integer *, 
	    doublecomplex *, integer *);
    static doublereal rnorm0;
    extern /* Subroutine */ int second_(real *);
    static integer msglvl;
    extern /* Subroutine */ int zlarnv_(integer *, integer *, integer *, 
	    doublecomplex *);


/*     %----------------------------------------------------% */
/*     | Include files for debugging and timing information | */
/*     %----------------------------------------------------% */


/* \SCCS Information: @(#) */
/* FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 */

/*     %---------------------------------% */
/*     | See debug.doc for documentation | */
/*     %---------------------------------% */

/*     %------------------% */
/*     | Scalar Arguments | */
/*     %------------------% */

/*     %--------------------------------% */
/*     | See stat.doc for documentation | */
/*     %--------------------------------% */

/* \SCCS Information: @(#) */
/* FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 */



/*     %-----------------% */
/*     | Array Arguments | */
/*     %-----------------% */


/*     %------------% */
/*     | Parameters | */
/*     %------------% */


/*     %------------------------% */
/*     | Local Scalars & Arrays | */
/*     %------------------------% */


/*     %----------------------% */
/*     | External Subroutines | */
/*     %----------------------% */


/*     %--------------------% */
/*     | External Functions | */
/*     %--------------------% */


/*     %-----------------% */
/*     | Data Statements | */
/*     %-----------------% */

    /* Parameter adjustments */
    --workd;
    --resid;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --ipntr;

    /* Function Body */

/*     %-----------------------% */
/*     | Executable Statements | */
/*     %-----------------------% */


/*     %-----------------------------------% */
/*     | Initialize the seed of the LAPACK | */
/*     | random number generator           | */
/*     %-----------------------------------% */

    if (inits) {
	iseed[0] = 1;
	iseed[1] = 3;
	iseed[2] = 5;
	iseed[3] = 7;
	inits = FALSE_;
    }

    if (*ido == 0) {

/*        %-------------------------------% */
/*        | Initialize timing statistics  | */
/*        | & message level for debugging | */
/*        %-------------------------------% */

	second_(&t0);
	msglvl = debug_1.mgetv0;

	*ierr = 0;
	iter = 0;
	first = FALSE_;
	orth = FALSE_;

/*        %-----------------------------------------------------% */
/*        | Possibly generate a random starting vector in RESID | */
/*        | Use a LAPACK random number generator used by the    | */
/*        | matrix generation routines.                         | */
/*        |    idist = 1: uniform (0,1)  distribution;          | */
/*        |    idist = 2: uniform (-1,1) distribution;          | */
/*        |    idist = 3: normal  (0,1)  distribution;          | */
/*        %-----------------------------------------------------% */

	if (! (*initv)) {
	    idist = 2;
	    zlarnv_(&idist, iseed, n, &resid[1]);
	}

/*        %----------------------------------------------------------% */
/*        | Force the starting vector into the range of OP to handle | */
/*        | the generalized problem when B is possibly (singular).   | */
/*        %----------------------------------------------------------% */

	second_(&t2);
	if (*(unsigned char *)bmat == 'G') {
	    ++timing_1.nopx;
	    ipntr[1] = 1;
	    ipntr[2] = *n + 1;
	    zcopy_(n, &resid[1], &c__1, &workd[1], &c__1);
	    *ido = -1;
	    goto L9000;
	}
    }

/*     %----------------------------------------% */
/*     | Back from computing B*(initial-vector) | */
/*     %----------------------------------------% */

    if (first) {
	goto L20;
    }

/*     %-----------------------------------------------% */
/*     | Back from computing B*(orthogonalized-vector) | */
/*     %-----------------------------------------------% */

    if (orth) {
	goto L40;
    }

    second_(&t3);
    timing_1.tmvopx += t3 - t2;

/*     %------------------------------------------------------% */
/*     | Starting vector is now in the range of OP; r = OP*r; | */
/*     | Compute B-norm of starting vector.                   | */
/*     %------------------------------------------------------% */

    second_(&t2);
    first = TRUE_;
    if (*(unsigned char *)bmat == 'G') {
	++timing_1.nbx;
	zcopy_(n, &workd[*n + 1], &c__1, &resid[1], &c__1);
	ipntr[1] = *n + 1;
	ipntr[2] = 1;
	*ido = 2;
	goto L9000;
    } else if (*(unsigned char *)bmat == 'I') {
	zcopy_(n, &resid[1], &c__1, &workd[1], &c__1);
    }

L20:

    if (*(unsigned char *)bmat == 'G') {
	second_(&t3);
	timing_1.tmvbx += t3 - t2;
    }

    first = FALSE_;
    if (*(unsigned char *)bmat == 'G') {
	zdotc_(&z__1, n, &resid[1], &c__1, &workd[1], &c__1);
	cnorm.r = z__1.r, cnorm.i = z__1.i;
	d__1 = cnorm.r;
	d__2 = d_imag(&cnorm);
	rnorm0 = sqrt(dlapy2_(&d__1, &d__2));
    } else if (*(unsigned char *)bmat == 'I') {
	rnorm0 = dznrm2_(n, &resid[1], &c__1);
    }
    *rnorm = rnorm0;

/*     %---------------------------------------------% */
/*     | Exit if this is the very first Arnoldi step | */
/*     %---------------------------------------------% */

    if (*j == 1) {
	goto L50;
    }

/*     %---------------------------------------------------------------- */
/*     | Otherwise need to B-orthogonalize the starting vector against | */
/*     | the current Arnoldi basis using Gram-Schmidt with iter. ref.  | */
/*     | This is the case where an invariant subspace is encountered   | */
/*     | in the middle of the Arnoldi factorization.                   | */
/*     |                                                               | */
/*     |       s = V^{T}*B*r;   r = r - V*s;                           | */
/*     |                                                               | */
/*     | Stopping criteria used for iter. ref. is discussed in         | */
/*     | Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   | */
/*     %---------------------------------------------------------------% */

    orth = TRUE_;
L30:

    i__1 = *j - 1;
    zgemv_("C", n, &i__1, &c_b1, &v[v_offset], ldv, &workd[1], &c__1, &c_b2, &
	    workd[*n + 1], &c__1, (ftnlen)1);
    i__1 = *j - 1;
    z__1.r = -1., z__1.i = -0.;
    zgemv_("N", n, &i__1, &z__1, &v[v_offset], ldv, &workd[*n + 1], &c__1, &
	    c_b1, &resid[1], &c__1, (ftnlen)1);

/*     %----------------------------------------------------------% */
/*     | Compute the B-norm of the orthogonalized starting vector | */
/*     %----------------------------------------------------------% */

    second_(&t2);
    if (*(unsigned char *)bmat == 'G') {
	++timing_1.nbx;
	zcopy_(n, &resid[1], &c__1, &workd[*n + 1], &c__1);
	ipntr[1] = *n + 1;
	ipntr[2] = 1;
	*ido = 2;
	goto L9000;
    } else if (*(unsigned char *)bmat == 'I') {
	zcopy_(n, &resid[1], &c__1, &workd[1], &c__1);
    }

L40:

    if (*(unsigned char *)bmat == 'G') {
	second_(&t3);
	timing_1.tmvbx += t3 - t2;
    }

    if (*(unsigned char *)bmat == 'G') {
	zdotc_(&z__1, n, &resid[1], &c__1, &workd[1], &c__1);
	cnorm.r = z__1.r, cnorm.i = z__1.i;
	d__1 = cnorm.r;
	d__2 = d_imag(&cnorm);
	*rnorm = sqrt(dlapy2_(&d__1, &d__2));
    } else if (*(unsigned char *)bmat == 'I') {
	*rnorm = dznrm2_(n, &resid[1], &c__1);
    }

/*     %--------------------------------------% */
/*     | Check for further orthogonalization. | */
/*     %--------------------------------------% */

    if (msglvl > 2) {
	dvout_(&debug_1.logfil, &c__1, &rnorm0, &debug_1.ndigit, "_getv0: re"
		"-orthonalization ; rnorm0 is", (ftnlen)38);
	dvout_(&debug_1.logfil, &c__1, rnorm, &debug_1.ndigit, "_getv0: re-o"
		"rthonalization ; rnorm is", (ftnlen)37);
    }

    if (*rnorm > rnorm0 * .717f) {
	goto L50;
    }

    ++iter;
    if (iter <= 1) {

/*        %-----------------------------------% */
/*        | Perform iterative refinement step | */
/*        %-----------------------------------% */

	rnorm0 = *rnorm;
	goto L30;
    } else {

/*        %------------------------------------% */
/*        | Iterative refinement step "failed" | */
/*        %------------------------------------% */

	i__1 = *n;
	for (jj = 1; jj <= i__1; ++jj) {
	    i__2 = jj;
	    resid[i__2].r = 0., resid[i__2].i = 0.;
/* L45: */
	}
	*rnorm = 0.;
	*ierr = -1;
    }

L50:

    if (msglvl > 0) {
	dvout_(&debug_1.logfil, &c__1, rnorm, &debug_1.ndigit, "_getv0: B-no"
		"rm of initial / restarted starting vector", (ftnlen)53);
    }
    if (msglvl > 2) {
	zvout_(&debug_1.logfil, n, &resid[1], &debug_1.ndigit, "_getv0: init"
		"ial / restarted starting vector", (ftnlen)43);
    }
    *ido = 99;

    second_(&t1);
    timing_1.tgetv0 += t1 - t0;

L9000:
    return 0;

/*     %---------------% */
/*     | End of zgetv0 | */
/*     %---------------% */

} /* zgetv0_ */

