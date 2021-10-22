/* ../FORTRAN/ARPACK/SRC/ssconv.f -- translated by f2c (version 20100827).
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

static doublereal c_b3 = .66666666666666663;

/* ----------------------------------------------------------------------- */
/* \BeginDoc */

/* \Name: ssconv */

/* \Description: */
/*  Convergence testing for the symmetric Arnoldi eigenvalue routine. */

/* \Usage: */
/*  call ssconv */
/*     ( N, RITZ, BOUNDS, TOL, NCONV ) */

/* \Arguments */
/*  N       Integer.  (INPUT) */
/*          Number of Ritz values to check for convergence. */

/*  RITZ    Real array of length N.  (INPUT) */
/*          The Ritz values to be checked for convergence. */

/*  BOUNDS  Real array of length N.  (INPUT) */
/*          Ritz estimates associated with the Ritz values in RITZ. */

/*  TOL     Real scalar.  (INPUT) */
/*          Desired relative accuracy for a Ritz value to be considered */
/*          "converged". */

/*  NCONV   Integer scalar.  (OUTPUT) */
/*          Number of "converged" Ritz values. */

/* \EndDoc */

/* ----------------------------------------------------------------------- */

/* \BeginLib */

/* \Routines called: */
/*     second  ARPACK utility routine for timing. */
/*     slamch  LAPACK routine that determines machine constants. */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \SCCS Information: @(#) */
/* FILE: sconv.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2 */

/* \Remarks */
/*     1. Starting with version 2.4, this routine no longer uses the */
/*        Parlett strategy using the gap conditions. */

/* \EndLib */

/* ----------------------------------------------------------------------- */

/* Subroutine */ int ssconv_(integer *n, real *ritz, real *bounds, real *tol, 
	integer *nconv)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__;
    static real t0, t1, eps23, temp;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int second_(real *);


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


/*     %---------------% */
/*     | Local Scalars | */
/*     %---------------% */


/*     %-------------------% */
/*     | External routines | */
/*     %-------------------% */

/*     %---------------------% */
/*     | Intrinsic Functions | */
/*     %---------------------% */


/*     %-----------------------% */
/*     | Executable Statements | */
/*     %-----------------------% */

    /* Parameter adjustments */
    --bounds;
    --ritz;

    /* Function Body */
    second_(&t0);

    eps23 = slamch_("Epsilon-Machine", (ftnlen)15);
    d__1 = (doublereal) eps23;
    eps23 = pow_dd(&d__1, &c_b3);

    *nconv = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        %-----------------------------------------------------% */
/*        | The i-th Ritz value is considered "converged"       | */
/*        | when: bounds(i) .le. TOL*max(eps23, abs(ritz(i)))   | */
/*        %-----------------------------------------------------% */

/* Computing MAX */
	r__2 = eps23, r__3 = (r__1 = ritz[i__], dabs(r__1));
	temp = dmax(r__2,r__3);
	if (bounds[i__] <= *tol * temp) {
	    ++(*nconv);
	}

/* L10: */
    }

    second_(&t1);
    timing_1.tsconv += t1 - t0;

    return 0;

/*     %---------------% */
/*     | End of ssconv | */
/*     %---------------% */

} /* ssconv_ */

