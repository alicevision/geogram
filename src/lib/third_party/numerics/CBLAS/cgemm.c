/* cgemm.f -- translated by f2c (version 20061008).
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
#include "blaswrap.h"

/* Subroutine */ int cgemm_(char *transa, char *transb, integer *m, integer *
	n, integer *k, complex *alpha, complex *a, integer *lda, complex *b, 
	integer *ldb, complex *beta, complex *c__, integer *ldc)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4, i__5, i__6;
    complex q__1, q__2, q__3, q__4;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);

    /* Local variables */
    integer i__, j, l, info;
    logical nota, notb;
    complex temp;
    logical conja, conjb;
    integer ncola;
    extern logical lsame_(char *, char *);
    integer nrowa, nrowb;
    extern /* Subroutine */ int xerbla_(char *, integer *);

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  CGEMM  performs one of the matrix-matrix operations */

/*     C := alpha*op( A )*op( B ) + beta*C, */

/*  where  op( X ) is one of */

/*     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ), */

/*  alpha and beta are scalars, and A, B and C are matrices, with op( A ) */
/*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix. */

/*  Arguments */
/*  ========== */

/*  TRANSA - CHARACTER*1. */
/*           On entry, TRANSA specifies the form of op( A ) to be used in */
/*           the matrix multiplication as follows: */

/*              TRANSA = 'N' or 'n',  op( A ) = A. */

/*              TRANSA = 'T' or 't',  op( A ) = A'. */

/*              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ). */

/*           Unchanged on exit. */

/*  TRANSB - CHARACTER*1. */
/*           On entry, TRANSB specifies the form of op( B ) to be used in */
/*           the matrix multiplication as follows: */

/*              TRANSB = 'N' or 'n',  op( B ) = B. */

/*              TRANSB = 'T' or 't',  op( B ) = B'. */

/*              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ). */

/*           Unchanged on exit. */

/*  M      - INTEGER. */
/*           On entry,  M  specifies  the number  of rows  of the  matrix */
/*           op( A )  and of the  matrix  C.  M  must  be at least  zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry,  N  specifies the number  of columns of the matrix */
/*           op( B ) and the number of columns of the matrix C. N must be */
/*           at least zero. */
/*           Unchanged on exit. */

/*  K      - INTEGER. */
/*           On entry,  K  specifies  the number of columns of the matrix */
/*           op( A ) and the number of rows of the matrix op( B ). K must */
/*           be at least  zero. */
/*           Unchanged on exit. */

/*  ALPHA  - COMPLEX         . */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is */
/*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise. */
/*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k */
/*           part of the array  A  must contain the matrix  A,  otherwise */
/*           the leading  k by m  part of the array  A  must contain  the */
/*           matrix A. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then */
/*           LDA must be at least  max( 1, m ), otherwise  LDA must be at */
/*           least  max( 1, k ). */
/*           Unchanged on exit. */

/*  B      - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is */
/*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise. */
/*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n */
/*           part of the array  B  must contain the matrix  B,  otherwise */
/*           the leading  n by k  part of the array  B  must contain  the */
/*           matrix B. */
/*           Unchanged on exit. */

/*  LDB    - INTEGER. */
/*           On entry, LDB specifies the first dimension of B as declared */
/*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then */
/*           LDB must be at least  max( 1, k ), otherwise  LDB must be at */
/*           least  max( 1, n ). */
/*           Unchanged on exit. */

/*  BETA   - COMPLEX         . */
/*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is */
/*           supplied as zero then C need not be set on input. */
/*           Unchanged on exit. */

/*  C      - COMPLEX          array of DIMENSION ( LDC, n ). */
/*           Before entry, the leading  m by n  part of the array  C must */
/*           contain the matrix  C,  except when  beta  is zero, in which */
/*           case C need not be set on entry. */
/*           On exit, the array  C  is overwritten by the  m by n  matrix */
/*           ( alpha*op( A )*op( B ) + beta*C ). */

/*  LDC    - INTEGER. */
/*           On entry, LDC specifies the first dimension of C as declared */
/*           in  the  calling  (sub)  program.   LDC  must  be  at  least */
/*           max( 1, m ). */
/*           Unchanged on exit. */


/*  Level 3 Blas routine. */

/*  -- Written on 8-February-1989. */
/*     Jack Dongarra, Argonne National Laboratory. */
/*     Iain Duff, AERE Harwell. */
/*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/*     Sven Hammarling, Numerical Algorithms Group Ltd. */


/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Parameters .. */
/*     .. */

/*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not */
/*     conjugated or transposed, set  CONJA and CONJB  as true if  A  and */
/*     B  respectively are to be  transposed but  not conjugated  and set */
/*     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A */
/*     and the number of rows of  B  respectively. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */
    nota = lsame_(transa, "N");
    notb = lsame_(transb, "N");
    conja = lsame_(transa, "C");
    conjb = lsame_(transb, "C");
    if (nota) {
	nrowa = *m;
	ncola = *k;
    } else {
	nrowa = *k;
	ncola = *m;
    }
    if (notb) {
	nrowb = *k;
    } else {
	nrowb = *n;
    }

/*     Test the input parameters. */

    info = 0;
    if (! nota && ! conja && ! lsame_(transa, "T")) {
	info = 1;
    } else if (! notb && ! conjb && ! lsame_(transb, "T")) {
	info = 2;
    } else if (*m < 0) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*k < 0) {
	info = 5;
    } else if (*lda < max(1,nrowa)) {
	info = 8;
    } else if (*ldb < max(1,nrowb)) {
	info = 10;
    } else if (*ldc < max(1,*m)) {
	info = 13;
    }
    if (info != 0) {
	xerbla_("CGEMM ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || (alpha->r == 0.f && alpha->i == 0.f || *k == 0) 
	    && (beta->r == 1.f && beta->i == 0.f)) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (alpha->r == 0.f && alpha->i == 0.f) {
	if (beta->r == 0.f && beta->i == 0.f) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i__3 = i__ + j * c_dim1;
		    c__[i__3].r = 0.f, c__[i__3].i = 0.f;
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i__3 = i__ + j * c_dim1;
		    i__4 = i__ + j * c_dim1;
		    q__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4].i, 
			    q__1.i = beta->r * c__[i__4].i + beta->i * c__[
			    i__4].r;
		    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
/* L30: */
		}
/* L40: */
	    }
	}
	return 0;
    }

/*     Start the operations. */

    if (notb) {
	if (nota) {

/*           Form  C := alpha*A*B + beta*C. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (beta->r == 0.f && beta->i == 0.f) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = i__ + j * c_dim1;
			c__[i__3].r = 0.f, c__[i__3].i = 0.f;
/* L50: */
		    }
		} else if (beta->r != 1.f || beta->i != 0.f) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = i__ + j * c_dim1;
			i__4 = i__ + j * c_dim1;
			q__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, q__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
			c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
/* L60: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    i__3 = l + j * b_dim1;
		    if (b[i__3].r != 0.f || b[i__3].i != 0.f) {
			i__3 = l + j * b_dim1;
			q__1.r = alpha->r * b[i__3].r - alpha->i * b[i__3].i, 
				q__1.i = alpha->r * b[i__3].i + alpha->i * b[
				i__3].r;
			temp.r = q__1.r, temp.i = q__1.i;
			i__3 = *m;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    i__4 = i__ + j * c_dim1;
			    i__5 = i__ + j * c_dim1;
			    i__6 = i__ + l * a_dim1;
			    q__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i, 
				    q__2.i = temp.r * a[i__6].i + temp.i * a[
				    i__6].r;
			    q__1.r = c__[i__5].r + q__2.r, q__1.i = c__[i__5]
				    .i + q__2.i;
			    c__[i__4].r = q__1.r, c__[i__4].i = q__1.i;
/* L70: */
			}
		    }
/* L80: */
		}
/* L90: */
	    }
	} else if (conja) {

/*           Form  C := alpha*conjg( A' )*B + beta*C. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp.r = 0.f, temp.i = 0.f;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			r_cnjg(&q__3, &a[l + i__ * a_dim1]);
			i__4 = l + j * b_dim1;
			q__2.r = q__3.r * b[i__4].r - q__3.i * b[i__4].i, 
				q__2.i = q__3.r * b[i__4].i + q__3.i * b[i__4]
				.r;
			q__1.r = temp.r + q__2.r, q__1.i = temp.i + q__2.i;
			temp.r = q__1.r, temp.i = q__1.i;
/* L100: */
		    }
		    if (beta->r == 0.f && beta->i == 0.f) {
			i__3 = i__ + j * c_dim1;
			q__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				q__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
			c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
		    } else {
			i__3 = i__ + j * c_dim1;
			q__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				q__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
			i__4 = i__ + j * c_dim1;
			q__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, q__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
			q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
			c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
		    }
/* L110: */
		}
/* L120: */
	    }
	} else {

/*           Form  C := alpha*A'*B + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp.r = 0.f, temp.i = 0.f;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			i__4 = l + i__ * a_dim1;
			i__5 = l + j * b_dim1;
			q__2.r = a[i__4].r * b[i__5].r - a[i__4].i * b[i__5]
				.i, q__2.i = a[i__4].r * b[i__5].i + a[i__4]
				.i * b[i__5].r;
			q__1.r = temp.r + q__2.r, q__1.i = temp.i + q__2.i;
			temp.r = q__1.r, temp.i = q__1.i;
/* L130: */
		    }
		    if (beta->r == 0.f && beta->i == 0.f) {
			i__3 = i__ + j * c_dim1;
			q__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				q__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
			c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
		    } else {
			i__3 = i__ + j * c_dim1;
			q__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				q__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
			i__4 = i__ + j * c_dim1;
			q__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, q__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
			q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
			c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
		    }
/* L140: */
		}
/* L150: */
	    }
	}
    } else if (nota) {
	if (conjb) {

/*           Form  C := alpha*A*conjg( B' ) + beta*C. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (beta->r == 0.f && beta->i == 0.f) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = i__ + j * c_dim1;
			c__[i__3].r = 0.f, c__[i__3].i = 0.f;
/* L160: */
		    }
		} else if (beta->r != 1.f || beta->i != 0.f) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = i__ + j * c_dim1;
			i__4 = i__ + j * c_dim1;
			q__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, q__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
			c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
/* L170: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    i__3 = j + l * b_dim1;
		    if (b[i__3].r != 0.f || b[i__3].i != 0.f) {
			r_cnjg(&q__2, &b[j + l * b_dim1]);
			q__1.r = alpha->r * q__2.r - alpha->i * q__2.i, 
				q__1.i = alpha->r * q__2.i + alpha->i * 
				q__2.r;
			temp.r = q__1.r, temp.i = q__1.i;
			i__3 = *m;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    i__4 = i__ + j * c_dim1;
			    i__5 = i__ + j * c_dim1;
			    i__6 = i__ + l * a_dim1;
			    q__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i, 
				    q__2.i = temp.r * a[i__6].i + temp.i * a[
				    i__6].r;
			    q__1.r = c__[i__5].r + q__2.r, q__1.i = c__[i__5]
				    .i + q__2.i;
			    c__[i__4].r = q__1.r, c__[i__4].i = q__1.i;
/* L180: */
			}
		    }
/* L190: */
		}
/* L200: */
	    }
	} else {

/*           Form  C := alpha*A*B'          + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (beta->r == 0.f && beta->i == 0.f) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = i__ + j * c_dim1;
			c__[i__3].r = 0.f, c__[i__3].i = 0.f;
/* L210: */
		    }
		} else if (beta->r != 1.f || beta->i != 0.f) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = i__ + j * c_dim1;
			i__4 = i__ + j * c_dim1;
			q__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, q__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
			c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
/* L220: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    i__3 = j + l * b_dim1;
		    if (b[i__3].r != 0.f || b[i__3].i != 0.f) {
			i__3 = j + l * b_dim1;
			q__1.r = alpha->r * b[i__3].r - alpha->i * b[i__3].i, 
				q__1.i = alpha->r * b[i__3].i + alpha->i * b[
				i__3].r;
			temp.r = q__1.r, temp.i = q__1.i;
			i__3 = *m;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    i__4 = i__ + j * c_dim1;
			    i__5 = i__ + j * c_dim1;
			    i__6 = i__ + l * a_dim1;
			    q__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i, 
				    q__2.i = temp.r * a[i__6].i + temp.i * a[
				    i__6].r;
			    q__1.r = c__[i__5].r + q__2.r, q__1.i = c__[i__5]
				    .i + q__2.i;
			    c__[i__4].r = q__1.r, c__[i__4].i = q__1.i;
/* L230: */
			}
		    }
/* L240: */
		}
/* L250: */
	    }
	}
    } else if (conja) {
	if (conjb) {

/*           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp.r = 0.f, temp.i = 0.f;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			r_cnjg(&q__3, &a[l + i__ * a_dim1]);
			r_cnjg(&q__4, &b[j + l * b_dim1]);
			q__2.r = q__3.r * q__4.r - q__3.i * q__4.i, q__2.i = 
				q__3.r * q__4.i + q__3.i * q__4.r;
			q__1.r = temp.r + q__2.r, q__1.i = temp.i + q__2.i;
			temp.r = q__1.r, temp.i = q__1.i;
/* L260: */
		    }
		    if (beta->r == 0.f && beta->i == 0.f) {
			i__3 = i__ + j * c_dim1;
			q__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				q__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
			c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
		    } else {
			i__3 = i__ + j * c_dim1;
			q__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				q__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
			i__4 = i__ + j * c_dim1;
			q__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, q__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
			q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
			c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
		    }
/* L270: */
		}
/* L280: */
	    }
	} else {

/*           Form  C := alpha*conjg( A' )*B' + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp.r = 0.f, temp.i = 0.f;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			r_cnjg(&q__3, &a[l + i__ * a_dim1]);
			i__4 = j + l * b_dim1;
			q__2.r = q__3.r * b[i__4].r - q__3.i * b[i__4].i, 
				q__2.i = q__3.r * b[i__4].i + q__3.i * b[i__4]
				.r;
			q__1.r = temp.r + q__2.r, q__1.i = temp.i + q__2.i;
			temp.r = q__1.r, temp.i = q__1.i;
/* L290: */
		    }
		    if (beta->r == 0.f && beta->i == 0.f) {
			i__3 = i__ + j * c_dim1;
			q__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				q__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
			c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
		    } else {
			i__3 = i__ + j * c_dim1;
			q__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				q__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
			i__4 = i__ + j * c_dim1;
			q__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, q__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
			q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
			c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
		    }
/* L300: */
		}
/* L310: */
	    }
	}
    } else {
	if (conjb) {

/*           Form  C := alpha*A'*conjg( B' ) + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp.r = 0.f, temp.i = 0.f;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			i__4 = l + i__ * a_dim1;
			r_cnjg(&q__3, &b[j + l * b_dim1]);
			q__2.r = a[i__4].r * q__3.r - a[i__4].i * q__3.i, 
				q__2.i = a[i__4].r * q__3.i + a[i__4].i * 
				q__3.r;
			q__1.r = temp.r + q__2.r, q__1.i = temp.i + q__2.i;
			temp.r = q__1.r, temp.i = q__1.i;
/* L320: */
		    }
		    if (beta->r == 0.f && beta->i == 0.f) {
			i__3 = i__ + j * c_dim1;
			q__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				q__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
			c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
		    } else {
			i__3 = i__ + j * c_dim1;
			q__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				q__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
			i__4 = i__ + j * c_dim1;
			q__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, q__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
			q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
			c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
		    }
/* L330: */
		}
/* L340: */
	    }
	} else {

/*           Form  C := alpha*A'*B' + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp.r = 0.f, temp.i = 0.f;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			i__4 = l + i__ * a_dim1;
			i__5 = j + l * b_dim1;
			q__2.r = a[i__4].r * b[i__5].r - a[i__4].i * b[i__5]
				.i, q__2.i = a[i__4].r * b[i__5].i + a[i__4]
				.i * b[i__5].r;
			q__1.r = temp.r + q__2.r, q__1.i = temp.i + q__2.i;
			temp.r = q__1.r, temp.i = q__1.i;
/* L350: */
		    }
		    if (beta->r == 0.f && beta->i == 0.f) {
			i__3 = i__ + j * c_dim1;
			q__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				q__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
			c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
		    } else {
			i__3 = i__ + j * c_dim1;
			q__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				q__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
			i__4 = i__ + j * c_dim1;
			q__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, q__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
			q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
			c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
		    }
/* L360: */
		}
/* L370: */
	    }
	}
    }

    return 0;

/*     End of CGEMM . */

} /* cgemm_ */
