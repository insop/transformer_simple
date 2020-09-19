#ifndef _SMLCOMP_H_
#define _SMLCOMP_H_

/* smlcomp.h: Small Matrix library computations  */
/***********************************************************/

/* configuration, data structures and access */
#include "smlbase.h"


REAL sml_hypot (REAL a, REAL b);   /*   sqrt(a*a + b*b)    */
int  Matmini (int a, int b);
int  Matmaxi (int a, int b);
size_t Matminu (size_t a, size_t b);
size_t Matmaxu (size_t a, size_t b);
REAL Matminr (REAL a, REAL b);
REAL Matmaxr (REAL a, REAL b);
double  Matmind (double a, double b);
double  Matmaxd (double a, double b);

void MatSort (MATRIX A);

int MatIsSquare  (MATRIX a);
int MatIsSymmetric(MATRIX a);

void MatI (MATRIX a, size_t n, REAL c); /* a = diag(c), n by n */
void MatJ (MATRIX a, size_t m, size_t n, REAL c); /* a = m by n, filled with c */

void MatFill   (MATRIX a, REAL x);       /* a[,] = x */
void MatFillAD (MATRIX a, REAL x);       /* above diag(a) = x */
void MatFillBD (MATRIX a, REAL x);       /* below diag(a) = x */
void MatFillD  (MATRIX a, REAL x);       /*       diag(a) = x */
void MatBFill  (MATRIX a, size_t i1, size_t i2, size_t j1, size_t j2, REAL x);
                                  /* a[i1:i2, j1:j2] = x */
void MatFillCol  (MATRIX a, size_t j, REAL x);    /* a[, j] <- x */
void MatFillRow  (MATRIX a, size_t i, REAL x);    /* a[i, ] <- x */

void MatCopy    (MATRIX a, MATRIX b);     /* a = b  */
void MatSetDiag (MATRIX a, MATRIX c);     /* diag(a) = c  */
void MatGetDiag (MATRIX a, MATRIX c);     /* a = diag(c)  */
void MatBCopy   (MATRIX a, size_t ii1, size_t jj1,
                 MATRIX b, size_t i1, size_t i2, size_t j1, size_t j2);
                     /* a[ii1:ii2, jj1:jj2] = b[i1:i2, j1:j2] */
void MatTran    (MATRIX a, MATRIX b);     /* a = b' */
void MatL2U     (MATRIX a);        /* lower -> upper 1/2  */
void MatU2L     (MATRIX a);        /* upper -> lower 1/2  */

void MatAdd       (MATRIX a, MATRIX b, MATRIX c);     /* a = b + c */
void MatSub       (MATRIX a, MATRIX b, MATRIX c);     /* a = b - c */
void MatMul       (MATRIX a, MATRIX b, MATRIX c);     /* a = b * c */

void MatWMul (const char* op, MATRIX a, MATRIX b, MATRIX w, MATRIX c);
/*
Weighted matrix multiplication.
The string op[] specifies the operation as follows:
     if      op = "t  "     a <- b' * c
     else if op = "  t"     a <- b  * c'
     else if op = "t t"     a <- b' * c'
     else if op = " w "     a <- b  * diag(w) * c
     else if op = "tw "     a <- b' * diag(w) * c
     else if op = " wt"     a <- b  * diag(w) * c'
     else if op = "twt"     a <- b' * diag(w) * c'
     else                   a <- b  * c
*/

void MatAddDiag   (MATRIX a, MATRIX b, MATRIX c);     /* a = b + diag(c) */
void MatAddScalar (MATRIX a, MATRIX b, REAL   c);     /* a = b + c */
void MatMulScalar (MATRIX a, MATRIX b, REAL   c);     /* a = b * c */
void MatPreLower  (MATRIX a, MATRIX l, MATRIX c);     /* a = l * c */
void MatPreDiag   (MATRIX a, MATRIX b, MATRIX c);     /* a = diag(b) * c */

void MatMulE    (MATRIX a, MATRIX b, MATRIX c); /* a = b * c elementwise */
void MatDivE    (MATRIX a, MATRIX b, MATRIX c); /* a = b / c elementwise */
void MatDivE0 (MATRIX a, MATRIX b, MATRIX c); /* a = b / c   el-wise, 1/0 = 0*/
void MatRecipE  (MATRIX a, MATRIX b);   /* a = 1/b elementwise  */
void MatRecipE0 (MATRIX a, MATRIX b);   /* a = 1/b elementwise, x/0=0 */
void MatSqrtE   (MATRIX a, MATRIX b);     /* a = b^0.5 elementwise  */
void MatExpE    (MATRIX a, MATRIX b);     /* a = exp(b) elementwise  */
void MatLogE    (MATRIX a, MATRIX b);     /* a = log(b) elementwise  */
void MatApplyE  (MATRIX a, REAL (*f)(REAL), MATRIX b); /* a=f(b) elementwise */


REAL   MatMin    (MATRIX a);             /*  min element of a */
REAL   MatMax    (MATRIX a);             /*  max element of a */
REAL   MatMinAbs (MATRIX a);          /*  min |element| of a */
REAL   MatMaxAbs (MATRIX a);          /*  max |element| of a */
REAL2  MatSum    (MATRIX a);           /*  sum of elements of a */
REAL2  MatSumAbs (MATRIX a);        /*  sum of |elements| of a */
REAL2  MatSS     (MATRIX a);           /*  sum of squares of a, ssq(a) */
REAL2  MatSS2    (MATRIX r, MATRIX a); /*  ssq(a) == r[1,1]*r[1,1]*r[1,2] */
REAL2  MatTrace  (MATRIX r);           /*  sum of diagonals */
size_t   Matnzd    (MATRIX a);        /*  # of non-zero diagonals */

void   MatColSum (MATRIX a, MATRIX b);    /*  col sums of b */
void   MatRowSum (MATRIX a, MATRIX b);    /*  row sums of b */
REAL2  MatCSum   (MATRIX y, MATRIX x);/* y <- cumulative sum  of x */
REAL2  MatCMean  (MATRIX y, MATRIX x);/* y <- cumulative mean of x */

REAL2  MatRowxRow  (MATRIX x, size_t i, MATRIX y, size_t j, MATRIX w, int use_w);
  /* x[i,] * y[j,]' if use_w==0, x[i,] * diag(w) * y[j,]' if use_w==1 */

REAL2  MatColxCol  (MATRIX x, size_t i, MATRIX y, size_t j, MATRIX w, int use_w);
  /* x[,i]' * y[,j] if use_w==0, x[,i]' * diag(w) * y[,j] if use_w==1 */

REAL2  MatRowxCol  (MATRIX x, size_t i, MATRIX y, size_t j, MATRIX w, int use_w);
  /* x[i,] * y[,j] if use_w==0, x[i,] * diag(w) * y[,j] if use_w==1 */

void Matpwp (MATRIX a, MATRIX b);
  /* a <- pairwise products of b, separately within each column */


int    Matdaxpy_c (MATRIX z, size_t i, REAL   a, MATRIX x, size_t j,
                   MATRIX y, size_t k); /* z[,i] <- a x[,j] + y[,k] */

/* vector norms */
REAL2 MatCol1norm      (MATRIX A, size_t i); /* ||A[,i]||_1 */
REAL2 MatCol1normalize (MATRIX A, size_t i);
REAL2 MatCol2norm      (MATRIX A, size_t i); /* ||A[,i]||_2 */
REAL2 MatCol2normalize (MATRIX A, size_t i);
REAL2 MatColinorm      (MATRIX A, size_t i); /* ||A[,i]||_infty */
REAL2 MatColinormalize (MATRIX A, size_t i);


/* QR factoring of A, A = Q R, by Householder reflectors */
int MatQR       (MATRIX Q,  MATRIX R, MATRIX A);
int MatQRgetQR  (MATRIX QR, MATRIX Rdiag, MATRIX A);
int MatQRgetQ   (MATRIX Q,  MATRIX QR);
int MatQRgetR   (MATRIX R,  MATRIX QR, MATRIX Rdiag);
int MatQRgetH   (MATRIX H,  MATRIX QR);
int MatQRsolve  (MATRIX X,  MATRIX B, MATRIX QR, MATRIX Rdiag);

/* The Singular Value Decomposition of A, A = U S V'  */
int  MatSVD (MATRIX U, MATRIX S, MATRIX V, MATRIX A);

/*   Eigenvectors and Eigenvalues, A = V D V'  */
int  MatEigen (MATRIX V, MATRIX D_real, MATRIX D_imag, MATRIX A);

int  MatCompanion (MATRIX companion, MATRIX poly_coeff);

/* SML  - Cholesky      */

int  MatCholIP (MATRIX A);               /* Cholesky, in-situ */
int  MatChol   (MATRIX L, MATRIX A);     /* L <- root(A), A = L L' */
int  MatSInvIP (MATRIX A);
                           /* g-inverts a symmetric non-negative definite matrix in-situ */
int  MatSHInvIP(MATRIX A);  /* lower A <- g-inverse of the Cholesky root of A  */
int  MatLtLIP  (MATRIX L);  /*  L'L in-situ,  L is a square lower triangular matrix */
int  MatLLtIP  (MATRIX L);  /*  L'L in-situ,  L is a square lower triangular matrix */
int  MatInvLIP (MATRIX L);  /* g-invert a square lower triangular matrix L in situ */
int  MatCholSolve (MATRIX X, MATRIX L, MATRIX B);
    /*  Solve a linear system A*X = B, using the previously computed
            cholesky factorization of A, A = L L'  */


#endif
