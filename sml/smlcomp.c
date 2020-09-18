/* smlcomp.c: Small Matrix Library Computations  */
/* clean version, use accessor functions only */
/* REAL and REAL2 are the only floating point types used */
/**********************************************************/

/* stdio.h not needed */
#include <stdio.h>

#include <math.h>
#include <string.h>

#include "smlcomp.h"

#define SML_NAN  ((REAL)0.0/(REAL)0.0)

/************************************************************/

#if 0

REAL sml_fabs (REAL a)
{
  return  (a >= 0) ? a : -a;
}

#endif

/************************************************************/
REAL sml_hypot (REAL a, REAL b)  // sqrt(a*a + b*b)
/*
Effective precision is that of sml_sqrt(REAL)
*/

{
  REAL aa = (a >= 0) ? a : -a;    // |a|
  REAL bb = (b >= 0) ? b : -b;    // |b|

  if (aa > bb) {
    REAL r = bb / aa;
    return aa*sml_sqrt(1.0 + r*r);
  }
  else if (bb != 0.0 ) {
    REAL r = aa / bb;
    return bb*sml_sqrt(1.0 + r*r);
  }

  return 0.0; // a == b == 0
}
/*----------------------------------------------*/
int  MatIsSquare  (MATRIX a)  {return (MatRows(a)==MatCols(a));}
int  Matmini (int a, int b)   {return ((a) < (b) ? (a) : (b));}
int  Matmaxi (int a, int b)   {return ((a) > (b) ? (a) : (b));}
size_t Matminu (size_t a, size_t b) {return ((a) < (b) ? (a) : (b));}
size_t Matmaxu (size_t a, size_t b) {return ((a) > (b) ? (a) : (b));}
REAL Matminr (REAL a, REAL b) {return ((a) < (b) ? (a) : (b));}
REAL Matmaxr (REAL a, REAL b) {return ((a) > (b) ? (a) : (b));}
double Matmind (double a, double b) {return ((a) < (b) ? (a) : (b));}
double Matmaxd (double a, double b) {return ((a) > (b) ? (a) : (b));}
/*----------------------------------------------*/
void MatI (MATRIX a, size_t n, REAL c)  /* a = diag(c), n by n */
{
  MatReDim(a, n, n);
  if (MatErr()) return;
  MatFill (a, 0);
  MatFillD(a, c);
}
/*----------------------------------------------*/
void MatJ (MATRIX a, size_t m, size_t n, REAL c)  /* a = m by n, filled with c */
{
  MatReDim(a, m, n);
  if (MatErr()) return;
  MatFill(a, c);
}
/*----------------------------------------------*/


int  MatIsSymmetric(MATRIX a)
  /* returns 1 if square and symmetric, 0 otherwise */
{

  size_t i, j, m = MatRows(a), n = MatCols(a);

  if (m != n)  return 0;

  for (i = 0; i != m; i++)
    for (j = 0; j != i; j++)
      if (MatGet0(a, i, j) != MatGet0(a, j, i)) return 0;

  return 1;

}

/*----------------------------------------------*/

size_t  Matnzd (MATRIX a)
/*  # of non-zero diagonals,  a is a square matrix */
{
  size_t i,  result = 0;

  if (!MatIsSquare(a))  {
    MatErrThrow("Matnzd: not square matrix.");
    return  -1;
  }
  for (i = 0; i != MatRows(a); ++i)  result += (MatGet0(a, i, i) != 0);

  return result;
}

/*----------------------------------------------*/
REAL2  MatTrace  (MATRIX a)
{
  size_t i, j, m = MatRows(a), n = MatCols(a);

  REAL2  s=0;

  if (m==0 || n==0) return 0;
  n = n <= m ? n : m;

  for (i = 0; i != n; ++i)
      s += MatGet0(a, i, i);

  return s;
}

/*----------------------------------------------*/
void   MatRowSum (MATRIX a, MATRIX b)

/*  a <- row sums of b.   a == b allowed */

{
  size_t i, j, brows = MatRows(b), bcols = MatCols(b);
  MATRIX t = MatDim(brows, 1);

  if (MatErr()) return;

  for (i=0; i != brows; ++i)  {
    REAL2  s=0;
    for (j=0; j != bcols; ++j)   s += MatGet0(b, i, j);
    MatSet0(t, i, 0, s);
  }
  MatSwap(a, t);
  MatUnDim(t);

}

/*----------------------------------------------*/
void   MatColSum (MATRIX a, MATRIX b)

/*  a <- col sums of b.   a == b allowed */
{
  size_t i, j, brows = MatRows(b), bcols = MatCols(b);
  MATRIX t = MatDim(1, bcols);

  if (MatErr()) return;

  for (j=0; j != bcols; ++j)  {
    REAL2  s=0;
    for (i=0; i != brows; ++i)   s += MatGet0(b, i, j);
    MatSet0(t, 0, j, s);
  }
  MatSwap(a, t);
  MatUnDim(t);
}

/*----------------------------------------------*/
void MatFillAD(MATRIX a, REAL x)
/* above and right of diag(a) = x */
{
  size_t i, j, cols = MatCols(a), r;

  r = Matminu(MatRows(a), cols);

  for (i = 0; i != r; ++i)
    for (j = i+1; j != cols; ++j)
      MatSet0(a, i, j, x);
}
/*----------------------------------------------*/
void MatFillBD(MATRIX a, REAL x)
/* below diag(a) = x */
{
  size_t i, j, rows = MatRows(a), r;

  r = Matminu(MatCols(a), rows);
  for (j = 0; j != r; ++j)
    for (i = j+1; i != rows; ++i)
      MatSet0(a, i, j, x);
}
/*----------------------------------------------*/
void MatFillD (MATRIX a, REAL x)        /*       diag(a) = x */
{
  size_t i, n;

  n = Matminu (MatRows(a), MatCols(a));
  for (i=0; i != n; ++i)   MatSet0(a, i, i, x);
}
/*----------------------------------------------*/

void MatBFill  (MATRIX a, size_t i1, size_t i2, size_t j1, size_t j2, REAL x)
    /* a[i1:i2, j1:j2] <- x */
   /* fill block rows i1-i2, cols j1-j2 of a with x */
{
  size_t i, j, error;

  error = (i1 > i2) || (j1 > j2) ||
          (i1 < 1) || (i2 > MatRows(a)) ||
          (j1 < 1) || (j2 > MatCols(a));

  if (error)  {
    MatErrThrow("MatBFill: invalid arguments.");
    return;
  }

  for (i=i1; i <= i2; i++)
    for (j=j1; j <= j2; j++)
      MatSet(a, i, j, x);

}

/*----------------------------------------------*/

void MatFillCol  (MATRIX a, size_t j, REAL x)
    /* a[, j] <- x,   fill column j with value x */
{
  size_t i, rows, error;

  error =  (j < 1 || j > MatCols(a));

  if (error)  {
    MatErrThrow("MatFillCol: invalid arguments.");
    return;
  }

  --j;
  rows = MatRows(a);
  for (i=0; i != rows; i++)  MatSet0(a, i, j, x);

}

/*----------------------------------------------*/

void MatFillRow  (MATRIX a, size_t i, REAL x)
    /* a[i, ] <- x,   fill row i with value x */
{
  size_t j, cols, error;

  error =  (i < 1 || i > MatRows(a));

  if (error)  {
    MatErrThrow("MatFillRow: invalid arguments.");
    return;
  }

  --i;
  cols = MatCols(a);
  for (j=0; j != cols; j++)  MatSet0(a, i, j, x);

}

/*----------------------------------------------*/

void MatSetDiag   (MATRIX a, MATRIX c)      /* diag(a) <- c */
  /* a must be a square matrix */
  /* c must be a column vector */
{
  size_t i, r;

  if (!MatIsSquare(a)) {
    MatErrThrow("MatSetDiag: nonsquare matrix.");
    return;
  }

  if ((MatRows(c) != MatRows(a)) || (MatCols(c) != 1)) {
    MatErrThrow("MatSetDiag: invalid dims.");
    return;
  }

  if (a == c) {
    MatErrThrow("MatSetDiag: same matrix in both arguments.");
    return;
  }

  r = MatRows(c);
  for (i = 0; i != r ; ++i)
      MatSet0(a, i, i, MatGet0(c, i, 0));
}

/*----------------------------------------------*/

void MatGetDiag   (MATRIX a, MATRIX c)      /* a <- diag(c) */
  /* c must be a square matrix */
  /* a <- column vector */
{
  size_t i, r;
  MATRIX tmp;

  if (!MatIsSquare(c)) {
    MatErrThrow("MatGetDiag: nonsquare matrix.");
    return;
  }

  r = MatRows(c);
  tmp = MatDim(r, 1);
  if (MatErr()) return;

  for (i = 0; i != r; ++i)
      MatSet0(tmp, i, 0, MatGet0(c, i, i));

  MatSwap (a, tmp);
  MatUnDim (tmp);
}

/*----------------------------------------------*/

void Matpwp (MATRIX a, MATRIX b)
  /* a <- pairwise products of b, separately within each column */
  /* products are ordered, e.g. n=4: 1*2, 1*3, 1*4, 2*3, 2*4, 3*4 */
  /* if b has n rows, a will have "n choose 2" rows */
{
  size_t i1, i2, j, k, n, p;
  MATRIX tmp;

  n = MatRows(b);
  p = MatCols(b);

  if (n < 2 || p == 0) {
    MatReDim(a, 0, 0);
    return;
  }

  /* Now: n >= 2 && p >= 1 */
  tmp = MatDim( (n*(n-1))/2, p);
  if (MatErr()) return;

  k = 1;
  for (i1 = 0; i1 <= n-2; ++i1)
    for (i2 = 1 + i1; i2 < n ; ++i2, ++k)
      for (j = 0; j < p; ++j)
        MatSet0(tmp, k, j, MatGet0(b, i1, j) * MatGet0(b, i2, j) );

  MatSwap (a, tmp);
  MatUnDim(tmp);
}

/*----------------------------------------------*/

void MatBCopy (MATRIX a, size_t ii1, size_t jj1,
               MATRIX b, size_t i1,  size_t i2, size_t j1, size_t j2)
/*
 Matrix Block Copy
   a[ii1:?, jj1:?] <- b[i1:i2, j1:j2]
 ? are determined so as to copy as much as possible
without running out of the bounds of matrix a
a == b raises the error flag
*/
{
  size_t i, ii2, jj2, error;

  if (a == b)  {
    MatErrThrow("MatBCopy: target == destination not allowed.");
    return;
  }

  error = (i1 > i2) || (j1 > j2) ||
          (i1 < 1) || (i2 > MatRows(b)) ||
          (j1 < 1) || (j2 > MatCols(b)) ||
          (ii1 < 1) || (ii1 > MatRows(a)) ||
          (jj1 < 1) || (jj1 > MatCols(a));

  if (error)  {
    MatErrThrow("MatBCopy: invalid arguments.");
    return;
  }

  ii2 = Matminu(ii1+i2-i1, MatRows(a));
  jj2 = Matminu(jj1+j2-j1, MatCols(a));

  --ii1; --ii2;
  --jj1; --jj2;
  --i1; --i2;
  --j1; --j2; /* 0-based loop */
  for (i = ii1; i <= ii2; ++i, ++i1)  {
    size_t j = jj1, t = j1;
    for ( ; j <= jj2; ++j, ++t)
      MatSet0(a, i, j, MatGet0(b, i1, t));
  }

}

/*----------------------------------------------*/

void MatAddDiag   (MATRIX a, MATRIX b, MATRIX c)      /* a <- b + diag(c) */
  /* b must be a square matrix */
  /* c must be a column vector */
{
  size_t i, brows=MatRows(b);

  if (!MatIsSquare(b))   {
    MatErrThrow("MatAddDiag: nonsquare matrix.");
    return;
  }

  if ((MatRows(c) != brows) || (MatCols(c) != 1))   {
    MatErrThrow("MatAddDiag: invalid dims.");
    return;
  }

  if (a != b)  { MatCopy(a, b); if (MatErr()) return; }

  for (i = 0; i != brows; ++i)
        MatSetPE0(a, i, i, MatGet0(c, i, 0));
}

/*----------------------------------------------*/

static void MatMul2 (MATRIX a, MATRIX b, MATRIX c)      /* a = b * c */
  /* assumes a!=b, a!=c and cols(b)==rows(c) */
  /* assumes valid dimensions of a, b, c */
{
  size_t brows, bcols, ccols, i, k, j;

  brows = MatRows(b);
  bcols = MatCols(b);
  ccols = MatCols(c);

  if (b != c)
  for (i = 0; i != brows; ++i)
    for (j = 0; j != ccols; ++j) {
      REAL2  s = 0;
      for (k = 0; k != bcols; ++k)
        s += MatGet0(b, i, k) * MatGet0(c, k, j);
      MatSet0(a, i, j, s);
    }

  else
  for (i = 0; i != brows; ++i)
    for (j = 0; j != ccols; ++j) {
      REAL2  s = 0;
      for (k = 0; k != bcols; ++k)
        s += MatGet0(b, i, k) * MatGet0(b, k, j);
      MatSet0(a, i, j, s);
    }

}

/*----------------------------------------------*/

void MatMul (MATRIX a, MATRIX b, MATRIX c)      /* a = b * c */
  /* a==b or a==c are OK, e.g.  MatMul(a,a,a) */
{
  MATRIX tmp;

  if (MatCols(b) !=  MatRows(c)) {
    MatErrThrow("MatMul: unconformable matrices.");
    return;
  }

  tmp = MatDim(MatRows(b), MatCols(c));
  if (MatErr()) return;

  MatMul2 (tmp, b, c);
  MatSwap (a, tmp);
  MatUnDim(tmp);
} /* MatMul */

/*----------------------------------------------*/

static void MatPreLower2 (MATRIX a, MATRIX b, MATRIX c)  /* a = b * c */
  /* assumes a<>b, a<>c and cols(b)=rows(c) */
  /* assumes valid dimensions of a, b, c */
  /* assumes b is square lower triangular */
{
  size_t brows, ccols, i, k, j;

  brows = MatRows(b);
  ccols = MatCols(c);
  for (i=0; i != brows; ++i)
    for (j=0; j != ccols; ++j) {
      REAL2  s = 0;
      for (k=0; k <= i; k++) s += MatGet0(b, i, k) * MatGet0(c, k, j);
      MatSet0(a, i, j, s);
    }

}

/*----------------------------------------------*/

void MatPreLower (MATRIX a, MATRIX b, MATRIX c)      /* a = b * c */
  /* may be called even if a==b or a==c, e.g.  MatMul(a,a,a) */
  /* b must be square and is assumed lower triangular */
{
  MATRIX tmp;

  if (MatCols(b) !=  MatRows(c)) {
    MatErrThrow("MatPreLower: unconformable matrices.");
    return;
  }

  if (!MatIsSquare(b))  {
    MatErrThrow("MatPreLower: not square matrix.");
    return;
  }

  tmp = MatDim(MatRows(b), MatCols(c));
  if (MatErr()) return;

  MatPreLower2 (tmp, b, c);
  MatSwap (a, tmp);
  MatUnDim(tmp);
}

/*----------------------------------------------*/

static void MatPreDiag2 (MATRIX a, MATRIX b, MATRIX c)
  /* a = diag(b) * c */
  /* assumes valid dims */
{
  size_t  brows, ccols, i, j;

  brows = MatRows(b);
  ccols = MatCols(c);
  if (a == c) {
    for (i=0; i != brows; ++i) {
      REAL w = MatGet0(b, i, 0);
      for (j=0; j != ccols; ++j)
        MatSet0(a, i, j, w * MatGet0(c, i, j));
    }
  }
  else  {
    for (i=0; i != brows; ++i) {
      REAL w = MatGet0(b, i, 0);
      for (j=0; j != ccols; ++j)
        MatSetTE0(a, i, j, w);
    }
  }

}
/*----------------------------------------------*/

void MatPreDiag (MATRIX a, MATRIX b, MATRIX c)   /* a <- diag(b) * c */
  /* b must be a column vector */
{
  MATRIX tmp;
  int ok = (MatRows(b)==MatRows(c) && MatCols(b)==1);

  if (!ok) {
    MatErrThrow("MatPreDiag: invalid dims.");
    return;
  }

  tmp = MatDim(MatRows(c), MatCols(c));
  if (MatErr()) return;

  MatPreDiag2 (tmp, b, c);
  MatSwap(a, tmp);
  MatUnDim(tmp);
}

/*----------------------------------------------*/

static void MatPostDiag2 (MATRIX a, MATRIX b, MATRIX c)
  /* a = b * diag(c) */
  /* assumes valid dims */
{
  size_t brows, bcols, i, j;

  brows = MatRows(b);
  bcols = MatCols(b);
  if (a == b) {
    for (j=0; j != bcols; ++j) {
      REAL   w = MatGet0(c, j, 0);
      for (i=0; i != brows; ++i)
        MatSet0(a, i, j, w * MatGet0(b, i, j));
    }
  }
  else {
    for (j=0; j != bcols; ++j) {
      REAL   w = MatGet0(c, j, 0);
      for (i=0; i != brows; ++i)
        MatSetTE0(a, i, j, w);
    }
  }

}
/*----------------------------------------------*/

void MatPostDiag (MATRIX a, MATRIX b, MATRIX c)   /* a <- b * diag(c) */
  /* c must be a column vector */
{
  MATRIX tmp;

  if (MatCols(b) !=  MatRows(c)  ||  MatCols(c) !=  1 ) {
    MatErrThrow("MatPostDiag: invalid dims.");
    return;
  }

  tmp = MatDim(MatRows(b), MatCols(b));
  if (MatErr()) return;

  MatPostDiag2 (tmp, b, c);
  MatSwap(a, tmp);
  MatUnDim(tmp);
}

/*----------------------------------------------*/
void MatTran (MATRIX xt, MATRIX x)      /* xt = x'  */
   /* xt == x allowed */
{
  size_t i, j, xrows, xcols;
  MATRIX tmp;

  xcols = MatCols(x);
  xrows = MatRows(x);

  /* same square matrix */
  if (xt == x && xrows == xcols)  {
    for (i = 0; i != xrows; ++i)
      for (j = 0; j < i; ++j)  {
        REAL xij = MatGet0(x, i, j);
        MatSet0(x, i, j, MatGet0(x, j, i));
        MatSet0(x, j, i, xij);
      }
    return;
  }

  /* else  tmp copy then assign to xt  */
  tmp = MatDim(xcols, xrows);
  if (MatErr()) return;;
  for (i = 0; i != xrows; ++i)
    for (j = 0; j != xcols; ++j)
      MatSet0(tmp, j, i,  MatGet0(x, i, j));
  MatSwap (xt, tmp);
  MatUnDim(tmp);
}
/*----------------------------------------------*/

void MatU2L (MATRIX a)
  /* upper -> lower 1/2 , a must be symmetric */
{
  size_t  i, j, n;

  if (!MatIsSquare(a))  {
    MatErrThrow("MatU2L: not square matrix.");
    return;
  }
  n = MatRows(a);

  for (i=0; i != n; ++i)
    for (j=i; j != n; ++j)
      MatSet0(a, j, i,  MatGet0(a, i, j));


} /* MatU2L  */

/*----------------------------------------------*/

void MatL2U (MATRIX a)
  /* lower -> upper 1/2 , a must be symmetric */
{
  size_t i, j, n;

  if (!MatIsSquare(a))  {
    MatErrThrow("MatL2U: not square matrix.");
    return;
  }
  n = MatRows(a);

  for (i=0; i != n; ++i)
    for (j=i; j != n; ++j)
      MatSet0(a, i, j, MatGet0(a, j, i));

} /* MatL2U  */
/*----------------------------------------------*/
static REAL2 MatRowxRow1 (MATRIX x, size_t i, MATRIX y, size_t j, MATRIX w, int use_w)
/*
   returns inner product of i-th row of x
   and the j-th row of y, weighted with w if use_w==1,
   dims assumed correct (i and j are 1-based).
*/
{
  size_t k, xcols = MatCols(x);
  REAL2  s=0;


  if (i==0 || j==0) return SML_NAN;
  --i; --j;
  if (use_w)
    for (k=0; k != xcols; ++k)
      s += MatGet0(x, i, k)*MatGet0(y, j, k)*MatGet0(w, k, 0);
  else
    for (k=0; k != xcols; ++k)
      s += MatGet0(x, i, k)*MatGet0(y, j, k);

  return s;
}
/*----------------------------------------------*/
REAL2  MatRowxRow (MATRIX x, size_t i, MATRIX y, size_t j, MATRIX w, int use_w)
/*
x[i,] * y[j,]` if use_w==0,
x[i,] * diag(w) * y[j,]` if use_w==1.
if use_w==1, w should be a column vector
*/
{
  size_t err, xcols = MatCols(x);

  err = (MatCols(y)!= xcols);
  if (use_w)  err = (err ||  MatRows(w) != xcols || MatCols(w) != 1);

  if (err) {
    MatErrThrow("MatRowxRow: different vector lengths.");
    return 0;
  }

  return MatRowxRow1(x, i, y, j, w, use_w);
}
/*----------------------------------------------*/
static REAL2  MatColxCol1 (MATRIX x, size_t i, MATRIX y, size_t j, MATRIX w, int use_w)
/*
returns inner product of i-th column of x
and the j-th column of y, weighted with w if use_w==1,
dims assumed correct
*/
{
  size_t  k, xrows=MatRows(x);
  REAL2   s=0;

  if (i==0 || j==0) return SML_NAN;
  --i; --j;
  if (use_w)
    for (k=0; k != xrows; k++)
      s += MatGet0(x, k, i)*MatGet0(y, k, j)*MatGet0(w, k, 0);
  else
    for (k=0; k != xrows; k++)
      s += MatGet0(x, k, i)*MatGet0(y, k, j);

  return s;
}
/*----------------------------------------------*/
REAL2   MatColxCol (MATRIX x, size_t i, MATRIX y, size_t j, MATRIX w, int use_w)
/*
x[,i]' * y[,j] if use_w==0, x[,i]' * diag(w) * y[,j] if use_w==1.
w must be a column
*/
{
  size_t err, xrows = MatRows(x);

  err = (MatRows(y)!=xrows);
  if (use_w)  err = (err ||  MatRows(w) != xrows || MatCols(w) != 1);
  if (err) {
    MatErrThrow("MatColxCol: different vector lengths.");
    return 0;
  }

  return MatColxCol1(x, i, y, j, w, use_w);
}
/*----------------------------------------------*/
static REAL2 MatRowxCol1 (MATRIX x, size_t i, MATRIX y, size_t j, MATRIX w, int use_w)
/*
   returns inner product of i-th row of x
   and the j-th col of y, weighted with w if use_w==1,
   dims assumed correct
*/
{
  size_t k, xcols = MatCols(x);
  REAL2  s=0;

  if (i==0 || j==0) return SML_NAN;
  --i; --j;
  if (use_w)
    for (k=0; k != xcols; k++)
      s += MatGet0(x, i, k)*MatGet0(y, k, j)*MatGet0(w, k, 0);
  else
    for (k=0; k != xcols; k++)
      s += MatGet0(x, i, k)*MatGet0(y, k, j);

  return s;
}
/*----------------------------------------------*/
REAL2  MatRowxCol (MATRIX x, size_t i, MATRIX y, size_t j, MATRIX w, int use_w)
/*
x[i,] * y[,j]  if use_w==0,
x[i,] * diag(w) * y[,j]  if use_w==1.
if use_w==1, w should be a column vector
*/
{
  size_t err, xcols = MatCols(x);

  err = (MatRows(y) != xcols);
  if (use_w)  err = (err ||  MatRows(w) != xcols || MatCols(w) != 1);

  if (err) {
    MatErrThrow("MatRowxRow: different vector lengths.");
    return 0;
  }

  return MatRowxRow1(x, i, y, j, w, use_w);
}
/*-------------here---------------------------------*/
/*----------------------------------------------*/

static void MatXtWY1 (MATRIX a, MATRIX x, MATRIX w, MATRIX y, int use_w)
   /*  assumes  dims correct,  a != x, a != y */
   /*  handles  x == y efficiently  */
{
  int xcols, ycols;
  int i, j;

  ycols = MatCols(y);
  xcols = MatCols(x);

  if (x == y) {
    /* compute lower half only */
    for (i=1; i <= xcols; i++)
      for (j=1; j <= i; j++)
        MatSet (a, i, j, MatColxCol1 (x, i, x, j, w, use_w));
    MatL2U (a);   /* reflect to upper half */
  }
  else  {
    for (i=1; i <= xcols; i++)
      for (j=1; j <= ycols; j++)
        MatSet(a, i, j, MatColxCol1 (x, i, y, j, w, use_w));
  }
}
/*----------------------------------------------*/

static void MatXtWY (MATRIX a, MATRIX x, MATRIX w, MATRIX y, int use_w)
   /* if use_w==0  :   a = x' * y
      if use_w==1  :   a = x' * diag(w) * y
    w must be a column vector */
{
  size_t   xrows, xcols, ycols, ok;
  MATRIX tmp;

  xrows = MatRows(x);
  xcols = MatCols(x);
  ycols = MatCols(y);

  ok = (xrows == MatRows(y));
  if (use_w)  ok = ok &&
       (xrows == MatRows(w)  && MatCols(w) == 1);

  if (!ok)  {
    MatErrThrow("MatXtWY: not conformable.");
    return;
  }

  tmp = MatDim(xcols, ycols);
  if (MatErr()) return;

  MatXtWY1 (tmp, x, w, y, use_w);

  MatSwap (a, tmp);
  MatUnDim(tmp);
}
/*----------------------------------------------*/

static void MatXWYt1 (MATRIX a, MATRIX x, MATRIX w, MATRIX y, int use_w)
   /*  assumes  dims correct,  a != x, a != y */
   /*  handles  x == y efficiently  */
{
  size_t i, j, xrows, yrows;

  xrows = MatRows(x);
  yrows = MatRows(y);

  if (x == y) {
    /* compute lower half only */
    for (i=1; i <= xrows; i++)
      for (j=1; j <= i; j++)
        MatSet (a, i, j, MatRowxRow1 (x, i, x, j, w, use_w));

    MatL2U (a);   /* reflect to upper half */
  }
  else {
    for (i=1; i <= xrows; i++)
      for (j=1; j <= yrows; j++)
        MatSet(a, i, j, MatRowxRow1 (x, i, y, j, w, use_w));
  }
}
/*----------------------------------------------*/

static void MatXWYt (MATRIX a, MATRIX x, MATRIX w, MATRIX y, int use_w)
   /* if use_w==0  :   a = x * y'
      if use_w==1  :   a = x * diag(w) * y'
    w must be a column vector */
{
  size_t xrows, xcols, ycols, ok;
  MATRIX tmp;

  xrows = MatRows(x);
  xcols = MatCols(x);
  ycols = MatCols(y);

  ok = (xcols == ycols);
  if (use_w)  ok = ok &&
       (xcols == MatRows(w)   && MatCols(w) == 1);

  if (!ok)  {
    MatErrThrow("MatXWYt: not conformable.");
    return;
  }

  tmp = MatDim(xrows, MatRows(y));
  if (MatErr()) return;

  MatXWYt1 (tmp, x, w, y, use_w);

  MatSwap (a, tmp);
  MatUnDim(tmp);
}
/*----------------------------------------------*/

static void MatXWY1 (MATRIX a, MATRIX x, MATRIX w, MATRIX y, int use_w)
   /*  assumes  dims correct,  a != x, a != y */
{
  int xrows, ycols;
  int i, j;

  xrows = MatRows(x);
  ycols = MatCols(y);

  for (i=1; i <= xrows; i++)
    for (j=1; j <= ycols; j++)
      MatSet(a, i, j, MatRowxRow1 (x, i, y, j, w, use_w));
}
/*----------------------------------------------*/

static void MatXWY (MATRIX a, MATRIX x, MATRIX w, MATRIX y, int use_w)
   /* if use_w==0  :   a = x * y
      if use_w==1  :   a = x * diag(w) * y
    w must be a column vector */
{
  size_t   xrows, xcols, yrows, ok;
  MATRIX tmp;

  xrows = MatRows(x);
  xcols = MatCols(x);
  yrows = MatRows(y);

  ok = (xcols == yrows);
  if (use_w)  ok = ok &&
       (xcols == MatRows(w)  && MatCols(w) == 1);

  if (!ok)  {
    MatErrThrow("MatXWY: not conformable.");
    return;
  }

  tmp = MatDim(xrows, MatCols(y));
  if (MatErr()) return;

  MatXWY1 (tmp, x, w, y, use_w);

  MatSwap (a, tmp);
  MatUnDim(tmp);
}
/*----------------------------------------------*/

static void MatXtWYt (MATRIX a, MATRIX x, MATRIX w, MATRIX y, int use_w)
   /* if use_w==0  :   a = x' * y'
      if use_w==1  :   a = x' * diag(w) * y'
    w must be a column vector */
{
  MatXWY (a, y, w, x, use_w);
  if (MatErr()) return;
  MatTran (a, a);
}
/*----------------------------------------------*/

void MatWMul (const char* op, MATRIX a, MATRIX b, MATRIX w, MATRIX c)
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
{

  if (op == NULL || strlen(op) != 3) {MatMul(a, b, c); return;}

  if (op[1] == ' ')   {
    if (op[0] == ' ')  {
      if (op[2] == ' ')      {MatMul(a, b, c);        return;}
      else if (op[2] == 't') {MatXWYt(a, b, 0, c, 0); return;}
    }
    else if (op[0] == 't')  {
      if (op[2] == ' ')      {MatXtWY (a, b, 0, c, 0); return;}
      else if (op[2] == 't') {MatXtWYt(a, b, 0, c, 0); return;}
    }
  }
  else if (op[1] == 'w')   {
    if (op[0] == ' ')  {
      if (op[2] == ' ')      {MatXWY (a, b, w, c, 1); return;}
      else if (op[2] == 't') {MatXWYt(a, b, w, c, 1); return;}
    }
    else if (op[0] == 't')  {
      if (op[2] == ' ')      {MatXtWY (a, b, w, c, 1); return;}
      else if (op[2] == 't') {MatXtWYt(a, b, w, c, 1); return;}
    }
  }

  /* all other values of op[] */
  MatMul(a, b, c);

}
/*----------------------------------------------*/
REAL2  MatCol1norm (MATRIX a, size_t j)
/*
Returns the 1-norm of a[,j]
*/
{
  size_t   i;
  REAL2  s=0;


  if (j == 0) return SML_NAN;
  --j;
  for (i = 0; i != MatRows(a);  ++i) s += sml_fabs(MatGet0(a, i, j));

  return s;
}
/*----------------------------------------------*/
REAL2  MatCol1normalize (MATRIX a, size_t j)
/*
Normalizes a[,j] to length 1 (in the 1-norm), provided it is
a nonzero column.
Returns the 1-norm of the original column.
*/
{
  REAL2 s2;
  REAL  s;

  if (j == 0) return SML_NAN;
  s = s2 = MatCol1norm (a, j);
  --j;
  if (s>0) {
    size_t k;
    for (k = 0 ; k != MatRows(a); ++k)  MatSetDE0(a,k,j, s);
  }
  return s2;
}
/*----------------------------------------------*/
REAL2  MatCol2norm (MATRIX a, size_t j)
/*
Returns the 2-norm of a[,j] = sqrt(ssq(a[,j]))
*/
{
  REAL  scale=0;
  REAL2 ssq=1;
  size_t i, nrows;

  if (j == 0) return SML_NAN; else --j;
  nrows  = MatRows(a);
  if (nrows == 0) return 0;
  if (nrows == 1) return sml_fabs(MatGet0(a, 0, j));

  for (i = 0; i != nrows; ++i)  {
    REAL absai;
    if ((absai = sml_fabs(MatGet0(a, i, j))) != 0) {
      if (scale < absai) {
        REAL   c = scale/absai;
        ssq      = 1.0 + ssq*c*c;
        scale    = absai;
      }
      else  {
        REAL   c = absai/scale;
        ssq     += c*c;
      }
    }
  }
  return scale * (REAL2) sml_sqrt(ssq);
}
/*----------------------------------------------*/
REAL2  MatCol2normalize (MATRIX a, size_t j)
/*
Normalizes a[,j] to length 1 (in the 2-norm), provided it is
a nonzero column.
Returns the 2-norm of the original column.
*/
{
  REAL2  s2;
  REAL   s;

  if (j == 0) return SML_NAN;
  s = s2 = MatCol2norm (a, j);
  --j;
  if (s>0) {
    size_t k;
    for (k = 0 ; k != MatRows(a); ++k)  MatSetDE0(a,k,j, s);
  }
  return s2;
}

/*----------------------------------------------*/
REAL2  MatColinorm (MATRIX a, size_t j)
/*
Returns the infty-norm of a[,j]
*/
{
  size_t   i;
  REAL s=0;


  if (j == 0) return SML_NAN;
  --j;
  for (i=0; i != MatRows(a);  ++i)
     if (sml_fabs(MatGet0(a, i, j)) > s)  s = sml_fabs(MatGet0(a, i, j));

  return s;
}
/*----------------------------------------------*/
REAL2  MatColinormalize (MATRIX a, size_t j)
/*
Normalizes a[,j] to length 1 (in the infty-norm), provided it is
a nonzero column.
Returns the infty-norm of the original column.
*/
{
  REAL2  s2;
  REAL   s;

  if (j == 0) return SML_NAN;
  s = s2 = MatColinorm (a, j);
  --j;
  if (s>0) {
    size_t k;
    for (k = 0 ; k != MatRows(a); ++k)  MatSetDE0(a,k,j, s);
  }
  return s2;
}
/*----------------------------------------------*/
int  Matdaxpy_c (MATRIX z, size_t i, REAL a,
                 MATRIX x, size_t j,
                 MATRIX y, size_t k)
/*
z[,i] <- a x[,j] + y[,k]
z <- a x + y,  for columns of matrices.
*/
{
  size_t h, nrows = MatRows(z);

  if (MatRows(x)!=nrows || MatRows(y)!=nrows) {
    MatErrThrow("Matdaxpy_c: different vector lengths.");
    return 1;
  }

  if (i < 1 || i > MatCols(z) ||
      j < 1 || j > MatCols(x) ||
      k < 1 || k > MatCols(y) )  {
    MatErrThrow("Matdaxpy_c: out of range colum number(s).");
    return 2;
  }

  --i; --j; --k;        /* we'll run a 0-based loop */
  for (h=0; h < nrows; h++)
    MatSet0(z, h, i, a * MatGet0(x,h,j) + MatGet0(y,h,k));

  return 0;
}
/*----------------------------------------------*/
int  MatCompanion (MATRIX comp, MATRIX a)
/*
input: a = a column of n+1 polynomial coefficients a_0, ..., a_n
The polynomial is
        p(z) = a_0 + a_1 z + a_2 z^2 + ... + a_n z^n
comp <- the companion matrix (n by n)
*/
{

  size_t i, m, n;
  REAL an;

  if (comp == a) return 1;
  if (MatCols(a) != 1) return 2;

  m =  MatRows(a);
  if (m <= 1)  return 3;

  n = m - 1;
  MatReDim(comp, n, n);
  if (MatErr()) return 4;

  MatFill(comp, 0);
  an = MatGet0(a, n, 0);
  for (i = 1; i < n; i++) {
      MatSet0(comp, i, i-1, 1.0);
      MatSet0(comp, i, n-1, -MatGet0(a, i, 0)/an);
  }
  MatSet0(comp, 0, n-1, -MatGet0(a, 0, 0)/an);

  return 0;

}


/*----------------------------------------------*/

/* Bypass MatGet/MatSet if SML_1DIM is defined */

#if defined(SML_1DIM)
#include "smlcompv.c"
#else

/*----------------------------------------------*/
/*  Essentially 1-dim routines that use MatGet/MatSet -> */
/*----------------------------------------------*/

void MatCopy (MATRIX a, MATRIX b)      /* a = b  */
  /* deep copy with redimensioning, nothing done if a==b  */
{

  size_t i, j, m = MatRows(b), n = MatCols(b);

  if (a==b) return;
  MatReDim(a, m, n);
  if (MatErr()) return;

 for (i = 0; i != m; i++)
    for (j = 0; j != n; j++)
      MatSet0(a, i, j, MatGet0(b, i, j));

}

/*----------------------------------------------*/

REAL MatMin (MATRIX a)             /*  min element of a */
{
  size_t i, j, m = MatRows(a), n = MatCols(a);
  REAL   s;

  if (m==0 || n==0) return SML_NAN;

  s = MatGet0(a, 0, 0);
  for (i = 0; i != m; i++)
    for (j = 0; j != n; j++)
      if (MatGet0(a, i, j) < s)  s = MatGet0(a, i, j);

  return (s);
}

/*----------------------------------------------*/

REAL MatMax (MATRIX a)             /*  max element of a */
{
  size_t i, j, m = MatRows(a), n = MatCols(a);
  REAL   s;

  if (m==0 || n==0) return SML_NAN;

  s = MatGet0(a, 0, 0);
  for (i = 0; i != m; i++)
    for (j = 0; j != n; j++)
      if (MatGet0(a, i, j) > s)  s = MatGet0(a, i, j);

  return (s);
}

/*----------------------------------------------*/

REAL  MatMaxAbs (MATRIX a)             /*  max |element| of a */
{
  size_t i, j, m = MatRows(a), n = MatCols(a);
  REAL   s;

  if (m==0 || n==0) return SML_NAN;

  s = sml_fabs(MatGet0(a, 0, 0));
  for (i = 0; i != m; i++)
    for (j = 0; j != n; j++) {
      REAL   x;
      if ( (x=sml_fabs(MatGet0(a, i, j))) > s)  s = x;
    }

  return s;
}

/*----------------------------------------------*/

REAL     MatMinAbs (MATRIX a)             /*  min |element| of a */
{
  size_t i, j, m = MatRows(a), n = MatCols(a);
  REAL   s;

  if (m==0 || n==0) return SML_NAN;

  s = sml_fabs(MatGet0(a, 0, 0));
  for (i = 0; i != m; i++)
    for (j = 0; j != n; j++) {
      REAL   x;
      if ( (x=sml_fabs(MatGet0(a, i, j))) < s)  s = x;
    }

  return s;
}

/*----------------------------------------------*/

REAL2  MatSum (MATRIX a)             /*  sum of elements of a */
{
  size_t i, j, m = MatRows(a), n = MatCols(a);

  REAL2  s=0;

  if (m==0 || n==0) return 0;

  for (i = 0; i != m; i++)
    for (j = 0; j != n; j++)
      s += MatGet0(a, i, j);

  return s;
}

/*----------------------------------------------*/

REAL2  MatSumAbs (MATRIX a)             /*  sum of |elements| of a */
{
  size_t i, j, m = MatRows(a), n = MatCols(a);

  REAL2  s=0;

  if (m==0 || n==0) return 0;

  for (i = 0; i != m; i++)
    for (j = 0; j != n; j++)
      s += sml_fabs(MatGet0(a, i, j));

  return s;
}

/*----------------------------------------------*/

REAL2  MatSS2  (MATRIX r, MATRIX a)
/*
Returns sum of squares of a
r <- a 2x1 matrix that stores a scaled result,
   sum of squares of a = r[1,1] * r[1,1] * r[2,1]
r == a is OK

*/

{
  size_t i, j, rows, cols;
  REAL2  s=0, result, scale;

  rows=MatRows(a);
  cols=MatCols(a);

  scale = (rows == 0 || cols == 0) ? 0 : MatMaxAbs(a);

  if (scale > 0)  /* we add to s */
    for (i=0; i != rows; i++)
      for (j=0; j != cols; j++) {
        REAL2  x = MatGet0(a, i, j) / scale;
        s += x*x;
      }

  MatReDim(r, 2, 1);
  if (!MatErr()) {
    MatSet(r, 1, 1,  scale);
    MatSet(r, 2, 1,  s);
  }

  result = scale * scale * s;
  return result;
}
/*----------------------------------------------*/

REAL2  MatSS  (MATRIX a)             /*  sum of squares of a */
{
  REAL2  result;
  MATRIX r=MatDim(2,1);

  if (MatErr()) return SML_NAN;
  result = MatSS2 (r, a);
  MatUnDim(r);
  return  result;
}
/*----------------------------------------------*/
REAL2  MatCSum     (MATRIX y, MATRIX x)
/*
Returns sum of x
y <- matrix of cumulative sum of x, by row
y == x is OK
*/
{
  size_t i, j, rows, cols;
  REAL2  s=0;

  rows=MatRows(x);
  cols=MatCols(x);
  MatReDim(y, rows, cols);
  if (MatErr()) return SML_NAN;

  for (i = 0; i != rows; i++)
    for (j = 0; j != cols; j++)
      MatSet0(y, i, j, (s += MatGet0(x, i, j)));

  return s;
}
/*----------------------------------------------*/
REAL2   MatCMean    (MATRIX y, MATRIX x)
/*
Returns mean of x
y <- matrix of cumulative mean of x, by row
y == x is OK
*/
{
  size_t i, j, rows, cols;
  REAL2   s=0, k=0;

  rows=MatRows(x);
  cols=MatCols(x);
  MatReDim(y, rows, cols);
  if (MatErr()) return SML_NAN;

  for (i = 0; i != rows; ++i)
    for (j = 0; j != cols; ++j)
      MatSet0(y, i, j, (s += MatGet0(x, i, j))/(++k));

  return MatGet(y, rows, cols);
}
/*----------------------------------------------*/
void MatFill  (MATRIX a, REAL x)
   /* assign the value x to every element of a */
{
  size_t i, j, rows=MatRows(a), cols=MatCols(a);

  for (j=0; j != cols; ++j)
    for (i=0; i != rows; ++i)
      MatSet0(a, i, j, x);
}

/*----------------------------------------------*/

void MatExpE (MATRIX a, MATRIX b)      /* a = exp(b) elementwise  */
{
  size_t i, j, rows=MatRows(b), cols=MatCols(b);

  MatReDim(a, rows, cols);
  if (MatErr()) return;

  for (j=0; j != cols; ++j)
    for (i=0; i != rows; ++i)
      MatSet0(a, i, j, sml_exp(MatGet0(b, i, j)));
}

/*----------------------------------------------*/

void MatLogE (MATRIX a, MATRIX b)      /* a = log(b) elementwise  */
{
  size_t i, j, rows=MatRows(b), cols=MatCols(b);

  MatReDim(a, rows, cols);
  if (MatErr()) return;

  for (j=0; j != cols; ++j)
    for (i=0; i != rows; ++i)
      MatSet0(a, i, j, sml_log(MatGet0(b, i, j)));
}

/*----------------------------------------------*/

void MatSqrtE (MATRIX a, MATRIX b)      /* a = b^0.5 elementwise  */

{
  size_t i, j, rows=MatRows(b), cols=MatCols(b);

  MatReDim(a, rows, cols);
  if (MatErr()) return;

  for (j=0; j != cols; ++j)
    for (i=0; i != rows; ++i)
      MatSet0(a, i, j, sml_sqrt(MatGet0(b, i, j)));

}

/*----------------------------------------------*/

void MatRecipE (MATRIX a, MATRIX b)      /* a = 1/b elementwise  */
  /* 1/0 gives inf  */

{
  size_t i, j, rows=MatRows(b), cols=MatCols(b);

  MatReDim(a, rows, cols);
  if (MatErr()) return;

  for (j=0; j != cols; ++j)
    for (i=0; i != rows; ++i)
      MatSet0(a, i, j, 1.0/(MatGet0(b, i, j)));

}
/*----------------------------------------------*/

void MatRecipE0 (MATRIX a, MATRIX b)      /* a = 1/b elementwise  */
  /* anything/0 := 0 */

{
  size_t i, j, rows=MatRows(b), cols=MatCols(b);

  MatReDim(a, rows, cols);
  if (MatErr()) return;

  for (j=0; j != cols; ++j)
    for (i=0; i != rows; ++i) {
      REAL x = MatGet0(b, i, j);
      if (x != 0)
        MatSet0(a, i, j, 1.0/x);
      else
        MatSet0(a, i, j, 0);
    }
}

/*----------------------------------------------*/

void MatApplyE (MATRIX a, REAL (*f)(REAL), MATRIX b)

/* a=f(b) elementwise */
{
  size_t i, j, rows=MatRows(b), cols=MatCols(b);

  MatReDim(a, rows, cols);
  if (MatErr()) return;

  for (j=0; j != cols; ++j)
    for (i=0; i != rows; ++i)
      MatSet0(a, i, j, f(MatGet0(b, i, j)));
}
/*----------------------------------------------*/

void MatMulE (MATRIX a, MATRIX b, MATRIX c)  /* a = b * c elementwise */
{
  size_t i, j, rows=MatRows(b), cols=MatCols(b);

  MatReDim(a, rows, cols);
  if (MatErr()) return;

  for (j=0; j != cols; ++j)
    for (i=0; i != rows; ++i)
      MatSet0(a, i, j, MatGet0(b, i, j) * MatGet0(c, i, j));
}
/*----------------------------------------------*/
void   MatDivE (MATRIX a, MATRIX b, MATRIX c)  /* a = b / c   el-wise */
{
  size_t i, j, rows=MatRows(b), cols=MatCols(b);

  MatReDim(a, rows, cols);
  if (MatErr()) return;

  for (j=0; j != cols; ++j)
    for (i=0; i != rows; ++i)
      MatSet0(a, i, j, MatGet0(b, i, j) / MatGet0(c, i, j));
}
/*----------------------------------------------*/
void   MatDivE0 (MATRIX a, MATRIX b, MATRIX c)  /* a = b / c   el-wise */
{
  size_t i, j, rows=MatRows(b), cols=MatCols(b);

  MatReDim(a, rows, cols);
  if (MatErr()) return;

  for (j=0; j != cols; ++j)
    for (i=0; i != rows; ++i) {
      REAL x = MatGet0(c, i, j);
      if (x != 0)
        MatSet0(a, i, j, MatGet0(b, i, j) / x);
      else
        MatSet0(a, i, j, 0);
    }
}
/*----------------------------------------------*/

static void MatAddPE (MATRIX a, MATRIX b)      /* a = a + b */
/*
Called when: a != b
Assumes conformability
Special handling of (a == b)
*/
{
  size_t i, j, m=MatRows(a), n=MatCols(a);

  if (a != b) {
    for (j = 0; j != n; ++j)
      for (i = 0; i != m; ++i)
        MatSetPE0(a, i, j, MatGet0(b, i, j));
  }
  else {
    for (j = 0; j != n; ++j)
      for (i = 0; i != m; ++i)
        MatSetPE0(a, i, j, MatGet0(a, i, j));
  }
}
/*----------------------------------------------*/

static void MatAdd2 (MATRIX a, MATRIX b, MATRIX c)      /* a = b + c */
/*
Called when: a != b && a != c
Assumes conformability
Special handling of (b == c)
*/
{
  size_t i, j, m=MatRows(a), n=MatCols(a);

  if (b != c) {
    for (j = 0; j != n; ++j)
      for (i = 0; i != m; ++i)
        MatSet0(a, i, j, MatGet0(b, i, j) + MatGet0(c, i, j));
  }
  else {
    for (j = 0; j != n; ++j) {
      for (i = 0; i != m; ++i) {
        REAL x = MatGet0(b, i, j);
        MatSet0(a, i, j, x + x);
      }
    }
  }
}
/*----------------------------------------------*/

void MatAdd (MATRIX a, MATRIX b, MATRIX c)      /* a = b + c */

{
  size_t i, j;

  if (MatCols(b) !=  MatCols(c) || MatRows(b) !=  MatRows(c) )   {
    MatErrThrow("MatAdd: unconformable matrices.");
    return;
  }

  MatReDim(a, MatRows(b), MatCols(b));
  if (MatErr()) return;

  if (a == b)      MatAddPE(a, c);      /* a = a + c */
  else if (a == c) MatAddPE(a, b);      /* a = b + a */
  else MatAdd2(a, b, c);                /* a = b + c */
}
/*----------------------------------------------*/

static void MatSubME (MATRIX a, MATRIX b)      /* a = a - b */
/*
Called when: a != b
Assumes conformability
*/
{
  size_t i, j, m=MatRows(a), n=MatCols(a);

    for (j = 0; j != n; ++j)
      for (i = 0; i != m; ++i)
        MatSetME0(a, i, j, MatGet0(b, i, j));
}
/*----------------------------------------------*/

static void MatSub2 (MATRIX a, MATRIX b, MATRIX c)      /* a = b - c */
/*
Called when: a != b && b != c
Assumes conformability
Special handling of (a == c)
*/
{
  size_t i, j, m=MatRows(a), n=MatCols(a);


  if (a != c) {  /* a = b - c */
    for (j = 0; j != n; ++j)
      for (i = 0; i != m; ++i)
        MatSet0(a, i, j, MatGet0(b, i, j) - MatGet0(c, i, j));
  }
  else {         /* a = b - a */
    for (j = 0; j != n; ++j) {
      for (i = 0; i != m; ++i) {
        MatSet0(a, i, j, MatGet0(b, i, j) - MatGet0(a, i, j));
      }
    }
  }
}
/*----------------------------------------------*/

void MatSub (MATRIX a, MATRIX b, MATRIX c)      /* a = b - c */

{
  size_t i, j;

  if (MatCols(b) !=  MatCols(c) || MatRows(b) !=  MatRows(c) )   {
    MatErrThrow("MatSub: unconformable matrices.");
    return;
  }

  MatReDim(a, MatRows(b), MatCols(b));
  if (MatErr()) return;

  if (a == b) MatSubME(a, c);      /* a = a - c */
  else MatSub2(a, b, c);           /* a = b - c */
}
/*----------------------------------------------*/

void MatAddScalar (MATRIX a, MATRIX b, REAL c)      /* a = b + c */
{
  size_t i, j;

  MatReDim(a, MatRows(b), MatCols(b));
  if (MatErr()) return;

  for (i = 0; i != MatRows(a); ++i)
    for (j = 0; j != MatCols(a); ++j)
      MatSet0(a, i, j, MatGet0(b, i, j) + c);
}

/*----------------------------------------------*/

void MatMulScalar (MATRIX a, MATRIX b, REAL c)      /* a = b * c */

{
  size_t i, j;

  MatReDim(a, MatRows(b), MatCols(b));
  if (MatErr()) return;

  for (i = 0; i != MatRows(a); ++i)
    for (j = 0; j != MatCols(a); ++j)
      MatSet0(a, i, j, MatGet0(b, i, j) * c);
}

/*----------------------------------------------*/
void MatSort (MATRIX A)
/*
Sort A, result depends on layout,
Most useful if A is a (row or column) vector
*/
{
  size_t m = MatRows(A);
  size_t n = MatCols(A);
  size_t r = m * n;
  size_t h;

  if (r <= 1) return;
  MatReDim(A, r, 1);  /* access as a column vector */

  /* Shell sort */
  for (h = 1; h <= (r-1)/9; h = 3*h+1) ;
  for ( ; h > 0; h /= 3) {
    size_t i;
    for (i = h; i < r; ++i)   {
        size_t j = i; REAL v = MatGet0(A, i, 0);
        while (j >= h &&   v < MatGet0(A, j-h, 0))
            { MatSet0(A, j, 0, MatGet0(A, j-h, 0)); j -= h; }
        MatSet0(A, j, 0, v);
    }
  }

  MatReDim(A, m, n);  /* restore old dims */
}
/*----------------------------------------------*/
/*----------------------------------------------*/
/* <- Essentially 1-dim routines that use MatGet/MatSet  */
/*----------------------------------------------*/

#endif
