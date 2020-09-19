/* smlchol.c:  small matrix library Cholesky decomposition  */
/**********************************************************/

#include <float.h>
#include <math.h>

#include "smlcomp.h"

/************************************************************/

/*----------------------------------------------*/
int  MatChol (MATRIX L, MATRIX A)
/*
   Cholesky decomposition of a symmetric positive semidefinite
   matrix A = L L'.
   See MatCholIP for details.
   A == L is OK
*/
{
  int  rc;

  MatCopy(L, A);
  if (MatErr()) return 0;
  rc = MatCholIP(L);
  MatFillAD(L, 0.0);
  return rc;
}
/*----------------------------------------------*/
int  MatCholIP (MATRIX a)
/*
   In-Place Cholesky decomposition of a symmetric positive semidefinite
   matrix a = L L'.
   L is returned in the lower triangle of a.
   Symmetry of a is not checked (upper a is not accessed).
   Error, returns (0) if a is not square (nothing done).
   Returns (-i) if the leading i by i submatrix is not
     semidefinite positive
   Returns rank(a) otherwise.
*/
{
    int     rank;
    size_t    i, n;
    const REAL  eps = REAL_EPSILON;  /* relative 0 */

  if  (!MatIsSquare(a))  {
     MatErrThrow("MatCholIP: nonsquare matrix.");
     return 0;
  }

  rank = n = MatRows(a);
  for (i = 0; i != n; i++) {    /* i-loop: rows */
    REAL  aii, delta;
    REAL2 s = 0;
    size_t j;
    for (j = 0; j != i; j++) {   /* j-loop: columns, off diagonals */
      if (MatGet0(a, j, j) == 0)  MatSet0(a, i, j, 0);
      else {
        REAL2 t = 0;
        size_t k;
        for (k = 0; k != j; k++)
          t +=  MatGet0(a, j, k) * MatGet0(a, i, k);   /* inner product */
        MatSet0(a, i, j,  (MatGet0(a, i, j) - t) / MatGet0(a, j, j));
      }
    } /* end of j-loop: off diagonal */
    /* the diagonal element */
    for (j = 0; j != i; j++)
      s +=  MatGet0(a, i, j)*MatGet0(a, i, j); /* sum of squares */
    delta = eps*MatGet0(a, i, i);
    delta = (delta >= 0) ? delta : -delta;
    aii   = MatGet0(a, i, i) - s;
    if (aii < -delta)   return -i;   /* -ve def or indef */
    else if (aii > delta)  MatSet0(a, i, i, sml_sqrt(aii)); /* no problem */
    else  {
      MatSet0(a, i, i, 0);
      rank--;
    }
  }    /* end of i-loop: rows */
  return rank;
}
/*----------------------------------------------*/
/* begin  generalized inverse based on Cholesky */
/*----------------------------------------------*/
int  MatInvLIP (MATRIX L)
/*
  g-invert a square lower triangular matrix L in situ.
  Upper L not accessed.
  Returns 0 if successful.
  Returns 1 if L is not square (nothing done to L).
*/
{
  size_t i, n;

  if  (!MatIsSquare(L))   return 1;
  n = MatRows(L);
  for (i = 0; i != n; i++) {
    REAL w = MatGet0(L, i, i);
    if (w != 0.0) {
      size_t j;
      for (j = 0; j != i; j++) {
        REAL2 s = 0;
        size_t k;
        for (k = j; k != i; k++)
          s +=  MatGet0(L, i, k) * MatGet0(L, k, j);
        MatSet0(L ,i, j,  -s / w);
      }
      MatSet0(L, i, i,  1.0 / w);
    }
    else {  /* w == 0.0 */
      size_t j;
      for (j = 0; j != i; j++)    MatSet0(L , i, j,  0.0);
    }
  }
  return 0;
}
/*----------------------------------------------*/
int  MatLtLIP (MATRIX L)
/*
  L <-  L'L,  where L is square lower triangular,
  Returns 0 if successful.
  Returns 1 if L is not square (nothing done to L).
*/
{
  size_t  i, j, k, rows;

  rows = MatRows(L);
  if  (rows != MatCols(L)) return 1;

  for (i = 0; i != rows; i++) {
    for (j = 0; j != i; j++) {
      REAL2 s = 0;
      for (k = i; k != rows; k++)
        s += MatGet0(L, k, i)*MatGet0(L, k, j);
      MatSet0(L, i, j,  s);  /*  lower half */
      MatSet0(L, j, i,  s);  /*  upper half */
    }
  }
  return 0;
}
/*----------------------------------------------*/
int  MatLLtIP (MATRIX L)
/*
  L <-  L L',  where L is square lower triangular,
  Returns 0 if successful.
  Returns 1 if L is not square (nothing done to L).
*/
{
  int  i, j, k, rows;

  rows = MatRows(L);
  if  (rows != (int) MatCols(L)) return 1;

  for (i = rows-1; i >= 0; i--) {
    for (j = i; j >= 0; j--) {
      REAL2 s = 0;
      for (k = j; k >= 0; k--)
        s += MatGet0(L, i, k)*MatGet0(L, j, k);
      MatSet0(L, i, j,  s);  /*  lower half */
      MatSet0(L, j, i,  s);  /*  upper half */
    }
  }
  return 0;
}
/*----------------------------------------------*/
int  MatSHInvIP (MATRIX A)
/*
lower A <- g-inverse of the Cholesky root of A
returns the return code from MatCholIP()
*/
{
  int rc;

  rc = MatCholIP(A) ;  /* lower Cholesky root in-situ */
  if (rc >  0)  MatInvLIP (A) ;  /* invert in-situ */
  return rc;
}
/*----------------------------------------------*/
int  MatSInvIP (MATRIX A)
/*
g-inverts a symmetric non-negative definite matrix A in situ
returns the return code from MatChol()
*/
{
  int rc;

  rc = MatSHInvIP (A);
  if (rc > 0)  MatLtLIP (A);  /* A <- A'A */
  return rc;
}
/*----------------------------------------------*/
/* end    generalized inverse based on Cholesky */
/*----------------------------------------------*/

int MatCholSolve(MATRIX X, MATRIX L, MATRIX B)

/*

  Solve a linear system A*X = B, using the previously computed
  cholesky factorization of A: L*L'.

   @param  B   A Matrix with as many rows as A and any number of columns.
   @computes   X so that L*L'*X = B.
*/
{
  size_t n  = MatRows(L);
  size_t nx = MatCols(B);
  size_t i, j;
  int  k;

  if (X == L)            return 1;
  if (n != MatRows(B))   return 2;

  MatCopy(X, B);
  if (MatErr()) return 1;

      // Solve L*X = B;
      for (j=0; j != nx; j++)
        for (k = 0; k != n; k++) {
          for (i = 0; i != k; i++)
               MatSetME0(X,k,j, MatGet0(X,i,j)*MatGet0(L,k,i));
          MatSetDE0(X,k,j, MatGet0(L,k,k));
        }

      // Solve L'*X = Y;
     for (j=0; j != nx; j++)
        for (k = n-1; k >= 0; k--) {
          for (i = k+1; i != n; i++)
               MatSetME0(X,k,j, MatGet0(X,i,j)*MatGet0(L,i,k));
          MatSetDE0(X,k,j,  MatGet0(L,k,k));
        }

  return 0;
}


