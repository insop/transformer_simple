/* smlmem1.c: small matrix library memory management  */
/*            This version is for a->m declared as REAL*    */
/************************************************/
#include <stdlib.h>
#include "smlbase.h"
/************************************************/

#define matrows(a)         ((a)->rows)
#define matcols(a)         ((a)->cols)

/*----------------------------------------------*/
/* definitions */
/*----------------------------------------------*/
REAL* MatData (MATRIX x)

/*
Returns a pointer to matrix data, i.e. first element in the array
*/
{
  return x->m;
}

/*----------------------------------------------*/
/*
This function is here just for symmetry with another version that
interfaces to TNT (Template Numerical Toolkit) and NR (Numerical Recipes
- the book)
*/

REAL* MatDataNR1 (MATRIX x) {return x->m;}

/* REAL** MatDataNR2 (MATRIX x) : not available in this mode! */

/*----------------------------------------------*/
REAL* MatDataij0(MATRIX x, size_t i, size_t j)

/*
returns the address of element (i,j) in the matrix x,
0-based indexing, i=0:nrow-1, j=0:ncol-1
*/

{
  if (i >= MatRows(x) || j >= MatCols(x)) {
    MatErrThrow ("MatDataij0: out of range.");
    return NULL;
  }

#ifdef SML_BY_ROW
  return  (x->m) + (i)*matcols(x) + (j);
#endif

#ifdef SML_BY_COL
  return  (x->m) + (j)*matrows(x) + (i);
#endif
}
/*----------------------------------------------*/
REAL* MatDataij (MATRIX x, size_t i, size_t j)

/*
returns the address of element (i,j) in the matrix x.
*/

{ if (i == 0 || j == 0) return NULL;
  return  MatDataij0(x, i-1, j-1);
}

/*----------------------------------------------*/

MATRIX MatDim (size_t rows, size_t cols)
  /*
  Allocate memory for a matrix. Returns NULL if allocation fails.

  If the return value "a" is not NULL, we have a guarantee that either:
    a->m == NULL && a->rows == a->cols == 0
  or
    a->m is a valid pointer
         && a->rows > 0 && a->cols > 0.
  */
{
  MATRIX a;
  REAL*  data;

  a = (MATRIX) malloc((size_t) sizeof(*a));
  if (a == NULL) {
    MatErrThrow ("MatDim: insufficient memory.");
    return NULL;
  }

  /* handle empty matrices first */
  if (rows == 0 || cols == 0) {
    a->m = NULL;
    matrows(a) = matcols(a) = 0;
    return a;
  }

  /* Now we have: rows > 0 &&  cols > 0 */

  data = (REAL*) malloc((size_t) (rows*cols*sizeof(*data)));
  if (data == NULL)  {
      MatErrThrow ("MatDim: insufficient memory.");
      free(a);
      return  NULL;
  }
  a->m = data;
  matrows(a) = rows;
  matcols(a) = cols;
  return a;
}
/*----------------------------------------------*/

void MatReDim (MATRIX a, size_t r, size_t c)
/*
The behaviour of MatReDim() is crucial to many sml functions.
Operates as a "reshape". As much of the old contents as possible
is preserved.
Original contents are preserved if redimensioning fails.
*/
{
  if (a == NULL) return;

  if (r == 0 || c == 0) r = c = 0;

  /* Same dimensions, including both empty */
  if (matrows(a) == r && matcols(a) == c) return;

  /* Existing matrix is empty, or r*c==0 (but not both)*/
  if (matrows(a) == 0 || matcols(a)==0 || r == 0) {
    MATRIX tmp = MatDim(r, c);
    if (tmp == NULL) return;
    MatSwap(a, tmp);
    MatUnDim(tmp);
    return;
  }

  /* Now, Existing matrix is not empty, r*c > 0  */

  if (r*c != matrows(a)*matcols(a)) {
    REAL* tm =
      (REAL*) realloc(a->m, (size_t) (r*c*sizeof(*tm)));
    if (tm == NULL) {
      MatErrThrow ("MatReDim: insufficient memory.");
      return;
    }
    else {
      a->m = tm;
    }
  }

  matrows(a) = r;
  matcols(a) = c;
  return;
}
/*----------------------------------------------*/

void MatUnDim (MATRIX a)
/* release matrix a, including matrix spec. */
{
  if (a != NULL) {
    free(a->m);
    a->m = NULL; /* obliterate for safety */
    free (a);
  }
}
/*----------------------------------------------*/
void   MatSwap  (MATRIX A, MATRIX B)
/* swap matrices without copying their array data */
/* swaps only the matrix descriptors */
{
  MATRIXDESC t = *A; *A = *B; *B = t;
}
/*----------------------------------------------*/
IVECTOR IvecDim(size_t n)   /* allocates int vec[0:n], i.e. size n+1  */
{
  IVECTOR ivec;


  ivec = (IVECTOR) malloc((size_t) (n+1)*sizeof(ivec[0]));

  if (ivec == NULL) {
    MatErrThrow ("IvecDim: insufficient memory.");
  }

  return ivec;
}
/*----------------------------------------------*/
void IvecUnDim(IVECTOR ivec)  { free(ivec); }
/*----------------------------------------------*/
