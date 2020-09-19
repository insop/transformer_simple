/* smlmem2.c: small matrix library memory management  */
/*            This version is for a->m declared as REAL**   */
/************************************************/
#include <stdlib.h>
#include "smlbase.h"
/************************************************/

#define matrows(a)         ((a)->rows)
#define matcols(a)         ((a)->cols)

static REAL*  matrix (size_t r, size_t c);
static REAL** ptrs   (size_t dim1);
static void   get_dims(size_t* dim1, size_t* dim2, size_t rows, size_t cols);
static void   setup_ptrs(REAL** m, REAL* data, size_t dim1, size_t dim2);


/*----------------------------------------------*/
REAL* MatData (MATRIX x)

/*
Returns a pointer to matrix data, i.e. first element in the array
*/
{
  return x->m[0];
}

/*----------------------------------------------*/
/*
These two functions allow passing a MATRIX to code from the
TNT (Template Numerical Toolkit) and NR (Numerical Recipes - the book)
*/

REAL*  MatDataNR1 (MATRIX x) {return x->m[0];} /* 1-dim array */
REAL** MatDataNR2 (MATRIX x) {return x->m;}    /* 2-dim array */

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
  return  x->m[i] + j;
#endif

#ifdef SML_BY_COL
  return  x->m[j] + i;
#endif
}

/*----------------------------------------------*/
REAL* MatDataij (MATRIX x, size_t i, size_t j)

/*
returns the address of element (i,j) in the matrix x,
1-based indexing, i=1:nrow, j=1:ncol
*/

{
  if (i==0 || j == 0) return NULL;
  return  MatDataij0(x, i-1, j-1);
}

/*----------------------------------------------*/

MATRIX MatDim (size_t rows, size_t cols)
  /*
  Allocate memory for a matrix. Returns NULL if allocation fails.

  If the return value "a" is not NULL, we have a guarantee that either:
    a->m == NULL && a->rows == a->cols == 0
  or
    a->m and a->m[0:dim1-1] are valid pointers
         && a->rows > 0 && a->cols > 0.
  dim1 := a->rows if storage is by_row
  dim1 := a->cols if storage is by_col
  */
{
  MATRIX a;
  REAL** m;
  REAL*  data;
  size_t   dim1, dim2;

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

  get_dims(&dim1, &dim2, rows, cols);

  /* Now we have: dim1 > 0 and dim2 > 0 */

  /* pointers to dim1 */
  m = ptrs(dim1);
  if (m == NULL) {
      /* ptrs() raised the error flag, free space and return */
      free ((void *) a);
      return NULL;
  }

  /* actual data */
  data = matrix(dim1, dim2);

  if (data == NULL) {
      /* matrix() raised the error flag, free space and return */
      free ((void *) m);
      free ((void *) a);
      return  NULL;
  }

  setup_ptrs(m, data, dim1, dim2);
  a->m = m;
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
  size_t   dim1old, dim2old, dim1new, dim2new;

  if (a == NULL) return;

  if (r == 0 || c == 0) r = c = 0;

  /* Same dimensions, including both empty */
  if (matrows(a) == r && matcols(a) == c) return;

  /* Existing matrix is empty, or r*c==0 (but not both empty)*/
  if (matrows(a) == 0 || matcols(a)==0 || r == 0) {
    MATRIX tmp = MatDim(r, c);
    if (tmp == NULL) return;
    MatSwap(a, tmp);
    MatUnDim(tmp);
    return;
  }

  /* Now, Existing matrix is not empty, r*c > 0  */

  get_dims(&dim1old, &dim2old, matrows(a), matcols(a));
  get_dims(&dim1new, &dim2new, r, c);

  if (dim1new != dim1old) {
    REAL** tm =
      (REAL**) realloc(a->m, (size_t) (dim1new*sizeof(*tm)));
    if (tm == NULL) {
      MatErrThrow ("MatReDim: insufficient memory.");
      return;
    }
    else {
      a->m = tm;
    }
  }

  if (r*c != matrows(a)*matcols(a)) {
    REAL* td =
        (REAL*) realloc(a->m[0], (size_t) (r*c*sizeof(*td)));
    if (td == NULL) {
      MatErrThrow ("MatReDim: insufficient memory.");
      return;
    }
    else {
      a->m[0] = td;
    }
  }

  setup_ptrs(a->m, a->m[0], dim1new, dim2new);
  matrows(a) = r;
  matcols(a) = c;
  return;

}
/*----------------------------------------------*/

void MatUnDim (MATRIX a)
        /* free all space bound  to matrix a  */
{
  if (a != NULL) {
    if (a->m != NULL) {
      free (a->m[0]);
      free (a->m);
      a->m = NULL; /* obliterate for safety */
    }
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
IVECTOR IvecDim(size_t n)   /* allocates iv[0:n], i.e. size n+1  */
{
  IVECTOR ivec;


  ivec = (IVECTOR) malloc((size_t) (n+1)*sizeof(ivec[0]));

  if (ivec == NULL) {
    MatErrThrow ("IvecDim: insufficient memory.");
  }

  return ivec;
}
/*----------------------------------------------*/
void IvecUnDim(IVECTOR ivec)  { free((void*) ivec); }
/*----------------------------------------------*/
/* static code */
/*----------------------------------------------*/
static void  setup_ptrs(REAL** m, REAL* data, size_t dim1, size_t dim2)
/* m[] <- pointers to data[] */
{
  size_t  i;

  for (i = 0; i != dim1; ++i) {m[i] = data;  data += dim2;}
}
/*----------------------------------------------*/

static REAL* matrix (size_t r, size_t c)
/*
Allocates memory for  an array of r*c elements.
Returns NULL if r*c==0 or malloc() error.
*/

{
  REAL* data =
         (REAL*) malloc((size_t) (r*c*sizeof(*data)));

  if (data == NULL)
    MatErrThrow ("matrix: insufficient memory.");

  return data;
}
/*----------------------------------------------*/
static REAL** ptrs  (size_t dim1)
  /* allocate space for pointers to dim1 (rows or columns) */
{
  REAL** m =
      (REAL**) malloc((size_t) (dim1*sizeof(*m)));

  if (m == NULL)
    MatErrThrow ("ptrs: insufficient memory.");

  return m;
}
/*----------------------------------------------*/

static void get_dims(size_t* dim1, size_t* dim2, size_t rows, size_t cols)
{
#ifdef SML_BY_ROW
  *dim1 = rows; *dim2 = cols;
#endif

#ifdef SML_BY_COL
  *dim1 = cols; *dim2 = rows;
#endif
}
