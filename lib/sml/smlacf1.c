/* smlacf1.c: function access to matrix contents, +- bounds checking */
/*            This version is for a->m declared as REAL*    */
/**********************************************************/

#ifdef SML_AC_FUN

#include <stdio.h>
#include "smlbase.h"

#define SML_NAN  ((REAL)0.0/(REAL)0.0)

#define matrows(a)         ((a)->rows)
#define matcols(a)         ((a)->cols)
#define out(a,i,j)   ((i)>matrows(a) || (j)>matcols(a) || (i)<1 || (j)<1)
#define out0(a,i,j)  ((i)>=matrows(a)|| (j)>=matcols(a))

#ifdef SML_BY_ROW
  #define matel0(a,i,j)      ((a)->m[(i)*matcols(a)+(j)])
#endif

#ifdef SML_BY_COL
  #define matel0(a,i,j)      ((a)->m[(j)*matrows(a)+(i)])
#endif

#define matel(a,i,j)       matel0(a,i-1,j-1)

/***********************************************************/
size_t    MatRows(MATRIX a) { return    matrows(a); }
/***********************************************************/
size_t    MatCols(MATRIX a) { return    matcols(a); }
/***********************************************************/
REAL    MatGet(MATRIX a, size_t i, size_t j)
{

#ifdef SML_BOUNDS_CHK
  if (out(a, i, j)) {
    printf("MatGet: i=%u, j=%u, rows=%u, cols=%u\n", i, j, matrows(a), matcols(a));
    MatErrThrow("MatGet: out of bounds access.");
    return  SML_NAN;
  }
#endif

  return matel(a, i, j);
}
/***********************************************************/

REAL    MatSet(MATRIX a, size_t i, size_t j, REAL x)
{

#ifdef SML_BOUNDS_CHK
  if (out(a, i, j)) {
    printf("MatSet: i=%u, j=%u, rows=%u, cols=%u\n", i, j, matrows(a), matcols(a));
    MatErrThrow("MatSet: out of bounds access.");
    return  SML_NAN;
  }
#endif

  return matel (a, i, j) = x;
}

/***********************************************************/
REAL    MatSetPE(MATRIX a, size_t i, size_t j, REAL x)
{
#ifdef SML_BOUNDS_CHK
  if (out(a, i, j)) {
    printf("MatSetPE: i=%u, j=%u, rows=%u, cols=%u\n", i, j, matrows(a), matcols(a));
    MatErrThrow("MatSetPE: out of bounds access.");
    return  SML_NAN;
  }
#endif

  return matel(a, i, j) += x;
}
/***********************************************************/
REAL    MatSetME(MATRIX a, size_t i, size_t j, REAL x)
{
#ifdef SML_BOUNDS_CHK
  if (out(a, i, j)) {
    printf("MatSetME: i=%u, j=%u, rows=%u, cols=%u\n", i, j, matrows(a), matcols(a));
    MatErrThrow("MatSetME: out of bounds access.");
    return  SML_NAN;
  }
#endif

  return matel(a, i, j) -= x;
}
/***********************************************************/
REAL    MatSetTE(MATRIX a, size_t i, size_t j, REAL x)
{
#ifdef SML_BOUNDS_CHK
  if (out(a, i, j)) {
    printf("MatSetTE: i=%u, j=%u, rows=%u, cols=%u\n", i, j, matrows(a), matcols(a));
    MatErrThrow("MatSetTE: out of bounds access.");
    return  SML_NAN;
  }
#endif

  return matel(a, i, j) *= x;
}
/***********************************************************/
REAL    MatSetDE(MATRIX a, size_t i, size_t j, REAL x)
{
#ifdef SML_BOUNDS_CHK
  if (out(a, i, j)) {
    printf("MatSetDE: i=%u, j=%u, rows=%u, cols=%u\n", i, j, matrows(a), matcols(a));
    MatErrThrow("MatSetDE: out of bounds access.");
    return  SML_NAN;
  }
#endif

  return matel(a, i, j) /= x;
}
/***********************************************************/
/* versions for 0-based indexing */
/***********************************************************/
REAL    MatGet0(MATRIX a, size_t i, size_t j)
{

#ifdef SML_BOUNDS_CHK
  if (out0(a, i, j)) {
    printf("MatGet0: i=%u, j=%u, rows=%u, cols=%u\n", i, j, matrows(a), matcols(a));
    MatErrThrow("MatGet0: out of bounds access.");
    return  SML_NAN;
  }
#endif

  return matel0(a, i, j);
}
/***********************************************************/

REAL    MatSet0(MATRIX a, size_t i, size_t j, REAL x)
{

#ifdef SML_BOUNDS_CHK
  if (out0(a, i, j)) {
    printf("MatSet0: i=%u, j=%u, rows=%u, cols=%u\n", i, j, matrows(a), matcols(a));
    MatErrThrow("MatSet0: out of bounds access.");
    return  SML_NAN;
  }
#endif

  return matel0 (a, i, j) = x;
}

/***********************************************************/
REAL    MatSetPE0(MATRIX a, size_t i, size_t j, REAL x)
{
#ifdef SML_BOUNDS_CHK
  if (out0(a, i, j)) {
    printf("MatSetPE0: i=%u, j=%u, rows=%u, cols=%u\n", i, j, matrows(a), matcols(a));
    MatErrThrow("MatSetPE0: out of bounds access.");
    return  SML_NAN;
  }
#endif

  return matel0(a, i, j) += x;
}
/***********************************************************/
REAL    MatSetME0(MATRIX a, size_t i, size_t j, REAL x)
{
#ifdef SML_BOUNDS_CHK
  if (out0(a, i, j)) {
    printf("MatSetME0: i=%u, j=%u, rows=%u, cols=%u\n", i, j, matrows(a), matcols(a));
    MatErrThrow("MatSetME0: out of bounds access.");
    return  SML_NAN;
  }
#endif

  return matel0(a, i, j) -= x;
}
/***********************************************************/
REAL    MatSetTE0(MATRIX a, size_t i, size_t j, REAL x)
{
#ifdef SML_BOUNDS_CHK
  if (out0(a, i, j)) {
    printf("MatSetTE0: i=%u, j=%u, rows=%u, cols=%u\n", i, j, matrows(a), matcols(a));
    MatErrThrow("MatSetTE0: out of bounds access.");
    return  SML_NAN;
  }
#endif

  return matel0(a, i, j) *= x;
}
/***********************************************************/
REAL    MatSetDE0(MATRIX a, size_t i, size_t j, REAL x)
{
#ifdef SML_BOUNDS_CHK
  if (out0(a, i, j)) {
    printf("MatSetDE0: i=%u, j=%u, rows=%u, cols=%u\n", i, j, matrows(a), matcols(a));
    MatErrThrow("MatSetDE0: out of bounds access.");
    return  SML_NAN;
  }
#endif

  return matel0(a, i, j) /= x;
}
/***********************************************************/
#endif
