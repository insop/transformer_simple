/* smlacf2.c: Function access to matrix contents, +- bounds checking */
/*            This version is for a->m declared as REAL**   */
/**********************************************************/

#ifdef SML_AC_FUN

#include "smlbase.h"

#define SML_NAN  ((REAL)0.0/(REAL)0.0)

#define matrows(a)         ((a)->rows)
#define matcols(a)         ((a)->cols)
#define out(a,i,j)   ((i)>matrows(a) || (j)>matcols(a) || (i)<1 || (j)<1)
#define out0(a,i,j)  ((i)>=matrows(a)|| (j)>=matcols(a))

#ifdef SML_BY_ROW
  #define matel0(a,i,j)      ((a)->m[i][j])
#endif

#ifdef SML_BY_COL
  #define matel0(a,i,j)      ((a)->m[j][i])
#endif

#define matel(a,i,j)       matel0((a),(i)-1,(j)-1)

/***********************************************************/
size_t    MatRows(MATRIX a) { return    matrows(a); }
/***********************************************************/
size_t    MatCols(MATRIX a) { return    matcols(a); }
/***********************************************************/
REAL    MatGet(MATRIX a, size_t i, size_t j)
{

#ifdef SML_BOUNDS_CHK
  if (out(a, i, j)) {
    printf("i=%u, j=%u, rows=%u, cols=%u\n", i, j, matrows(a), matcols(a));
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
    printf("i=%u, j=%u, rows=%u, cols=%u\n", i, j, matrows(a), matcols(a));
    MatErrThrow("MatSet: out of bounds access.");
    return SML_NAN;
  }
#endif

  return matel(a, i, j) = x;
}

/***********************************************************/
REAL    MatSetPE(MATRIX a, size_t i, size_t j, REAL x)
{
#ifdef SML_BOUNDS_CHK
  if (out(a, i, j)) {
    MatErrThrow("MatSetPE: out of bounds access.");
    return SML_NAN;
  }
#endif

  return matel(a, i, j) += x;
}
/***********************************************************/
REAL    MatSetME(MATRIX a, size_t i, size_t j, REAL x)
{
#ifdef SML_BOUNDS_CHK
  if (out(a, i, j)) {
    MatErrThrow("MatSetME: out of bounds access.");
    return SML_NAN;
  }
#endif

  return matel(a, i, j) -= x;
}
/***********************************************************/
REAL    MatSetTE(MATRIX a, size_t i, size_t j, REAL x)
{
#ifdef SML_BOUNDS_CHK
  if (out(a, i, j)) {
    MatErrThrow("MatSetTE: out of bounds access.");
    return SML_NAN;
  }
#endif

  return matel(a, i, j) *= x;
}
/***********************************************************/
REAL    MatSetDE(MATRIX a, size_t i, size_t j, REAL x)
{
#ifdef SML_BOUNDS_CHK
  if (out(a, i, j)) {
    MatErrThrow("MatSetDE: out of bounds access.");
    return SML_NAN;
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
    MatErrThrow("MatSet0: out of bounds access.");
    return SML_NAN;
  }
#endif

  return matel0(a, i, j) = x;
}

/***********************************************************/
REAL    MatSetPE0(MATRIX a, size_t i, size_t j, REAL x)
{
#ifdef SML_BOUNDS_CHK
  if (out0(a, i, j)) {
    MatErrThrow("MatSetPE0: out of bounds access.");
    return SML_NAN;
  }
#endif

  return matel0(a, i, j) += x;
}
/***********************************************************/
REAL    MatSetME0(MATRIX a, size_t i, size_t j, REAL x)
{
#ifdef SML_BOUNDS_CHK
  if (out0(a, i, j)) {
    MatErrThrow("MatSetME0: out of bounds access.");
    return SML_NAN;
  }
#endif

  return matel0(a, i, j) -= x;
}
/***********************************************************/
REAL    MatSetTE0(MATRIX a, size_t i, size_t j, REAL x)
{
#ifdef SML_BOUNDS_CHK
  if (out0(a, i, j)) {
    MatErrThrow("MatSetTE0: out of bounds access.");
    return SML_NAN;
  }
#endif

  return matel0(a, i, j) *= x;
}
/***********************************************************/
REAL    MatSetDE0(MATRIX a, size_t i, size_t j, REAL x)
{
#ifdef SML_BOUNDS_CHK
  if (out0(a, i, j)) {
    MatErrThrow("MatSetDE0: out of bounds access.");
    return SML_NAN;
  }
#endif

  return matel0(a, i, j) /= x;
}
/***********************************************************/
#endif
