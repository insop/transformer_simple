#ifndef _SMLBASE_H_
#define _SMLBASE_H_

/* smlbase2.h: Small Matrix Library configuration and basic structures */

/*
Adds support for long double.

Requires C99 for "float".

Some non-C99 compilers support long double versions of the math
functions ending with "l" (ell), e.g. sqrtl, logl, etc.  On such
compilers, "long double" should work fine.

Doesn't use <tgmath.h>, so that "double" works on non-C99 compilers.

*/
/***********************************************************/

/* size_t  */
#include <stddef.h>

/* numeric constants */
#include <float.h>

/* SML  Configuration */
#include "sml.cfg"

/**********************************************************************/
/* SML - check the defs */

#if defined(SML_BY_ROW) && defined(SML_BY_COL)
  #error "SML_BY_ROW and SML_BY_COL defined, choose one."
#endif

#if !defined(SML_BY_ROW) && !defined(SML_BY_COL)
  #error "define either SML_BY_ROW or SML_BY_COL"
#endif

#if defined(SML_REAL_DBL) && defined(SML_REAL_FLT)
    #error "SML_REAL_DBL and SML_REAL_FLT defined, choose one"
#endif

#if defined(SML_REAL_LDBL) && defined(SML_REAL_FLT)
    #error "SML_REAL_LDBL and SML_REAL_FLT defined, choose one"
#endif

#if defined(SML_REAL_LDBL) && defined(SML_REAL_DBL)
    #error "SML_REAL_LDBL and SML_REAL_DBL defined, choose one"
#endif

#if !defined(SML_REAL_FLT) && !defined(SML_REAL_DBL) && !defined(SML_REAL_LDBL)
  #error "define SML_REAL_FLT or SML_REAL_DBL or SML_REAL_LDBL"
#endif

#if defined(SML_AC_FUN) && defined(SML_AC_MAC)
  #error "SML_AC_FUN and SML_AC_MAC defined, choose one."
#endif

#if !defined(SML_AC_FUN) && !defined(SML_AC_MAC)
  #error "define either SML_AC_FUN or SML_AC_MAC"
#endif
/**********************************************************************/

#ifdef SML_REAL_FLT
  typedef  float REAL;
#endif
#ifdef SML_REAL_DBL
  typedef  double REAL;
#endif
#ifdef SML_REAL_LDBL
  typedef  long double REAL;
#endif


typedef struct  {
  size_t rows, cols;     /* dimensions */
  REAL** m;              /* pointers to rows or columns */
 }  MATRIXDESC;          /* matrix descriptor structure */

typedef MATRIXDESC* MATRIX;

typedef int*        IVECTOR;

/**********************************************************************/
/* other types,  format strings  */

#ifdef SML_REAL_FLT
  #define REAL2        double
  #define REAL_EPSILON FLT_EPSILON
  #define REAL_MIN     FLT_MIN
  #define REAL_MAX     FLT_MAX
  #define SML_SCNREAL  "%f"
  #define SML_PRTREAL  ""
  #define SML_PRTREAL2 ""
#endif

#ifdef SML_REAL_DBL
  #define REAL2        long double
  #define REAL_EPSILON DBL_EPSILON
  #define REAL_MIN     DBL_MIN
  #define REAL_MAX     DBL_MAX
  #define SML_SCNREAL  "%lf"
  #define SML_PRTREAL  ""
  #define SML_PRTREAL2 "L"
#endif

#ifdef SML_REAL_LDBL
  #define REAL2        long double
  #define REAL_EPSILON LDBL_EPSILON
  #define REAL_MIN     LDBL_MIN
  #define REAL_MAX     LDBL_MAX
  #define SML_SCNREAL  "%Lf"
  #define SML_PRTREAL  "L"
  #define SML_PRTREAL2 "L"
#endif

/**********************************************************************/

/* some math functions  */

#if defined(__STDC_VERSION__) && __STDC_VERSION__>=199901L

      /* C99 */
#ifdef SML_REAL_FLT
  #define sml_fabs(x)  fabsf(x)
  #define sml_sqrt(x)  sqrtf(x)
  #define sml_log(x)   logf(x)
  #define sml_exp(x)   expf(x)
#endif

#ifdef SML_REAL_DBL
  #define sml_fabs(x)  fabs(x)
  #define sml_sqrt(x)  sqrt(x)
  #define sml_log(x)   log(x)
  #define sml_exp(x)   exp(x)
#endif

#ifdef SML_REAL_LDBL
  #define sml_fabs(x)  fabsl(x)
  #define sml_sqrt(x)  sqrtl(x)
  #define sml_log(x)   logl(x)
  #define sml_exp(x)   expl(x)
#endif

#else
      /* not C99 */
  #define sml_fabs(x)  fabs(x)
  #define sml_sqrt(x)  sqrt(x)
  #define sml_log(x)   log(x)
  #define sml_exp(x)   exp(x)

#endif

/**********************************************************************/

/* macro access to matrix contents */
#ifdef  SML_AC_MAC
#include "smlacm.h"
#endif

/*  function access to matrix contents */
#ifdef  SML_AC_FUN
#include "smlacf.h"
#endif

/*  memory management */
#include "smlmem.h"

/*  error handling */
#include "smlerr.h"

#endif
