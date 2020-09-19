/* marmem.c: Matrix array (MAR) library memory management  */
/**********************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mar.h"
/************************************************/

static MATRIX* mar (size_t  n);
static void marundim2 (MATRIX* a, size_t i, size_t j);

/*----------------------------------------------*/
/* definitions */
/*----------------------------------------------*/

MAR MarDim (size_t rows)
  /*
  Otherwise if rows=0 memory will be allocated only for the
  matrix descriptor, but not the data array, no error condition raised,
  and rows will be set to 0
  */
{
  MAR a;

  a = (MAR) malloc((size_t) sizeof(*a));
  if (a == NULL) {
    MatErrThrow ("MarDim: insufficient memory.");
    return (MAR) NULL;
  }

  a->m = mar(rows);  /* mar(0)  ok */

  if (MatErr()) {
      /* mar() raised the error flag, free space and return */
      free ((void *) a);
      return (MAR) NULL;
  }

  MarRows(a) = rows;

  return a;
}

/*----------------------------------------------*/

void   MarReDimAll (MAR a, MATRIX rows, size_t cols)
{
  size_t i;
  if (MarRows(a) != MatRows(rows)) {
    MatErrThrow ("MarReDimAll: invalid dims.");
    return;
  }

  for (i=1; i <= MatRows(rows); i++)
    MatReDim (MarEl(a, i), (size_t) MatGet(rows, i, 1), cols);
}
/*----------------------------------------------*/

void   MarReDim1   (MAR a, size_t i, size_t rows, size_t cols)
  /* Redimension one matrix */
{
   if (i < 1 || i > MarRows(a))  {
     MatErrThrow ("MarReDim1: invalid arguments.");
     return;
   }

  MatReDim (MarEl(a, i), rows, cols);
}

/*----------------------------------------------*/

void   MarReDim    (MAR a, size_t rows)
  /* reallocate if different size */
  /* follows the MarDim() behaviour */
{


  if (MarRows(a) != rows)   { /* different size */
    free ((void *) (a->m)); /* array */
    a->m = mar(rows);
  }

  if (MatErr())  MarRows(a) = 0; else MarRows(a) = rows;

} /* MarReDim */
/*----------------------------------------------*/

void MarUnDim (MAR a)
        /* release a, including all matrices */
{
  if (a == NULL) return;

  if (MarRows(a) > 0) {
    marundim2 (a->m, 0, MarRows(a)-1);   /* matrices */
    free ((void *) (a->m)); /* array */
  }

  free ((void *) a);    /* mar descriptor */

} /* MarUnDim */
/*----------------------------------------------*/
/* static code */
/*----------------------------------------------*/

static MATRIX* mar (size_t  n)
/*
  Allocate memory for  array[0:n-1] of matrices,
  each matrix initialized with 0 rows and 0 columns.
  Either all or none allocation, no allocated memory left hanging.
*/

{
  size_t i;
  MATRIX* a;


  if (n == 0) return NULL;

  a = (MATRIX*) malloc((size_t) (n*sizeof(MATRIX)) );
  if (a == NULL) MatErrThrow ("mar: insufficient memory.");
  else {
    for (i=0; i < n; i++) {
      a[i] = MatDim(0, 0);
      if (MatErr()) break;     /* a[i] failed */
    }
    if (MatErr()) {            /* all or none */
      if (i > 0) marundim2(a, 0, i-1);       /* free matrices a[0:i-1] */
      free ((void*) a);
      a =  (MATRIX*) NULL;
      MatErrThrow ("mar: insufficient memory.");
    }
  }

  return  a;
}

/*----------------------------------------------*/

static void marundim2 (MATRIX* a, size_t i, size_t j)
  /* undim matrices a[i:j] */
{

  if (i > j) return;

  /* i==0 possible, so i <= j is always true because size_t */
  while (i < j)     MatUnDim(a[j--]);

  MatUnDim(a[i]);

}
/*----------------------------------------------*/
