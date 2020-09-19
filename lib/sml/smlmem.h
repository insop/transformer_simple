#ifndef _SMLMEM_H_
#define _SMLMEM_H_

/* smlmem.h: Small Matrix Library memory management header file  */
/***********************************************************/

MATRIX MatDim   (size_t rows, size_t cols);
void   MatReDim (MATRIX a, size_t rows, size_t cols);
void   MatUnDim (MATRIX);
void   MatSwap  (MATRIX, MATRIX);

REAL*  MatData   (MATRIX x);                 /* 1-dim array */
REAL*  MatDataij (MATRIX x, size_t i, size_t j); /* address of x[i,j] */
REAL*  MatDataij0(MATRIX x, size_t i, size_t j); /* address of x[i,j], 0-based */

REAL*  MatDataNR1(MATRIX x);    /* 1-dim array for NR and TNT */
REAL** MatDataNR2(MATRIX x);    /* 2-dim array for NR and TNT */

IVECTOR IvecDim(size_t n);  /* allocates iv[0:n], i.e. size n+1  */
void    IvecUnDim(IVECTOR ivec);


#endif
