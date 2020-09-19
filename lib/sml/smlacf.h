#ifndef _SMLACF_H_
#define _SMLACF_H_

/* smlacf.h: Small Matrix Library access by functions */
/***********************************************************/

size_t    MatRows(MATRIX a);
size_t    MatCols(MATRIX a);
REAL    MatGet(MATRIX a, size_t i, size_t j);
REAL    MatSet(MATRIX a, size_t i, size_t j, REAL x);
REAL    MatSetPE(MATRIX a, size_t i, size_t j, REAL x);
REAL    MatSetME(MATRIX a, size_t i, size_t j, REAL x);
REAL    MatSetTE(MATRIX a, size_t i, size_t j, REAL x);
REAL    MatSetDE(MATRIX a, size_t i, size_t j, REAL x);

REAL    MatGet0(MATRIX a, size_t i, size_t j);
REAL    MatSet0(MATRIX a, size_t i, size_t j, REAL x);
REAL    MatSetPE0(MATRIX a, size_t i, size_t j, REAL x);
REAL    MatSetME0(MATRIX a, size_t i, size_t j, REAL x);
REAL    MatSetTE0(MATRIX a, size_t i, size_t j, REAL x);
REAL    MatSetDE0(MATRIX a, size_t i, size_t j, REAL x);

#endif
