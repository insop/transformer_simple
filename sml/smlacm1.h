#ifndef _SMLACM_H_
#define _SMLACM_H_

/* smlacm.h: Small Matrix Library access by macros */
/***********************************************************/

#define MatRows(a)         ((a)->rows)
#define MatCols(a)         ((a)->cols)

/* 0-based indexing, append "0" to macro or function name */
#ifdef  SML_BY_ROW
#define MatEl0(a,i,j)      ((a)->m[(i)*MatCols(a)+(j)])
#endif

#ifdef  SML_BY_COL
#define MatEl0(a,i,j)      ((a)->m[(j)*MatRows(a)+(i)])
#endif

#define MatEl(a,i,j)        MatEl0(a,(i)-1,(j)-1)

#define MatGet(a,i,j)       MatEl(a,i,j)
#define MatSet(a,i,j,x)     MatEl(a,i,j)=(x)
#define MatSetPE(a,i,j,x)   MatEl(a,i,j)+=(x)
#define MatSetME(a,i,j,x)   MatEl(a,i,j)-=(x)
#define MatSetTE(a,i,j,x)   MatEl(a,i,j)*=(x)
#define MatSetDE(a,i,j,x)   MatEl(a,i,j)/=(x)

#define MatGet0(a,i,j)      MatEl0(a,i,j)
#define MatSet0(a,i,j,x)    MatEl0(a,i,j)=(x)
#define MatSetPE0(a,i,j,x)  MatEl0(a,i,j)+=(x)
#define MatSetME0(a,i,j,x)  MatEl0(a,i,j)-=(x)
#define MatSetTE0(a,i,j,x)  MatEl0(a,i,j)*=(x)
#define MatSetDE0(a,i,j,x)  MatEl0(a,i,j)/=(x)

#endif
