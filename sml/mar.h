#ifndef _MAR_H_
#define _MAR_H_

/* mar.h: Matrix Array library header file  */
/***********************************************************/
/* MAR - include files */

#include "smlbase.h"
#include "smlio.h"

/* MAR  - general defs  */

#define MarRows(a) ((a)->rows)
#define MarCols(a) ((a)->rows)
#define MarEl(a,i) ((a)->m[(i)-1])

typedef struct
 {
  size_t     rows;          /* dimension  */
  MATRIX*  m;             /* pointer to actual contents [0:rows-1] */
 }  MARDESC;         /* matrix array descriptor */

typedef MARDESC* MAR;

/* MAR  - memory management */

MAR    MarDim   (size_t n);    /* allocate an array of n matrices */
void   MarReDim    (MAR a, size_t n);
void   MarReDim1   (MAR a, size_t i, size_t rows, size_t cols);
void   MarReDimAll (MAR a, MATRIX rows, size_t cols);
void   MarUnDim (MAR a);

/* MAR  - input/output  */

void MarRead (FILE* f, MAR a);
void MarReadFile (const char* fname, MAR a);
void MarWrite  (FILE* f, MAR a, const char *fmt, const char *title);
void MarWriteFile (const char* fname, MAR a, const char *fmt, const char* title);

#endif
