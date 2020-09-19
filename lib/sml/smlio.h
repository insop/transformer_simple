#ifndef _SMLIO_H_
#define _SMLIO_H_

/* smlio.h: Small Matrix library  input/output  header file  */
/***********************************************************/

#include <stdio.h>
#include "smlbase.h"

FILE *open_file (const char *fname, const char *mode);
int  close_file (FILE *f);

size_t  MatRead      (FILE* f,   MATRIX a);
size_t  MatReadFile  (const char fname[], MATRIX a);
size_t  MatWrite     (FILE* f,   MATRIX a, const char *fmt, const char *label);
size_t  MatWriteFile (const char fname[], MATRIX a, const char *fmt, const char *label);

#endif
