/* smlio.c:  Small Matrix Library input/output  */
/************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "smlio.h"
/************************************************/
FILE *open_file (const char *fname, const char *mode)
{
  FILE *f;
  errno = 0;
  f = fopen(fname, mode);
  if (f == NULL)
    fprintf (stderr, "Open of %s with mode=%s failed: %s\n",
             fname, mode, strerror(errno));
  return (f);
}
/************************************************/
int  close_file (FILE *f)
{
  int s;
  if (f == NULL) return (0);
  errno = 0;
  s = fclose(f);
  if (s == EOF) perror("Close failed");
  return s;
}
/************************************************/
size_t   MatRead (FILE* f, MATRIX a)
/*
read matrix a from file f (by row)
returns # elements read
*/
{
  size_t  numread=0, rows, cols, i, j;

  rows = MatRows(a);
  cols = MatCols(a);
  for (i=1; i <= rows; i++)
    for (j=1; j <= cols; j++)  {
      REAL x;
      int  rc = fscanf (f, SML_SCNREAL, &x);
      if (rc == 1)  {
        MatSet(a, i, j, x); numread++;
      }
      else {
        if (rc == EOF) MatErrThrow("MatRead: end of file error.");
        else           MatErrThrow("MatRead: input error.");
        return  numread;
      }
    }

  return  numread;

} /* MatRead */
/************************************************/
size_t    MatReadFile (const char fname[], MATRIX a)
/* returns # elements read */
{
  FILE*  InFile;
  size_t n;

  InFile = open_file (fname, "r");
  if (InFile == NULL)  {
    MatErrThrow("MatReadFile: IO error.");
    return 0;
  }

  n = MatRead (InFile, a);
  close_file (InFile);
  return  n;
}
/************************************************/
size_t  MatWrite  (FILE* f, MATRIX a, const char* fmt, const char label[])
  /* print matrix a to file f, with rows labelled */
{
  const  size_t CPP = 12; /* max columns per page */
  const  char s[2][2]  = {"", "s"};
  char   ifmt[10] = "[%4d] ";
  size_t   i, j, column, rows, cols, nwritten=0, rc;

  fprintf (f, "%s\n", label);
  rows = MatRows(a);
  cols = MatCols(a);
  fprintf (f, "%d row%s, %d column%s\n",
      rows, s[rows>1], cols, s[cols>1]);
  for (i = 1; i <= rows; i++) {
    column = 0;
    fprintf (f, ifmt, i);
    for (j = 1; j <= cols; j++) {
      if (++column > CPP) {      /* end of line? */
        column = 1;           /* new line */
        fprintf (f, "\n");
      }
      rc = fprintf (f, fmt, MatGet(a, i, j));
      if (rc != EOF) nwritten++;
      else {
        MatErrThrow("MatWrite: IO error.");
        return (nwritten);
      }
    }
    fprintf (f, "\n");  /* end of row */
  } /* all rows done */

  fprintf (f, "\n");
  return nwritten;

} /* MatWrite */
/************************************************/
size_t  MatWriteFile (const char fname[], MATRIX a, const char *fmt, const char *label)
/*
write matrix to file
returns # written
*/
{
  FILE*  OutFile;
  size_t nwritten;

  OutFile = open_file (fname, "w");
  if (OutFile == NULL)  {
    MatErrThrow("MatWriteFile: IO error.");
    return 0;
  }

  nwritten = MatWrite (OutFile, a, fmt, label);
  close_file (OutFile);
  return nwritten;
}
/************************************************/
