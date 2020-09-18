/* mario.c:  Matrix array (MAR) library input/output  */
/************************************************/
#include <stdlib.h>
#include "mar.h"
/************************************************/

void MarRead (FILE* f, MAR a)
  /* read matrix array a from file f (by row) */
  /* no i/o error checking  */
{
  size_t  i;
  for (i = 1; i <= MarRows(a); i++)    MatRead (f, MarEl(a, i)) ;
} /* MarRead */

/*----------------------------------------------*/
void   MarReadFile (const char* fname, MAR a)
{
  FILE*  InFile;

  InFile = open_file (fname, "r");
  if (InFile == NULL)  exit (2);

  MarRead (InFile, a);

  close_file (InFile);
}
/*----------------------------------------------*/

void MarWrite  (FILE* f, MAR a,const  char *fmt,const  char *title)
  /* print matrix array a to file f */
  /* no i/o error checking  */
{
  size_t  i;
  char  s[20];

  fprintf (f, "%s\n", title);

  for (i = 1; i <= MarRows(a); i++)  {
    sprintf  (s, "matrix %d:", i);
    MatWrite (f, MarEl(a, i), fmt, s) ;
  }

} /* MarWrite */

/*----------------------------------------------*/

void MarWriteFile (const char* fname, MAR a, const char *fmt,const char* title)
{
  FILE*  OutFile;

  OutFile = open_file (fname, "w");
  if (OutFile == NULL)  exit (2);

  MarWrite (OutFile, a, fmt, title);

  close_file (OutFile);
}
