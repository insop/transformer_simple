/* smlerr.c  - small matrix library error handling */
/*----------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include "smlbase.h"
/*----------------------------------------------*/

#define TRUE  1
#define FALSE 0

static int MatAbortFlag = TRUE; /* abort on error ? */
static int MatErrFlag   = FALSE; /* status of last Mat call */
static int MatErrN      = 0;  /* cumumaltive # errors */

/*----------------------------------------------*/

int MatSetAbort (void)  /* abort on error */
{
  int tmp = MatAbortFlag;
  MatAbortFlag = TRUE;
  return tmp;
}
/*----------------------------------------------*/

int MatClrAbort (void)  /* do not abort on error */
{
  int tmp = MatAbortFlag;
  MatAbortFlag = FALSE;
  return tmp;
}
/*----------------------------------------------*/

int MatClrErr (void)  /* clear error flag */
{
  int tmp = MatErrFlag;
  MatErrFlag = FALSE;
  MatErrN    = 0;
  return tmp;
}
/*----------------------------------------------*/

int MatErr (void)  /* return status */
{   return MatErrFlag; }
/*----------------------------------------------*/

int MatErrCount (void)  /* return # errors */
{   return MatErrN; }
/*----------------------------------------------*/

void MatErrThrow (const char msg[])
{
  fprintf (stderr, "SML error: %s\n", msg);
  MatErrN++;
  MatErrFlag = TRUE;
  if (MatAbortFlag) exit(1);
}
/*----------------------------------------------*/
