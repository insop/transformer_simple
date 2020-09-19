#ifndef _SMLERR_H_
#define _SMLERR_H_

/* smlerr.h: Small Matrix Library Error Handling header file  */
/***********************************************************/

/* functions that set or clear values return old values */
int  MatSetAbort (void); /* abort on error */
int  MatClrAbort (void); /* do not abort on error */
int  MatClrErr (void);   /* clear error flag */
int  MatErr (void);      /* return status */
int  MatErrCount (void); /* return # errors */
void MatErrThrow (const char* msg);

#endif
