/* Computer generated file. Don't edit. */

/* eg.c: demonstration of some SML facilities */
/****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sml.h"

int main ()
{
  MATRIX A, B, C, D, E, H, I, J, L, Q, R, S, T, U, V, W, X, Y, Z;
  size_t   n, i, j, k;
  long double x0, x1, x2, x3, x4, x5, x6, x7, x8, x9;

  const char* fmt = "%8" SML_PRTREAL "g ";

  printf("\n");

  printf("  /* initialize */\n");
  /* initialize */
  printf("  A = MatDim(0, 0);\n");
  A = MatDim(0, 0);
  printf("  B = MatDim(0, 0);\n");
  B = MatDim(0, 0);
  printf("  C = MatDim(0, 0);\n");
  C = MatDim(0, 0);
  printf("  D = MatDim(0, 0);\n");
  D = MatDim(0, 0);
  printf("  E = MatDim(0, 0);\n");
  E = MatDim(0, 0);
  printf("  H = MatDim(0, 0);\n");
  H = MatDim(0, 0);
  printf("  I = MatDim(0, 0);\n");
  I = MatDim(0, 0);
  printf("  J = MatDim(0, 0);\n");
  J = MatDim(0, 0);
  printf("  L = MatDim(0, 0);\n");
  L = MatDim(0, 0);
  printf("  Q = MatDim(0, 0);\n");
  Q = MatDim(0, 0);
  printf("  R = MatDim(0, 0);\n");
  R = MatDim(0, 0);
  printf("  S = MatDim(0, 0);\n");
  S = MatDim(0, 0);
  printf("  T = MatDim(0, 0);\n");
  T = MatDim(0, 0);
  printf("  U = MatDim(0, 0);\n");
  U = MatDim(0, 0);
  printf("  V = MatDim(0, 0);\n");
  V = MatDim(0, 0);
  printf("  W = MatDim(0, 0);\n");
  W = MatDim(0, 0);
  printf("  X = MatDim(0, 0);\n");
  X = MatDim(0, 0);
  printf("  Y = MatDim(0, 0);\n");
  Y = MatDim(0, 0);
  printf("  Z = MatDim(0 ,0);\n");
  Z = MatDim(0 ,0);
  printf("\n");

  printf("/*******************************************************/\n");
/*******************************************************/
  printf("\n");

  printf("  /* changing dimensions */\n");
  /* changing dimensions */
  printf("\n");

  printf("  MatReDim(X, 0, 0);\n");
  MatReDim(X, 0, 0);
  printf("  MatReDim(X, 0, 0);\n");
  MatReDim(X, 0, 0);
  printf("  MatReDim(X, 5, 5);\n");
  MatReDim(X, 5, 5);
  printf("  MatReDim(X, 2, 3);\n");
  MatReDim(X, 2, 3);
  printf("  MatReDim(X, 1, 6);\n");
  MatReDim(X, 1, 6);
  printf("  MatReDim(X, 6, 1);\n");
  MatReDim(X, 6, 1);
  printf("  MatReDim(X, 4, 4);\n");
  MatReDim(X, 4, 4);
  printf("  MatReDim(X, 4, 4);\n");
  MatReDim(X, 4, 4);
  printf("\n");

  printf("/*******************************************************/\n");
/*******************************************************/
  printf("\n");

  printf("  /* scalar functions */\n");
  /* scalar functions */
  printf("  n = 5;\n");
  n = 5;
  printf("  MatJ(X, n, n, 1);\n");
  MatJ(X, n, n, 1);
  printf("\n");

  printf("  x1 = MatMin(X);\n");
  x1 = MatMin(X);
  printf("  printf(\"%%Lg\\n\", x1);\n");
  printf("%Lg\n", x1);
  printf("\n");

  printf("  x1 = MatMax(X);\n");
  x1 = MatMax(X);
  printf("  printf(\"%%Lg\\n\", x1);\n");
  printf("%Lg\n", x1);
  printf("\n");

  printf("  x1 = MatMinAbs(X);\n");
  x1 = MatMinAbs(X);
  printf("  printf(\"%%Lg\\n\", x1);\n");
  printf("%Lg\n", x1);
  printf("\n");

  printf("  x1 = MatMaxAbs(X);\n");
  x1 = MatMaxAbs(X);
  printf("  printf(\"%%Lg\\n\", x1);\n");
  printf("%Lg\n", x1);
  printf("\n");

  printf("  x1 = MatSum(X);\n");
  x1 = MatSum(X);
  printf("  printf(\"%%Lg\\n\", x1);\n");
  printf("%Lg\n", x1);
  printf("\n");

  printf("  x1 = MatSumAbs(X);\n");
  x1 = MatSumAbs(X);
  printf("  printf(\"%%Lg\\n\", x1);\n");
  printf("%Lg\n", x1);
  printf("\n");

  printf("  x1 = MatSS(X);\n");
  x1 = MatSS(X);
  printf("  printf(\"%%Lg\\n\", x1);\n");
  printf("%Lg\n", x1);
  printf("\n");

  printf("/*******************************************************/\n");
/*******************************************************/
  printf("\n");

  printf(" /* swap */\n");
 /* swap */
  printf("\n");

  printf("  MatI(X, 2, -1);\n");
  MatI(X, 2, -1);
  printf("  MatI(Y, 4,  1);\n");
  MatI(Y, 4,  1);
  MatWrite (stdout , X, fmt,"X:");
  MatWrite (stdout , Y, fmt,"Y:");
  printf("  MatSwap(X, Y);\n");
  MatSwap(X, Y);
  MatWrite (stdout , X, fmt,"X:");
  MatWrite (stdout , Y, fmt,"Y:");
  printf("\n");

  printf("/*******************************************************/\n");
/*******************************************************/
  printf("\n");

  printf("  /*  sums and means */\n");
  /*  sums and means */
  printf("\n");

  printf("  MatReDim(X, 5, 5);\n");
  MatReDim(X, 5, 5);
  printf("  MatFill(X, -1.4);\n");
  MatFill(X, -1.4);
  MatWrite (stdout , X, fmt,"X:");
  printf("\n");

  printf("  x1 = MatSum(X);\n");
  x1 = MatSum(X);
  printf("  x2 = MatSumAbs(X);\n");
  x2 = MatSumAbs(X);
  printf("  x3 = MatSS(X);\n");
  x3 = MatSS(X);
  printf("  printf(\"%%Lg %%Lg %%Lg \\n\", x1, x2, x3);\n");
  printf("%Lg %Lg %Lg \n", x1, x2, x3);
  printf("\n");

  printf("\n");

  printf("  x1 = MatSS(X);\n");
  x1 = MatSS(X);
  printf("  x2 = MatSS2(R, X);\n");
  x2 = MatSS2(R, X);
  printf("  x3 = MatGet(R, 1, 1)*MatGet(R, 1, 1)*MatGet(R, 2, 1);\n");
  x3 = MatGet(R, 1, 1)*MatGet(R, 1, 1)*MatGet(R, 2, 1);
  printf("  printf(\"%%Lg %%Lg %%Lg \\n\", x1, x2, x3);\n");
  printf("%Lg %Lg %Lg \n", x1, x2, x3);
  printf("\n");

  printf("  MatColSum(Y, X);\n");
  MatColSum(Y, X);
  MatWrite (stdout , Y, fmt,"Y:");
  printf("\n");

  printf("  MatRowSum(Y, X);\n");
  MatRowSum(Y, X);
  MatWrite (stdout , Y, fmt,"Y:");
  printf("\n");

  printf("  MatJ(X, 4, 4, 1);\n");
  MatJ(X, 4, 4, 1);
  MatWrite (stdout , X, fmt,"X:");
  printf("\n");

  printf("  /* The \"cumulative\" functions reveal the layout (by-row or by-column)*/\n");
  /* The "cumulative" functions reveal the layout (by-row or by-column)*/
  printf("  /* If a matrix consists of 1 row or 1 column, layout makes no difference */\n");
  /* If a matrix consists of 1 row or 1 column, layout makes no difference */
  printf("  x1 = MatCSum (Y, X);  /* Y <- cumulative sum  of X */\n");
  x1 = MatCSum (Y, X);  /* Y <- cumulative sum  of X */
  MatWrite (stdout , Y, fmt,"Y:");
  printf("\n");

  printf("  x2 = MatCMean(Y, Y);  /* Y <- cumulative mean of Y */\n");
  x2 = MatCMean(Y, Y);  /* Y <- cumulative mean of Y */
  MatWrite (stdout , Y, fmt,"Y:");
  printf("\n");

  printf("  printf(\"%%Lg %%Lg \\n\", x1, x2);\n");
  printf("%Lg %Lg \n", x1, x2);
  printf("\n");

  printf("/*******************************************************/\n");
/*******************************************************/
  printf("\n");

  printf("  /* filling */\n");
  /* filling */
  printf("\n");

  printf("  MatReDim(X, 7, 3);\n");
  MatReDim(X, 7, 3);
  printf("  MatFill(X, 0);\n");
  MatFill(X, 0);
  printf("  MatFillRow(X, 3, 1); /* Fill a row */\n");
  MatFillRow(X, 3, 1); /* Fill a row */
  MatWrite (stdout , X, fmt,"X:");
  printf("  MatFillCol(X, 2,-1); /* Fill a column */\n");
  MatFillCol(X, 2,-1); /* Fill a column */
  MatWrite (stdout , X, fmt,"X:");
  printf("\n");

  printf("  MatFill(X, 0);\n");
  MatFill(X, 0);
  printf("  MatFillAD(X, 1);     /* Fill Above the Diagonal */\n");
  MatFillAD(X, 1);     /* Fill Above the Diagonal */
  MatWrite (stdout , X, fmt,"X:");
  printf("\n");

  printf("\n");

  printf("  MatFillBD(X, 2);     /* Fill Below the Diagonal */\n");
  MatFillBD(X, 2);     /* Fill Below the Diagonal */
  MatWrite (stdout , X, fmt,"X:");
  printf("\n");

  printf("\n");

  printf("  MatFillD(X, 3);      /* Fill in the Diagonal */\n");
  MatFillD(X, 3);      /* Fill in the Diagonal */
  MatWrite (stdout , X, fmt,"X:");
  printf("\n");

  printf("  MatReDim(X, 3, 7);\n");
  MatReDim(X, 3, 7);
  printf("  MatFill(X, 0);\n");
  MatFill(X, 0);
  printf("  MatFillAD(X, 1);\n");
  MatFillAD(X, 1);
  MatWrite (stdout , X, fmt,"X:");
  printf("\n");

  printf("  MatFillBD(X, 2);\n");
  MatFillBD(X, 2);
  MatWrite (stdout , X, fmt,"X:");
  printf("\n");

  printf("  MatFillD(X, 3);\n");
  MatFillD(X, 3);
  MatWrite (stdout , X, fmt,"X:");
  printf("\n");

  printf("\n");

  printf("  MatBFill(X, 1, 2, 4, 6, 0);     /* Block Fill */\n");
  MatBFill(X, 1, 2, 4, 6, 0);     /* Block Fill */
  MatWrite (stdout , X, fmt,"X:");
  printf("\n");

  printf("  n = 4;\n");
  n = 4;
  printf("  MatI (X, n, 3);          /* diagonal matrix */\n");
  MatI (X, n, 3);          /* diagonal matrix */
  MatWrite (stdout , X, fmt,"X:");
  printf("\n");

  printf("  MatJ (X, n, n+n, 3);\n");
  MatJ (X, n, n+n, 3);
  MatWrite (stdout , X, fmt,"X:");
  printf("\n");

  printf("/*******************************************************/\n");
/*******************************************************/
  printf("\n");

  printf("  /* sorting */\n");
  /* sorting */
  printf("\n");

  printf("  MatReDim(X, 1, 10);\n");
  MatReDim(X, 1, 10);
  printf("  for (i = 1; i <= 10; ++i) MatSet(X, 1, i, i * i %% 7);\n");
  for (i = 1; i <= 10; ++i) MatSet(X, 1, i, i * i % 7);
  MatWrite (stdout , X, fmt,"X:");
  printf("  MatSort(X);\n");
  MatSort(X);
  MatWrite (stdout , X, fmt,"X:");
  printf("/*******************************************************/\n");
/*******************************************************/
  printf("\n");

  printf("  /* elementwise */\n");
  /* elementwise */
  printf("\n");

  printf("  n = 4;\n");
  n = 4;
  printf("  MatJ(A, n, n, 1);\n");
  MatJ(A, n, n, 1);
  printf("  MatJ(B, n, n, 2);\n");
  MatJ(B, n, n, 2);
  printf("  MatJ(C, n, n, 3);\n");
  MatJ(C, n, n, 3);
  printf("\n");

  printf("  MatMulE(A, A, A);\n");
  MatMulE(A, A, A);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatMulE(A, A, C);\n");
  MatMulE(A, A, C);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatMulE(A, B, A);\n");
  MatMulE(A, B, A);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatMulE(A, B, C);\n");
  MatMulE(A, B, C);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("\n");

  printf("  MatDivE(A, A, A);\n");
  MatDivE(A, A, A);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatDivE(A, A, C);\n");
  MatDivE(A, A, C);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatDivE(A, B, A);\n");
  MatDivE(A, B, A);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatDivE(A, B, C);\n");
  MatDivE(A, B, C);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("\n");

  printf("  MatSqrtE(A, A);\n");
  MatSqrtE(A, A);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatSqrtE(A, B);\n");
  MatSqrtE(A, B);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatI(B, n,  exp((double)1.0));\n");
  MatI(B, n,  exp((double)1.0));
  MatWrite (stdout , B, fmt,"B:");
  printf("\n");

  printf("  MatRecipE(A, B);\n");
  MatRecipE(A, B);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatRecipE0(A, B);\n");
  MatRecipE0(A, B);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatLogE(A, B);\n");
  MatLogE(A, B);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatExpE(A, B);\n");
  MatExpE(A, B);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatI(B, n, 2*atan2((double)1, (double)1));\n");
  MatI(B, n, 2*atan2((double)1, (double)1));
  MatWrite (stdout , B, fmt,"B:");
  printf("  MatApplyE(A, sin, B);         /* A = sin(B) elementwise */\n");
  MatApplyE(A, sin, B);         /* A = sin(B) elementwise */
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("/*******************************************************/\n");
/*******************************************************/
  printf("\n");

  printf("  /* add, subtract */\n");
  /* add, subtract */
  printf("\n");

  printf("  n = 4;\n");
  n = 4;
  printf("  MatJ(A, n, n, 1);\n");
  MatJ(A, n, n, 1);
  printf("  MatJ(B, n, n, 2);\n");
  MatJ(B, n, n, 2);
  printf("  MatJ(C, n, n, 3);\n");
  MatJ(C, n, n, 3);
  printf("\n");

  printf("  MatAdd(A, B, C);\n");
  MatAdd(A, B, C);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatAdd(A, A, C);\n");
  MatAdd(A, A, C);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatAdd(A, B, A);\n");
  MatAdd(A, B, A);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatAdd(A, A, A);\n");
  MatAdd(A, A, A);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatSub(A, B, C);\n");
  MatSub(A, B, C);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatSub(A, A, C);\n");
  MatSub(A, A, C);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatSub(A, B, A);\n");
  MatSub(A, B, A);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatSub(A, A, A);\n");
  MatSub(A, A, A);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatSub(A, B, B);\n");
  MatSub(A, B, B);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatI(B, n, -1);\n");
  MatI(B, n, -1);
  printf("  MatLogE(A, B);\n");
  MatLogE(A, B);
  printf("  MatSub(A, A, A);      /* A - A != 0 (there are NAN's) */\n");
  MatSub(A, A, A);      /* A - A != 0 (there are NAN's) */
  MatWrite (stdout , A, fmt,"A:");
  printf("/*******************************************************/\n");
/*******************************************************/
  printf("\n");

  printf("  /* multiply */\n");
  /* multiply */
  printf("\n");

  printf("  n = 4;\n");
  n = 4;
  printf("  MatJ(A, n, n, 1);\n");
  MatJ(A, n, n, 1);
  printf("  MatJ(B, n, n, 2);\n");
  MatJ(B, n, n, 2);
  printf("  MatJ(C, n, n, 3);\n");
  MatJ(C, n, n, 3);
  printf("\n");

  printf("  MatMul(A, B, C);\n");
  MatMul(A, B, C);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatMul(A, A, C);\n");
  MatMul(A, A, C);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatMul(A, B, A);\n");
  MatMul(A, B, A);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("  MatJ(A, n, n, 1);\n");
  MatJ(A, n, n, 1);
  printf("  MatMul(A, A, A);\n");
  MatMul(A, A, A);
  MatWrite (stdout , A, fmt,"A:");
  printf("\n");

  printf("/*******************************************************/\n");
/*******************************************************/
  printf("\n");

  printf("  n = 5;\n");
  n = 5;
  printf("  MatI (I, n, 1);          /* identity matrix */\n");
  MatI (I, n, 1);          /* identity matrix */
  MatWrite (stdout , I, fmt,"I:");
  printf("  MatI (X, n, 1);\n");
  MatI (X, n, 1);
  printf("  for (i = 1; i <= n; i++) MatSet(X, i, n+1-i, 1);\n");
  for (i = 1; i <= n; i++) MatSet(X, i, n+1-i, 1);
  MatWrite (stdout , X, fmt,"X:");
  printf("  MatCopy (S, X);               /* save in S */\n");
  MatCopy (S, X);               /* save in S */
  printf("\n");

  printf("\n");

  printf("/*******************************************************/\n");
/*******************************************************/
  printf("\n");

  printf("  /* Householder  QR */\n");
  /* Householder  QR */
  printf("\n");

  printf("  MatCopy (X, S);               /* restore */\n");
  MatCopy (X, S);               /* restore */
  printf("  /* The unpacked form */\n");
  /* The unpacked form */
  printf("  MatQR (Q, R, X);\n");
  MatQR (Q, R, X);
  MatWrite (stdout , Q, fmt,"Q:");
  MatWrite (stdout , R, fmt,"R:");
  printf("  MatMul(Y, Q, R);     /* Y <- Q * R */\n");
  MatMul(Y, Q, R);     /* Y <- Q * R */
  printf("  MatSub(Y, X, Y);    /*  Y <- X - Y */\n");
  MatSub(Y, X, Y);    /*  Y <- X - Y */
  printf("  x1 = MatSumAbs(Y);\n");
  x1 = MatSumAbs(Y);
  printf("\n");

  printf("  MatWMul (\"t  \", Y, Q, NULL, Q);       /* Y <- Q' * Q  */\n");
  MatWMul ("t  ", Y, Q, NULL, Q);       /* Y <- Q' * Q  */
  MatWrite (stdout , Y, fmt,"Y:");
  printf("  MatSub(Y, Y, I);\n");
  MatSub(Y, Y, I);
  printf("  x2 = MatSumAbs(Y);\n");
  x2 = MatSumAbs(Y);
  printf("  printf(\"MatQR Error: %%Lg   %%Lg\\n\", x1, x2 );\n");
  printf("MatQR Error: %Lg   %Lg\n", x1, x2 );
  printf("\n");

  printf("  /* The packed form and how to unpack */\n");
  /* The packed form and how to unpack */
  printf("  MatQRgetQR(T, D, X);\n");
  MatQRgetQR(T, D, X);
  printf("  MatQRgetQ (Q, T);\n");
  MatQRgetQ (Q, T);
  printf("  MatQRgetR (R, T, D);\n");
  MatQRgetR (R, T, D);
  printf("  MatQRgetH (H, T);\n");
  MatQRgetH (H, T);
  MatWrite (stdout , T, fmt,"T:");
  MatWrite (stdout , Q, fmt,"Q:");
  MatWrite (stdout , R, fmt,"R:");
  MatWrite (stdout , H, fmt,"H:");
  printf("\n");

  printf("/*******************************************************/\n");
/*******************************************************/
  printf("\n");

  printf("  /* Cholesky  */\n");
  /* Cholesky  */
  printf("\n");

  printf("  MatWMul(\"  t\", X, S, NULL, S);       /* X <- S * S'  */\n");
  MatWMul("  t", X, S, NULL, S);       /* X <- S * S'  */
  printf("  MatChol(L, X);                       /* Cholesky: X = L * L' */\n");
  MatChol(L, X);                       /* Cholesky: X = L * L' */
  MatWrite (stdout , L, fmt,"L:");
  printf("  MatCopy   (A, L);\n");
  MatCopy   (A, L);
  printf("  MatLLtIP  (A);                       /* A <- L * L' */\n");
  MatLLtIP  (A);                       /* A <- L * L' */
  MatWrite (stdout , A, fmt,"A:");
  printf("  MatSub(Y, A, X);\n");
  MatSub(Y, A, X);
  printf("  x1 = MatSumAbs(Y);\n");
  x1 = MatSumAbs(Y);
  printf("\n");

  printf("  MatCopy(V, X);\n");
  MatCopy(V, X);
  printf("  MatCholIP(V);                        /* Cholesky in-place */\n");
  MatCholIP(V);                        /* Cholesky in-place */
  MatWrite (stdout , V, fmt,"V:");
  printf("\n");

  printf("  MatTran(R, R);                        /* convert R to lower */\n");
  MatTran(R, R);                        /* convert R to lower */
  printf("  MatMulScalar(R, R, -1.0);             /* change signs */\n");
  MatMulScalar(R, R, -1.0);             /* change signs */
  printf("  MatSub(Y, L, R);\n");
  MatSub(Y, L, R);
  printf("  x2 = MatSumAbs(Y);\n");
  x2 = MatSumAbs(Y);
  printf("  printf(\"MatChol Error: %%Lg  %%Lg\\n\", x1, x2 );\n");
  printf("MatChol Error: %Lg  %Lg\n", x1, x2 );
  printf("\n");

  printf("\n");

  printf("  MatCopy(Z, X);\n");
  MatCopy(Z, X);
  printf("  x1 = MatSInvIP(Z);                           /* g-inverse  */\n");
  x1 = MatSInvIP(Z);                           /* g-inverse  */
  printf("  printf(\"MatXInv rank(Z) = %%d \\n\", (int) x1  );\n");
  printf("MatXInv rank(Z) = %d \n", (int) x1  );
  printf("\n");

  printf("  MatMul (V, X, Z);\n");
  MatMul (V, X, Z);
  printf("  MatMul (V, V, X);\n");
  MatMul (V, V, X);
  MatWrite (stdout , V, fmt,"V:");
  printf("  MatMul (V, Z, X);\n");
  MatMul (V, Z, X);
  printf("  MatMul (V, V, Z);\n");
  MatMul (V, V, Z);
  MatWrite (stdout , V, fmt,"V:");
  printf("\n");

  printf("\n");

  printf("/*******************************************************/\n");
/*******************************************************/
  printf("\n");

  printf("  /* SVD  */\n");
  /* SVD  */
  printf("\n");

  printf("  MatCopy (A, S);               /* we'll use S */\n");
  MatCopy (A, S);               /* we'll use S */
  printf("  MatCopy (X, S);               /* restore */\n");
  MatCopy (X, S);               /* restore */
  MatWrite (stdout , X, fmt,"X:");
  printf("\n");

  printf("  x1 = MatSVD (U, S, V, X);  /*  SVD: X = U * diag(S) * V'   */\n");
  x1 = MatSVD (U, S, V, X);  /*  SVD: X = U * diag(S) * V'   */
  MatWrite (stdout , U, fmt,"U:");
  MatWrite (stdout , S, fmt,"S:");
  MatWrite (stdout , V, fmt,"V:");
  printf("\n");

  printf("  printf(\"Return code = %%d \\n\", (int) x1 );\n");
  printf("Return code = %d \n", (int) x1 );
  printf("\n");

  printf("  MatWMul (\" wt\", Y, U, S, V);   /*  Y <- U * diag(S) * V'   */\n");
  MatWMul (" wt", Y, U, S, V);   /*  Y <- U * diag(S) * V'   */
  MatWrite (stdout , Y, fmt,"Y:");
  printf("\n");

  printf("  MatSub(Y, X, Y);              /* Y <- Y - X   */\n");
  MatSub(Y, X, Y);              /* Y <- Y - X   */
  printf("  x1 = MatSumAbs(Y);\n");
  x1 = MatSumAbs(Y);
  printf("\n");

  printf("  MatWMul (\"t  \", Y, U, NULL, U);       /* Y <- U' * U  */\n");
  MatWMul ("t  ", Y, U, NULL, U);       /* Y <- U' * U  */
  printf("  MatSub(Y, Y, I);\n");
  MatSub(Y, Y, I);
  printf("  x2 = MatSumAbs(Y);\n");
  x2 = MatSumAbs(Y);
  printf("\n");

  printf("  MatWMul (\"t  \", Y, V, NULL, V);       /* Y <- V' * V  */\n");
  MatWMul ("t  ", Y, V, NULL, V);       /* Y <- V' * V  */
  printf("  MatSub(Y, Y, I);\n");
  MatSub(Y, Y, I);
  printf("  x3 = MatSumAbs(Y);\n");
  x3 = MatSumAbs(Y);
  printf("  x4 = sml_fabs(MatSS(X) - MatSS(S));\n");
  x4 = sml_fabs(MatSS(X) - MatSS(S));
  printf("\n");

  printf("  printf(\"MatSVD Error: %%Lg  %%Lg  %%Lg  %%Lg\\n\", x1, x2, x3, x4 );\n");
  printf("MatSVD Error: %Lg  %Lg  %Lg  %Lg\n", x1, x2, x3, x4 );
  printf("\n");

  printf("\n");

  printf("\n");

  printf("  MatRecipE0(S, S);              /*  S <- 1/S  */\n");
  MatRecipE0(S, S);              /*  S <- 1/S  */
  printf("  MatWMul (\" wt\", Y, V, S, U);   /*  Y <- ginv(X) = V * diag(1/S) * S' */\n");
  MatWMul (" wt", Y, V, S, U);   /*  Y <- ginv(X) = V * diag(1/S) * S' */
  MatWrite (stdout , Y, fmt,"Y:");
  printf("  MatMul (Z, X, Y);\n");
  MatMul (Z, X, Y);
  printf("  MatMul (Z, Z, X);\n");
  MatMul (Z, Z, X);
  MatWrite (stdout , Z, fmt,"Z:");
  printf("\n");

  printf("  MatCopy (S, A);               /* restore S */\n");
  MatCopy (S, A);               /* restore S */
  printf("/*******************************************************/\n");
/*******************************************************/
  printf("\n");

  printf("  /* Eigen  */\n");
  /* Eigen  */
  printf("\n");

  printf("  MatCopy (X, S);               /* restore */\n");
  MatCopy (X, S);               /* restore */
  printf("\n");

  printf("  MatWMul (\"t  \", S, X, NULL, X);       /* s <- x' * x  */\n");
  MatWMul ("t  ", S, X, NULL, X);       /* s <- x' * x  */
  MatWrite (stdout , S, fmt,"S:");
  printf("\n");

  printf("  x1 = MatEigen (V, D, E, X);  /*  Eigen: X = V * diag(D,E) * V'   */\n");
  x1 = MatEigen (V, D, E, X);  /*  Eigen: X = V * diag(D,E) * V'   */
  printf("  /* eigenvalues = (D,E), D = the real part, E = the imaginary part */\n");
  /* eigenvalues = (D,E), D = the real part, E = the imaginary part */
  MatWrite (stdout , V, fmt,"V:");
  MatWrite (stdout , D, fmt,"D:");
  MatWrite (stdout , E, fmt,"E:");
  printf("\n");

  printf("  printf(\"Return code = %%Lg \\n\", x1 );\n");
  printf("Return code = %Lg \n", x1 );
  printf("  MatWMul (\" wt\", Y, V, D, V);   /*  y <- v * diag(d) * v'   */\n");
  MatWMul (" wt", Y, V, D, V);   /*  y <- v * diag(d) * v'   */
  MatWrite (stdout , Y, fmt,"Y:");
  printf("\n");

  printf("\n");

  printf("  MatSub(Y, X, Y);              /* y <- y - x   */\n");
  MatSub(Y, X, Y);              /* y <- y - x   */
  printf("  x1 = MatSumAbs(Y);\n");
  x1 = MatSumAbs(Y);
  printf("  printf(\"MatEigen Error: %%Lg \\n\", x1 );\n");
  printf("MatEigen Error: %Lg \n", x1 );
  printf("\n");

  printf("/*******************************************************/\n");
/*******************************************************/
  printf("\n");

  printf("  /* Roots of Polynomials */\n");
  /* Roots of Polynomials */
  printf("\n");

  printf("  /* p(x) = -12 + 22 x - 12 x^2 + 2 x^3 */\n");
  /* p(x) = -12 + 22 x - 12 x^2 + 2 x^3 */
  printf("  /* roots: 1, 2, 3 */\n");
  /* roots: 1, 2, 3 */
  printf("  MatReDim(B, 4, 1);\n");
  MatReDim(B, 4, 1);
  printf("  MatSet0(B, 0, 0,-12);\n");
  MatSet0(B, 0, 0,-12);
  printf("  MatSet0(B, 1, 0, 22);\n");
  MatSet0(B, 1, 0, 22);
  printf("  MatSet0(B, 2, 0,-12);\n");
  MatSet0(B, 2, 0,-12);
  printf("  MatSet0(B, 3, 0,  2);\n");
  MatSet0(B, 3, 0,  2);
  MatWrite (stdout , B, fmt,"B:");
  printf("\n");

  printf("  x1 = MatCompanion(C, B);\n");
  x1 = MatCompanion(C, B);
  MatWrite (stdout , C, fmt,"C:");
  printf("  printf(\"Return code = %%Lg \\n\", x1 );\n");
  printf("Return code = %Lg \n", x1 );
  printf("\n");

  printf("  x1 = MatEigen (V, D, E, C);  /*  Eigen: C = V * diag(D,E) * V'   */\n");
  x1 = MatEigen (V, D, E, C);  /*  Eigen: C = V * diag(D,E) * V'   */
  MatWrite (stdout , D, fmt,"D:");
  MatWrite (stdout , E, fmt,"E:");
  printf("  printf(\"Return code = %%Lg \\n\", x1 );\n");
  printf("Return code = %Lg \n", x1 );
  printf("/*******************************************************/\n");
/*******************************************************/
  printf("\n");

  printf("/* Roots of Polynomials */\n");
/* Roots of Polynomials */
  printf("\n");

  printf("  /* p(x) = -26 + 34 x - 10 x^2 + 2 x^3 */\n");
  /* p(x) = -26 + 34 x - 10 x^2 + 2 x^3 */
  printf("  /* roots: 1, 2 + 3 i, 2 - 3 i */\n");
  /* roots: 1, 2 + 3 i, 2 - 3 i */
  printf("  MatReDim(B, 4, 1);\n");
  MatReDim(B, 4, 1);
  printf("  MatSet0(B, 0, 0,-26);\n");
  MatSet0(B, 0, 0,-26);
  printf("  MatSet0(B, 1, 0, 34);\n");
  MatSet0(B, 1, 0, 34);
  printf("  MatSet0(B, 2, 0,-10);\n");
  MatSet0(B, 2, 0,-10);
  printf("  MatSet0(B, 3, 0,  2);\n");
  MatSet0(B, 3, 0,  2);
  MatWrite (stdout , B, fmt,"B:");
  printf("\n");

  printf("  x1 = MatCompanion(C, B);\n");
  x1 = MatCompanion(C, B);
  MatWrite (stdout , C, fmt,"C:");
  printf("  printf(\"Return code = %%Lg \\n\", x1 );\n");
  printf("Return code = %Lg \n", x1 );
  printf("\n");

  printf("  x1 = MatEigen (V, D, E, C);  /*  Eigen: C = V * diag(D,E) * V'   */\n");
  x1 = MatEigen (V, D, E, C);  /*  Eigen: C = V * diag(D,E) * V'   */
  MatWrite (stdout , D, fmt,"D:");
  MatWrite (stdout , E, fmt,"E:");
  printf("  printf(\"Return code = %%Lg \\n\", x1 );\n");
  printf("Return code = %Lg \n", x1 );
  printf("\n");

  printf("/*******************************************************/\n");
/*******************************************************/
  printf("\n");

  printf("  /* free  memory */\n");
  /* free  memory */
  printf("\n");

  printf("  MatUnDim(Z);\n");
  MatUnDim(Z);
  printf("  MatUnDim(Y);\n");
  MatUnDim(Y);
  printf("  MatUnDim(X);\n");
  MatUnDim(X);
  printf("  MatUnDim(W);\n");
  MatUnDim(W);
  printf("  MatUnDim(V);\n");
  MatUnDim(V);
  printf("  MatUnDim(U);\n");
  MatUnDim(U);
  printf("  MatUnDim(T);\n");
  MatUnDim(T);
  printf("  MatUnDim(S);\n");
  MatUnDim(S);
  printf("  MatUnDim(R);\n");
  MatUnDim(R);
  printf("  MatUnDim(Q);\n");
  MatUnDim(Q);
  printf("  MatUnDim(L);\n");
  MatUnDim(L);
  printf("  MatUnDim(J);\n");
  MatUnDim(J);
  printf("  MatUnDim(I);\n");
  MatUnDim(I);
  printf("  MatUnDim(H);\n");
  MatUnDim(H);
  printf("  MatUnDim(E);\n");
  MatUnDim(E);
  printf("  MatUnDim(D);\n");
  MatUnDim(D);
  printf("  MatUnDim(C);\n");
  MatUnDim(C);
  printf("  MatUnDim(B);\n");
  MatUnDim(B);
  printf("  MatUnDim(A);\n");
  MatUnDim(A);
  return 0;
}