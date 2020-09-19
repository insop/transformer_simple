/*
smlcompv.c: 1-dim versions that bypass MatGet/MatSet
*/
/******************************************************************/

#define SML_NAN       ((REAL)0.0/(REAL)0.0)
#define VEC(A)        (MatData(A), MatRows(A) *  MatCols(A))
#define VECARG(p, n)  (REAL *p, size_t n)

/*----------------------------------------------*/
/*
Functions that compute a scalar function of a matrix ->
*/
/*----------------------------------------------*/

static REAL  VecMax VECARG(p, n)            /*  max  element  */
{
  size_t i;
  REAL   s;


  if (n==0) return SML_NAN;

  s = p[0];
  for (i=1; i != n; ++i) {
      if ( p[i] > s)  s = p[i];
  }

  return s;
}

/*----------------------------------------------*/

REAL  MatMax (MATRIX A)             /*  max element  */
{
  return VecMax VEC(A);
}

/*----------------------------------------------*/

static REAL  VecMin VECARG(p, n)            /*  min  element  */
{
  size_t i;
  REAL   s;


  if (n==0) return SML_NAN;

  s = p[0];
  for (i=1; i != n; ++i) {
      if ( p[i] < s)  s = p[i];
  }

  return s;
}

/*----------------------------------------------*/

REAL  MatMin (MATRIX A)             /*  min element  */
{
  return VecMin VEC(A);
}

/*----------------------------------------------*/

static REAL  VecMaxAbs VECARG(p, n)            /*  max |element| */
{
  size_t i;
  REAL   s;


  if (n==0) return SML_NAN;

  s = p[0];
  for (i=1; i != n; ++i) {
      REAL   x;
      if ( (x=sml_fabs(p[i]))  > s)  s = x;
    }

  return s;
}

/*----------------------------------------------*/

REAL  MatMaxAbs (MATRIX A)             /*  max |element|  */
{
  return VecMaxAbs VEC(A);
}

/*----------------------------------------------*/

static REAL  VecMinAbs VECARG(p, n)            /*  min |element| */
{
  size_t i;
  REAL   s;


  if (n==0) return SML_NAN;

  s = p[0];
  for (i=1; i != n; ++i) {
      REAL   x;
      if ( (x=sml_fabs(p[i])) < s)  s = x;
    }

  return s;
}

/*----------------------------------------------*/

REAL  MatMinAbs (MATRIX A)             /*  min |element|  */
{
  return VecMinAbs VEC(A);
}

/*----------------------------------------------*/

static REAL2  VecSum  VECARG(p, n)   /*  sum of p[0:n-1] */
{
  REAL2 s=0;
  size_t i;

  for (i=0; i != n; s += p[i++]);

  return s;
}

/*----------------------------------------------*/

REAL2 MatSum (MATRIX A)             /*  sum of elements */
{
  return VecSum VEC(A);
}
/*----------------------------------------------*/
static REAL2  VecSumAbs VECARG(p, n)   /*  sum of |p[0:n-1]| */
{
  REAL2 s=0;
  size_t i;

  for (i=0; i != n; s += sml_fabs(p[i++]));

  return s;
}

/*----------------------------------------------*/

REAL2  MatSumAbs (MATRIX A)             /*  sum of |elements| */
{
  return VecSumAbs VEC(A);
}

/*----------------------------------------------*/

static REAL2  VecSS2  (REAL *r, REAL *p, size_t n)
/*
Returns sum of squares of p[0:n-1]
r <- a 2x1 matrix that stores a scaled result,
   sum of squares of a = r[0] * r[0] * r[1]

*/

{
  size_t i;
  REAL2  s=0;
  REAL   scale;

  scale = (n == 0) ? 0 : VecMaxAbs(p, n);

  if (scale > 0)  /* we add to s */
    for (i=0; i != n; ++i)  {
        REAL  x = p[i] / scale;
        s += x*x;
      }

  r[0] = scale;
  r[1] = s;

  return scale * scale * s;
}
/*----------------------------------------------*/

REAL2  MatSS  (MATRIX A)             /*  sum of squares */
{
  REAL    t[2];
  return  VecSS2 (t, MatData(A), MatRows(A) *  MatCols(A));
}
/*----------------------------------------------*/

REAL2  MatSS2  (MATRIX R, MATRIX A)
/*
Returns sum of squares of A
R <- a 2x1 matrix that stores a scaled result,
   sum of squares of A = R[1,1] * R[1,1] * R[2,1]
R == A is OK

*/

{
  REAL  t[2];
  REAL2 s = VecSS2 (t, MatData(A), MatRows(A) *  MatCols(A));

  MatReDim(R, 2, 1);
  if (!MatErr()) {
    MatSet0(R, 0, 0,  t[0]);
    MatSet0(R, 1, 0,  t[1]);
  }

  return s;
}

/*----------------------------------------------*/

static void VecFill(REAL *p, size_t n, REAL x)
{
  size_t i;

  for (i=0; i != n; p[i++] = x);

}

/*----------------------------------------------*/

void MatFill  (MATRIX A, REAL x)
   /* assign the value x to every element of A */
{
  VecFill (MatData(A), MatRows(A) *  MatCols(A), x);
}

/*----------------------------------------------*/
/*
Functions that compute one matrix as a function of another ->
*/
/*----------------------------------------------*/

void MatCopy (MATRIX A, MATRIX B)      /* A = B  */
  /* deep copy with redimensioning, nothing done if A==B  */
{

  size_t m = MatRows(B), n = MatCols(B);

  if (A==B) return;
  MatReDim(A, m, n);
  if (MatErr()) return;

  memmove (MatData(A), MatData(B), m * n * sizeof(REAL));
}

/*----------------------------------------------*/

static REAL2  VecCSum (REAL* a, size_t n, REAL* b) /* a = cusum of b[0:n-1] */
{
  REAL2 s=0;
  size_t i;

  if (a == b)
    for (i=0; i != n; ++i)  a[i] = (s += a[i]);
  else
    for (i=0; i != n; ++i)  a[i] = (s += b[i]);

  return s;
}

/*----------------------------------------------*/

REAL2 MatCSum (MATRIX A, MATRIX B)   /*  cum sum of elements */
/*
This function makes most sense for matrices of either 1 column or 1 row.
For other types of matrices, the result (A) depends on
the storage, by-row or by-column.
*/
{
  MatReDim(A, MatRows(B), MatCols(B));
  if (MatErr()) return SML_NAN;
  return VecCSum (MatData(A), MatRows(A) * MatCols(A), MatData(B));
}
/*----------------------------------------------*/

static REAL2 VecCMean (REAL* a, size_t n, REAL* b) /* a = cumean of b[0:n-1] */
{
  REAL2 s=0;
  size_t i;

  if (a == b)
    for (i=0; i != n; ++i)  a[i] = (s += a[i]) / (i+1);
  else
    for (i=0; i != n; ++i)  a[i] = (s += b[i]) / (i+1);

  return s / n;
}

/*----------------------------------------------*/

REAL2 MatCMean (MATRIX A, MATRIX B)   /*  cum mean of elements */
{
  MatReDim(A, MatRows(B), MatCols(B));
  if (MatErr()) return SML_NAN;
  return VecCMean (MatData(A), MatRows(A) * MatCols(A), MatData(B));
}
/*----------------------------------------------*/

static void VecApplyE (REAL* a, size_t n, REAL (*f)(REAL), REAL* b)

/* a[0:n-1] = f(b[0:n-1])  */
{
  size_t i;

  if (a == b)
    for (i=0; i != n; ++i)  a[i] = f(a[i]);
  else
    for (i=0; i != n; ++i)  a[i] = f(b[i]);
}
/*----------------------------------------------*/

void MatApplyE (MATRIX A, REAL (*f)(REAL), MATRIX B)

/* A=f(B) elementwise */
{
  MatReDim(A, MatRows(B), MatCols(B));
  if (MatErr()) return;

  VecApplyE(MatData(A), MatRows(A) * MatCols(A), f, MatData(B));
}
/*----------------------------------------------*/

static void VecExpE  (REAL* a, size_t n, REAL* b)

/* a[0:n-1] = exp(b[0:n-1])  */
{
  size_t i;

  if (a == b)
    for (i=0; i != n; ++i)  a[i] = exp(a[i]);
  else
    for (i=0; i != n; ++i)  a[i] = exp(b[i]);
}
/*----------------------------------------------*/

void MatExpE (MATRIX A, MATRIX B)

/* A=exp(B) elementwise */
{

  MatReDim(A, MatRows(B), MatCols(B));
  if (MatErr()) return;

  VecExpE(MatData(A), MatRows(A) * MatCols(A), MatData(B));
}
/*----------------------------------------------*/

static void VecLogE  (REAL* a, size_t n, REAL* b)

/* a[0:n-1] = log(b[0:n-1])  */
{
  size_t i;

  if (a == b)
    for (i=0; i != n; ++i)  a[i] = log(a[i]);
  else
    for (i=0; i != n; ++i)  a[i] = log(b[i]);
}
/*----------------------------------------------*/

void MatLogE (MATRIX A, MATRIX B)

/* A=log(B) elementwise */
{

  MatReDim(A, MatRows(B), MatCols(B));
  if (MatErr()) return;

  VecLogE(MatData(A), MatRows(A) * MatCols(A), MatData(B));
}
/*----------------------------------------------*/

static void VecRecipE  (REAL* a, size_t n, REAL* b)

/* a[0:n-1] = 1/(b[0:n-1])  */
{
  size_t i;


  if (a == b)
    for (i=0; i != n; ++i)  a[i] = 1.0/(a[i]);
  else
    for (i=0; i != n; ++i)  a[i] = 1.0/(b[i]);
}
/*----------------------------------------------*/

void MatRecipE (MATRIX A, MATRIX B)

 /* A = 1 / B el-wise */
{
  MatReDim(A, MatRows(B), MatCols(B));
  if (MatErr()) return;

  VecRecipE(MatData(A), MatRows(A) * MatCols(A), MatData(B));
}
/*----------------------------------------------*/

static void VecRecipE0  (REAL* a, size_t n, REAL* b)

/* a[0:n-1] = 1/(b[0:n-1]), 1/0 = 0  */
{
  size_t i;

  if (a == b) for (i=0; i != n; ++i) {
    if (a[i] != 0.0) a[i] = 1.0 / a[i];
    else a[i] = 0.0;
  }
  else  for (i=0; i != n; ++i) {
    if (b[i] != 0.0) a[i] = 1.0 / b[i];
    else a[i] = 0.0;
  }
}
/*----------------------------------------------*/

void MatRecipE0 (MATRIX A, MATRIX B)

 /* A = 1 / B el-wise, defining  1/0 = 0 */
{

  MatReDim(A, MatRows(B), MatCols(B));
  if (MatErr()) return;

  VecRecipE0(MatData(A), MatRows(A) * MatCols(A), MatData(B));
}
/*----------------------------------------------*/

static void VecSqrtE  (REAL* a, size_t n, REAL* b)

/* a[0:n-1] = sqrt(b[0:n-1])  */
{
  size_t i;

  if (a == b)
    for (i=0; i != n; ++i)  a[i] = sqrt(a[i]);
  else
    for (i=0; i != n; ++i)  a[i] = sqrt(b[i]);
}
/*----------------------------------------------*/

void MatSqrtE (MATRIX A, MATRIX B)

/* A=sqrt(B) elementwise */
{

  MatReDim(A, MatRows(B), MatCols(B));
  if (MatErr()) return;

  VecSqrtE(MatData(A), MatRows(A) * MatCols(A), MatData(B));
}
/*----------------------------------------------*/

static void VecMulE  (REAL* a, size_t n, REAL* b, REAL* c)

/* a[0:n-1] = b[0:n-1] * c[0:n-1] */
{
  size_t i;

  if (a == b && a == c)
    for (i=0; i != n; ++i)  a[i] *= a[i];
  else if (a == b)
    for (i=0; i != n; ++i)  a[i] *= c[i];
  else if (a == c)
    for (i=0; i != n; ++i)  a[i] *= b[i];
  else if (b == c)
    for (i=0; i != n; ++i)  a[i] = b[i] * b[i];
  else
    for (i=0; i != n; ++i)  a[i] = b[i] * c[i];
}
/*----------------------------------------------*/

void MatMulE (MATRIX A, MATRIX B, MATRIX C)

/* A= B*C elementwise */
{

  if (MatCols(B) !=  MatCols(C) || MatRows(B) !=  MatRows(C)) {
    MatErrThrow("MatMulE: unconformable matrices.");
    return;
  }

  MatReDim(A, MatRows(B), MatCols(B));
  if (MatErr()) return;

  VecMulE(MatData(A), MatRows(A) * MatCols(A), MatData(B), MatData(C));
}
/*----------------------------------------------*/

static void VecDivE (REAL* a, size_t n, REAL* b, REAL* c)

/* a[0:n-1] = b[0:n-1] / c[0:n-1] */
{
  size_t i;

  if (a == b && a == c)
    for (i=0; i != n; ++i)  a[i] /= a[i]; /* note: 0 / 0 != 1 */
  else if (a == b)
    for (i=0; i != n; ++i)  a[i] /= c[i];
  else if (a == c)
    for (i=0; i != n; ++i)  a[i] = b[i] / a[i];
  else if (b == c)
    for (i=0; i != n; ++i)  a[i] = b[i] / b[i]; /* note: 0 / 0 != 1 */
  else
    for (i=0; i != n; ++i)  a[i] = b[i] / c[i];
}
/*----------------------------------------------*/

void MatDivE (MATRIX A, MATRIX B, MATRIX C)      /* A = B / C el-wise*/

{
  if (MatCols(B) !=  MatCols(C) || MatRows(B) !=  MatRows(C) )   {
    MatErrThrow("MatDivE: unconformable matrices.");
    return;
  }

  MatReDim(A, MatRows(B), MatCols(B));
  if (MatErr()) return;

  VecDivE(MatData(A), MatRows(A) * MatCols(A), MatData(B), MatData(C));
}
/*----------------------------------------------*/

void MatDivE0 (MATRIX A, MATRIX B, MATRIX C)
 /* A = B / C el-wise, defining  1/0 = 0 */

{
  MatRecipE0 (A, C);
  MatMulE    (A, A, B);
}
/*----------------------------------------------*/

static void VecAdd  (REAL* a, size_t n, REAL* b, REAL* c)

/* a[0:n-1] = b[0:n-1] + c[0:n-1] */
{
  size_t i;

  if (a == b && a == c)
    for (i=0; i != n; ++i)  a[i] += a[i];
  else if (a == b)
    for (i=0; i != n; ++i)  a[i] += c[i];
  else if (a == c)
    for (i=0; i != n; ++i)  a[i] += b[i];
  else if (b == c)
    for (i=0; i != n; ++i)  a[i] = b[i] + b[i];
  else
    for (i=0; i != n; ++i)  a[i] = b[i] + c[i];
}
/*----------------------------------------------*/

void MatAdd (MATRIX A, MATRIX B, MATRIX C)      /* A = B + C */

{
  if (MatCols(B) !=  MatCols(C) || MatRows(B) !=  MatRows(C) )   {
    MatErrThrow("MatAdd: unconformable matrices.");
    return;
  }

  MatReDim(A, MatRows(B), MatCols(C));
  if (MatErr()) return;

  VecAdd(MatData(A), MatRows(A) * MatCols(A), MatData(B), MatData(C));
}
/*----------------------------------------------*/

static void VecSub  (REAL* a, size_t n, REAL* b, REAL* c)

/* a[0:n-1] = b[0:n-1] - c[0:n-1] */
{
  size_t i;

  if (a == b && a == c)
    for (i=0; i != n; ++i)  a[i] -= a[i]; /* note: nan - nan != 0 */
  else if (a == b)
    for (i=0; i != n; ++i)  a[i] -= c[i];
  else if (a == c)
    for (i=0; i != n; ++i)  a[i] = b[i] - a[i];
  else if (b == c)
    for (i=0; i != n; ++i)  a[i] = b[i] - b[i]; /* note: nan - nan != 0 */
  else
    for (i=0; i != n; ++i)  a[i] = b[i] - c[i];
}
/*----------------------------------------------*/

void MatSub (MATRIX A, MATRIX B, MATRIX C)      /* A = B - C */

{
  if (MatCols(B) !=  MatCols(C) || MatRows(B) !=  MatRows(C) )   {
    MatErrThrow("MatSub: unconformable matrices.");
    return;
  }

  MatReDim(A, MatRows(B), MatCols(B));
  if (MatErr()) return;

  VecSub(MatData(A), MatRows(A) * MatCols(A), MatData(B), MatData(C));
}
/*----------------------------------------------*/

static void VecAddScalar (REAL* a, size_t n, REAL* b, REAL x)

/* a[0:n-1] = b[0:n-1] + x */
{
  size_t i;

  if (a == b)
    for (i=0; i != n; ++i)  a[i] += x;
  else
    for (i=0; i != n; ++i)  a[i] = b[i] + x;
}
/*----------------------------------------------*/
void MatAddScalar (MATRIX A, MATRIX B, REAL x)      /* A = B + x */
{
  MatReDim(A, MatRows(B), MatCols(B));
  if (MatErr()) return;
  VecAddScalar(MatData(A), MatRows(A) * MatCols(A), MatData(B), x);
}
/*----------------------------------------------*/

static void VecMulScalar (REAL* a, size_t n, REAL* b, REAL x)

/* a[0:n-1] = b[0:n-1] * x */
{
  size_t i;

  if (a == b)
    for (i=0; i != n; ++i)  a[i] *= x;
  else
    for (i=0; i != n; ++i)  a[i] = b[i] * x;
}
/*----------------------------------------------*/
void MatMulScalar (MATRIX A, MATRIX B, REAL x)      /* A = B + x */
{
  MatReDim(A, MatRows(B), MatCols(B));
  if (MatErr()) return;
  VecMulScalar(MatData(A), MatRows(A) * MatCols(A), MatData(B), x);
}
/*----------------------------------------------*/
static void sort(REAL a[], size_t n)
/*   shellsort a[0:n-1] */
{
  size_t h;

  /* up to? 1, 4, 13, 40, 121, ...*/
  for (h = 1; h <= (n-1)/9; h = 3*h+1);

  for ( ; h > 0; h /= 3)  {      /* ..., 1039, 364, 121, 40, 13, 4, 1 */
    size_t i;
    for (i = h; i < n; ++i)   {  /* insertion sort */
        size_t j = i; REAL v = a[i];
        while (j >= h &&   v < a[j-h])
            { a[j] = a[j-h]; j -= h; }
        a[j] = v;
    }
  }
}
/*----------------------------------------------*/
void MatSort (MATRIX A)
/*
Sort A, result depends on layout,
Most useful if A is a (row or column) vector
*/
{
  if (MatRows(A)*MatCols(A) <= 1) return;
  sort(MatData(A), MatRows(A)*MatCols(A));
}
/*----------------------------------------------*/
