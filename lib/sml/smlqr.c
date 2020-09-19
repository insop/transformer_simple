/* smlqr: QR using Householder reflectors */
/*****************************************************************/

#include <math.h>
#include "smlcomp.h"

/*****************************************************************/

/**
<p>
    Classical QR Decompisition:
   for an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
   orthogonal matrix Q and an n-by-n upper triangular matrix R so that
   A = Q*R.
<P>
   The QR decompostion always exists, even if the matrix does not have
   full rank, so the constructor will never fail.  The primary use of the
   QR decomposition is in the least squares solution of nonsquare systems
   of simultaneous linear equations.  This will fail if isFullRank()
   returns 0 (false).

<p>
    The Q and R factors can be retrived via the getQ() and getR()
    methods. Furthermore, a solve() method is provided to find the
    least squares solution of Ax=b using the QR factors.

   <p>
    (Adapted from JAMA, a Java Matrix Library, developed by jointly
    by the Mathworks and NIST; see  http://math.nist.gov/javanumerics/jama).
*/

   /** Array for internal storage of decomposition.
   @serial internal array storage.
   Array2D<REAL> QR;
   */


   /** Row and column dimensions.
   @serial column dimension.
   @serial row dimension.
   size_t m, n;
   */

   /** Array for internal storage of diagonal of R.
   @serial diagonal of R.
   Array1D<REAL> Rdiag;
   */

/*****************************************************************/

/**
 Generate and return the (economy-sized) orthogonal factor
@param     Q the (ecnomy-sized) orthogonal factor (Q*R=A).
 and the upper triangular factor, R, of the QR factorization
@param     R
*/

int MatQR (MATRIX Q,  MATRIX R, MATRIX A)
{
  MATRIX QR=MatDim(0,0);
  MATRIX Rdiag=MatDim(0,0);
  int rc = MatQRgetQR(QR, Rdiag, A);

  if (!rc) rc = MatQRgetQ (Q, QR);
  if (!rc) rc = MatQRgetR (R, QR, Rdiag);
  MatUnDim(QR);
  MatUnDim(Rdiag);
  return rc;
}

/**
    Create a QR factorization object for A.

    @param A rectangular (m>=n) matrix.
*/
int MatQRgetQR (MATRIX QR, MATRIX Rdiag_, MATRIX A)
{
   size_t m = MatRows(A);
   size_t n = MatCols(A);
   REAL*  Rdiag;
   size_t i, j, k;

   if (m < n) return 1;
   if (A  == Rdiag_ || QR == Rdiag_)  return 2;
   MatCopy(QR, A);        if (MatErr()) return 2;
   MatReDim(Rdiag_, n, 1); if (MatErr()) return 2;
   Rdiag  = MatData(Rdiag_);

   // Main loop.
   for (k = 0; k < n; k++) {
      // Compute 2-norm of k-th column without under/overflow.
      REAL nrm = 0;
      for (i = k; i < m; i++) {
         nrm = sml_hypot(nrm,MatGet0(QR,i,k));
      }

      if (nrm != 0.0) {
         // Form k-th Householder vector.
         if (MatGet0(QR,k,k) < 0) {
            nrm = -nrm;
         }
         for (i = k; i < m; i++) {
            MatSetDE0(QR,i,k, nrm);
         }
         MatSetPE0(QR,k,k,1.0);

         // Apply transformation to remaining columns.
         for (j = k+1; j < n; j++) {
            REAL s = 0.0;
            for (i = k; i < m; i++) {
               s += MatGet0(QR,i,k)*MatGet0(QR,i,j);
            }
            s = -s/MatGet0(QR,k,k);
            for (i = k; i < m; i++) {
               MatSetPE0(QR,i,j,  s*MatGet0(QR,i,k));
            }
         }
      }
      Rdiag[k] = -nrm;
   }
   return 0;
}


/**
    Flag to denote the matrix is of full rank.

    @return 1 if matrix is full rank, 0 otherwise.
*/
int isFullRank(MATRIX Rdiag_)
{
  REAL*  Rdiag = MatData(Rdiag_);
  size_t n = MatCols(Rdiag_);
  size_t j;
  for (j = 0; j != n; j++)  if (Rdiag[j] == 0)  return 0;
  return 1;
}




/**
 Generate and return the (economy-sized) orthogonal factor
@param     Q the (ecnomy-sized) orthogonal factor (Q*R=A).
*/

int MatQRgetQ (MATRIX Q, MATRIX QR)
 {
   size_t m = MatRows(QR);
   size_t n = MatCols(QR);
   size_t k = n;

   MatReDim(Q, m, n); if (MatErr()) return 2;

   if (!m || !n) return 0;

   // original code:  for (int k = n-1; k >= 0; k--) {
   do {
      size_t i, j;
      k--;
      for (i = 0; i < m; i++)    MatSet0(Q,i,k, 0.0);
      MatSet0(Q,k,k,  1.0);
      for (j = k; j < n; j++) {
         if (MatGet0(QR,k,k) != 0) {
            REAL s = 0.0;
            for (i = k; i < m; i++)
               s += MatGet0(QR,i,k)*MatGet0(Q,i,j);
            s = -s/MatGet0(QR,k,k);
            for (i = k; i < m; i++)
               MatSetPE0(Q,i,j,  s*MatGet0(QR,i,k));
         }
      }
   } while (k);
   return 0;
}



/** Return the upper triangular factor, R, of the QR factorization
@return     R
*/

int MatQRgetR (MATRIX R, MATRIX QR, MATRIX Rdiag_)
{
   size_t n = MatCols(QR);
   REAL*  Rdiag = MatData(Rdiag_);
   size_t i, j;

   MatReDim(R, n, n); if (MatErr()) return 2;

   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         if (i < j)
            MatSet0(R,i,j, MatGet0(QR,i,j));
         else if (i == j)
            MatSet0(R,i,j,  Rdiag[i]);
         else
            MatSet0(R,i,j,   0.0);
      }
   }
   return 0;
}


/**

Retreive the Householder vectors from QR factorization
@returns lower trapezoidal matrix whose columns define the reflections
*/

int  MatQRgetH (MATRIX H, MATRIX QR)
{
   size_t m = MatRows(QR);
   size_t n = MatCols(QR);
   size_t i, j;

   MatReDim(H, m, n); if (MatErr()) return 2;

   /* note: H is completely filled in by algorithm, so
      initializaiton of H is not necessary.
   */
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++)  {
         if (i >= j)
            MatSet0(H,i,j,  MatGet0(QR,i,j));
         else
            MatSet0(H,i,j, 0.0);
      }
   }
   return 0;
}



/** Least squares solution of A*x = b
@param B     m-length array (vector).
@return x    n-length array (vector) that minimizes the two norm of Q*R*X-B.
     If B is non-conformant, or if QR.isFullRank() is false,
                     the routine returns a null (0-length) vector.
*/

static int MatQRsolve1 (MATRIX x_, MATRIX b, MATRIX QR, MATRIX Rdiag_)
{
   size_t m = MatRows(QR);
   size_t n = MatCols(QR);
   REAL*  x;
   REAL*  Rdiag = MatData(Rdiag_);
   size_t  k;

   if (!m || !n) return 0;

   if (x_ == QR || x_ == Rdiag_) return 1;

   if (MatRows(b) != m)   return 1;      /* arrays must be conformant */

   if ( !isFullRank(Rdiag_) )    return 2;    /* matrix is rank deficient */

   MatCopy(x_, b);   if (MatErr()) return 2;
   x = MatData(x_);

   // Compute y = transpose(Q)*b
   for (k = 0; k < n; k++)  {
         size_t i;
         REAL s = 0.0;
         for (i = k; i < m; i++)
         {
            s += MatGet0(QR,i,k)*x[i];
         }
         s = -s/MatGet0(QR,k,k);
         for (i = k; i < m; i++)
         {
            x[i] += s*MatGet0(QR,i,k);
         }
   }
   // Solve R*x = y;
   // original code: for (int k = n-1; k >= 0; k--)  {
   k = n;    // we have checked that n > 0
   do {
      size_t i;
      k--;
      x[k] /= Rdiag[k];
      for (i = 0; i < k; i++) {
            x[i] -= x[k]*MatGet0(QR,i,k);
      }
   } while (k);


   /* return n x 1 portion of x */
   MatReDim(x_, n, 1);

   return 0;
}


/** Least squares solution of A*X = B
@param B     m x k Array (must conform).
@return X     n x k Array that minimizes the two norm of Q*R*X-B. If
                     B is non-conformant, or if QR.isFullRank() is false,
                     the routine returns a null (0x0) array.
*/

static int MatQRsolve2 (MATRIX X, MATRIX B, MATRIX QR, MATRIX Rdiag_)
{
   MATRIX XX;
   size_t m = MatRows(QR);
   size_t n = MatCols(QR);
   size_t nx= MatCols(B);
   REAL*  Rdiag = MatData(Rdiag_);
   size_t  i, j, k;

   if (!m || !n) return 0;

   if (X == QR || X == Rdiag_) return 1;

   if (MatRows(B) != m)   return 1;      /* arrays must be conformant */

   if ( !isFullRank(Rdiag_) )    return 2;    /* matrix is rank deficient */

   XX = MatDim(0,0);        if (MatErr()) return 2;
   MatCopy(XX, B);          if (MatErr()) return 2;
   MatReDim(X, n, nx);      if (MatErr()) return 2;


   // Compute Y = transpose(Q)*B
   for (k = 0; k < n; k++) {
      size_t j;
      for (j = 0; j < nx; j++) {
         size_t i;
         REAL s = 0.0;
         for (i = k; i < m; i++) {
            s += MatGet0(QR,i,k)*MatGet0(XX,i,j);
         }
         s = -s/MatGet0(QR,k,k);
         for (i = k; i < m; i++) {
            MatSetPE0(XX,i,j, s*MatGet0(QR,i,k));
         }
      }
   }
   // Solve R*X = Y;
   // original code: for (int k = n-1; k >= 0; k--) {
   k = n;       // we checked n > 0
   do {
      size_t i, j;
      k--;
      for (j = 0; j < nx; j++) {
         MatSetDE0(XX,k,j,  Rdiag[k]);
      }
      for (i = 0; i < k; i++) {
         for (j = 0; j < nx; j++) {
            MatSetME0(XX,i,j,  MatGet0(XX,k,j)*MatGet0(QR,i,k));
         }
      }
   }  while (k);


   /* return n x nx portion of X */
   for (i=0; i<n; i++)
     for (j=0; j<nx; j++)
         MatSet0(X,i,j, MatGet0(XX,i,j));

   MatUnDim(XX);
   return 0;
}


int MatQRsolve (MATRIX X, MATRIX B, MATRIX QR, MATRIX Rdiag_)
{
   if (MatCols(B) == 1) return MatQRsolve1(X, B, QR, Rdiag_);
   else                 return MatQRsolve2(X, B, QR, Rdiag_);
}
