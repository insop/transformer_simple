/* smlsvd: Singular Value Decomposition. */
/*****************************************************************/

#include <math.h>
#include "smlcomp.h"

/*****************************************************************/
/*
For an m-by-n matrix A with m >= n, the singular value decomposition is
an m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and
an n-by-n orthogonal matrix V so that A = U*S*V'.

The singular values, sigma[k] = S[k][k], are ordered so that
sigma[0] >= sigma[1] >= ... >= sigma[n-1].

The singular value decompostion always exists.
The matrix condition number and the effective numerical
rank can be computed from this decomposition.

(Adapted from JAMA, a Java Matrix Library, developed by jointly by the
Mathworks and NIST; see http://math.nist.gov/javanumerics/jama).
*/
/************************************************************/

int  MatSVD (MATRIX U, MATRIX S, MATRIX V, MATRIX AA)

{
  MATRIX  A, e_, work_ ;
  REAL  *s , *e , *work  ;

  int m  = MatRows(AA);
  int n  = MatCols(AA);
  int nu = Matmini(m,n);
  int nct = Matmini(m-1,n);
  int nrt = Matmaxi(0,Matmini(n-2,m));
  int p   = Matmini(n,m+1);
  int pp = p-1;
  int iter = 0;
  int wantu = 1;  /* boolean */
  int wantv = 1;  /* boolean */
  int i, j, k;
  const REAL   eps = REAL_EPSILON;


  if (U == S || U == V || U == AA || S == V || S == AA || V == AA) {
    MatErrThrow("MatSVD: arguments must be distinct.");
    return 1;
  }

  MatReDim(S, Matmini(m+1,n), 1);
  if (MatErr()) return 1;
  s = MatData(S);

  MatReDim(U, m, nu);
  if (MatErr()) return 1;
  MatFill(U, 0);

  MatReDim(V, n, n);
  if (MatErr()) return 1;

  e_ = MatDim(n, 1);
  if (MatErr()) return 1;
  e  = MatData(e_);

  work_ = MatDim(m, 1);
  if (MatErr()) return 1;
  work  = MatData(work_);

  A = MatDim(m, n);
  if (MatErr()) return 1;
  MatCopy(A, AA);


  // Reduce A to bidiagonal form, storing the diagonal elements
  // in s and the super-diagonal elements in e.

  for (k = 0; k < Matmaxi(nct,nrt); k++) {
     if (k < nct) {

        // Compute the transformation for the k-th column and
        // place the k-th diagonal in s[k].
        // Compute 2-norm of k-th column without under/overflow.
        s[k] = 0;
        for (i = k; i < m; i++) {
           s[k] = sml_hypot(s[k],MatGet0(A,i,k));
        }
        if (s[k] != 0.0) {
           if (MatGet0(A,k,k) < 0.0) {
              s[k] = -s[k];
           }
           for (i = k; i < m; i++) {
              MatSetDE0(A,i,k, s[k]);
           }
           MatSetPE0(A,k,k, 1.0);
        }
        s[k] = -s[k];
     }
     for (j = k+1; j < n; j++) {
        if ((k < nct) && (s[k] != 0.0))  {

        // Apply the transformation.

           REAL2  t = 0;
           for (i = k; i < m; i++) {
              t += MatGet0(A,i,k)*MatGet0(A,i,j);
           }
           t = -t/MatGet0(A,k,k);
           for (i = k; i < m; i++) {
              MatSetPE0(A,i,j, t*MatGet0(A,i,k));
           }
        }

        // Place the k-th row of A into e for the
        // subsequent calculation of the row transformation.

        e[j] = MatGet0(A,k,j);
     }
     if (wantu && (k < nct)) {

        // Place the transformation in U for subsequent back
        // multiplication.

        for (i = k; i < m; i++) {
           MatSet0(U,i,k, MatGet0(A,i,k));
        }
     }
     if (k < nrt) {

        // Compute the k-th row transformation and place the
        // k-th super-diagonal in e[k].
        // Compute 2-norm without under/overflow.
        e[k] = 0;
        for (i = k+1; i < n; i++) {
           e[k] = sml_hypot(e[k],e[i]);
        }
        if (e[k] != 0.0) {
           if (e[k+1] < 0.0) {
              e[k] = -e[k];
           }
           for (i = k+1; i < n; i++) {
              e[i] /= e[k];
           }
           e[k+1] += 1.0;
        }
        e[k] = -e[k];
        if ((k+1 < m) && (e[k] != 0.0)) {

        // Apply the transformation.

           for (i = k+1; i < m; i++) {
              work[i] = 0.0;
           }
           for (j = k+1; j < n; j++) {
              for (i = k+1; i < m; i++) {
                 work[i] += e[j]*MatGet0(A,i,j);
              }
           }
           for (j = k+1; j < n; j++) {
              REAL t = -e[j]/e[k+1];
              for (i = k+1; i < m; i++) {
                 MatSetPE0(A,i,j, t*work[i]);
              }
           }
        }
        if (wantv) {

        // Place the transformation in V for subsequent
        // back multiplication.

           for (i = k+1; i < n; i++) {
              MatSet0(V,i,k, e[i]);
           }
        }
     }
  }

  // Set up the final bidiagonal matrix or order p.

  if (nct < n) {
     s[nct] = MatGet0(A,nct,nct);
  }
  if (m < p) {
     s[p-1] = 0.0;
  }
  if (nrt+1 < p) {
     e[nrt] = MatGet0(A,nrt,p-1);
  }
  e[p-1] = 0.0;

  // If required, generate U.

  if (wantu) {
     for (j = nct; j < nu; j++) {
        for (i = 0; i < m; i++) {
           MatSet0(U,i,j, 0.0);
        }
        MatSet0(U,j,j, 1.0);
     }
     for (k = nct-1; k >= 0; k--) {
        if (s[k] != 0.0) {
           for (j = k+1; j < nu; j++) {
              REAL2 t = 0;
              for (i = k; i < m; i++) {
                 t += MatGet0(U,i,k)*MatGet0(U,i,j);
              }
              t = -t/MatGet0(U,k,k);
              for (i = k; i < m; i++) {
                 MatSetPE0(U,i,j,  t*MatGet0(U,i,k));
              }
           }
           for (i = k; i < m; i++ ) {
              MatSet0(U,i,k,   -MatGet0(U,i,k));
           }
           MatSet0(U,k,k, 1.0 + MatGet0(U,k,k));
           for (i = 0; i < k-1; i++) {
              MatSet0(U,i,k, 0.0);
           }
        } else {
           for (i = 0; i < m; i++) {
              MatSet0(U,i,k, 0.0);
           }
           MatSet0(U,k,k,  1.0);
        }
     }
  }

  // If required, generate V.

  if (wantv) {
     for (k = n-1; k >= 0; k--) {
        if ((k < nrt) && (e[k] != 0.0)) {
           for (j = k+1; j < nu; j++) {
              REAL2 t = 0;
              for (i = k+1; i < n; i++) {
                 t += MatGet0(V,i,k)*MatGet0(V,i,j);
              }
              t = -t/MatGet0(V,k+1,k);
              for (i = k+1; i < n; i++) {
                 MatSetPE0(V,i,j, t*MatGet0(V,i,k));
              }
           }
        }
        for (i = 0; i < n; i++) {
           MatSet0(V,i,k,  0.0);
        }
        MatSet0(V,k,k, 1.0);
     }
  }

  // Main iteration loop for the singular values.

  while (p > 0) {
     int k;
     int kase=0;

     // Here is where a test for too many iterations would go.

     // This section of the program inspects for
     // negligible elements in the s and e arrays.  On
     // completion the variables kase and k are set as follows.

     // kase = 1     if s(p) and e[k-1] are negligible and k<p
     // kase = 2     if s(k) is negligible and k<p
     // kase = 3     if e[k-1] is negligible, k<p, and
     //              s(k), ..., s(p) are not negligible (qr step).
     // kase = 4     if e(p-1) is negligible (convergence).

     for (k = p-2; k >= -1; k--) {
        if (k == -1) {
           break;
        }
        if (sml_fabs(e[k]) <= eps*(sml_fabs(s[k]) + sml_fabs(s[k+1]))) {
           e[k] = 0.0;
           break;
        }
     }
     if (k == p-2) {
        kase = 4;
     } else {
        REAL t;
        int ks;
        for (ks = p-1; ks >= k; ks--) {
           if (ks == k) {
              break;
           }
           t = (ks != p ? sml_fabs(e[ks]) : 0.) +
                      (ks != k+1 ? sml_fabs(e[ks-1]) : 0.);
           if (sml_fabs(s[ks]) <= eps*t)  {
              s[ks] = 0.0;
              break;
           }
        }
        if (ks == k) {
           kase = 3;
        } else if (ks == p-1) {
           kase = 1;
        } else {
           kase = 2;
           k = ks;
        }
     }
     k++;

     // Perform the task indicated by kase.

     switch (kase) {

        // Deflate negligible s(p).

        case 1: {
           REAL f = e[p-2];
           e[p-2] = 0.0;
           for (j = p-2; j >= k; j--) {
              REAL t = sml_hypot(s[j],f);
              REAL cs = s[j]/t;
              REAL sn = f/t;
              s[j] = t;
              if (j != k) {
                 f = -sn*e[j-1];
                 e[j-1] = cs*e[j-1];
              }
              if (wantv) {
                 for (i = 0; i < n; i++) {
                    t = cs*MatGet0(V,i,j) + sn*MatGet0(V,i,p-1);
                    MatSet0(V,i,p-1,  -sn*MatGet0(V,i,j) + cs*MatGet0(V,i,p-1));
                    MatSet0(V,i,j, t);
                 }
              }
           }
        }
        break;

        // Split at negligible s(k).

        case 2: {
           REAL f = e[k-1];
           e[k-1] = 0.0;
           for (j = k; j < p; j++) {
              REAL t = sml_hypot(s[j],f);
              REAL cs = s[j]/t;
              REAL sn = f/t;
              s[j] = t;
              f = -sn*e[j];
              e[j] = cs*e[j];
              if (wantu) {
                 for (i = 0; i < m; i++) {
                    t = cs*MatGet0(U,i,j) + sn*MatGet0(U,i,k-1);
                    MatSet0(U,i,k-1,   -sn*MatGet0(U,i,j) + cs*MatGet0(U,i,k-1));
                    MatSet0(U,i,j,  t);
                 }
              }
           }
        }
        break;

        // Perform one qr step.

        case 3: {

           // Calculate the shift.

           REAL scale = Matmaxr(Matmaxr(Matmaxr(Matmaxr(
                   sml_fabs(s[p-1]),sml_fabs(s[p-2])),sml_fabs(e[p-2])),
                   sml_fabs(s[k])),sml_fabs(e[k]));
           REAL sp = s[p-1]/scale;
           REAL spm1 = s[p-2]/scale;
           REAL epm1 = e[p-2]/scale;
           REAL sk = s[k]/scale;
           REAL ek = e[k]/scale;
           REAL b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
           REAL c = (sp*epm1)*(sp*epm1);
           REAL shift = 0.0;
           REAL f, g;
           if ((b != 0.0) | (c != 0.0)) {
              shift = sml_sqrt(b*b + c);
              if (b < 0.0) {
                 shift = -shift;
              }
              shift = c/(b + shift);
           }
           f = (sk + sp)*(sk - sp) + shift;
           g = sk*ek;

           // Chase zeros.

           for (j = k; j < p-1; j++) {
              REAL t = sml_hypot(f,g);
              REAL cs = f/t;
              REAL sn = g/t;
              if (j != k) {
                 e[j-1] = t;
              }
              f = cs*s[j] + sn*e[j];
              e[j] = cs*e[j] - sn*s[j];
              g = sn*s[j+1];
              s[j+1] = cs*s[j+1];
              if (wantv) {
                 for (i = 0; i < n; i++) {
                    t = cs*MatGet0(V,i,j) + sn*MatGet0(V,i,j+1);
                    MatSet0(V,i,j+1,  -sn*MatGet0(V,i,j) + cs*MatGet0(V,i,j+1));
                    MatSet0(V,i,j, t);
                 }
              }
              t = sml_hypot(f,g);
              cs = f/t;
              sn = g/t;
              s[j] = t;
              f = cs*e[j] + sn*s[j+1];
              s[j+1] = -sn*e[j] + cs*s[j+1];
              g = sn*e[j+1];
              e[j+1] = cs*e[j+1];
              if (wantu && (j < m-1)) {
                 for (i = 0; i < m; i++) {
                    t = cs*MatGet0(U,i,j) + sn*MatGet0(U,i,j+1);
                    MatSet0(U,i,j+1,  -sn*MatGet0(U,i,j) + cs*MatGet0(U,i,j+1));
                    MatSet0(U,i,j, t);
                 }
              }
           }
           e[p-2] = f;
           iter = iter + 1;
        }
        break;

        // Convergence.

        case 4: {

           // Make the singular values positive.

           if (s[k] <= 0.0) {
              s[k] = (s[k] < 0.0 ? -s[k] : 0.0);
              if (wantv) {
                 for (i = 0; i <= pp; i++) {
                    MatSet0(V,i,k,  -MatGet0(V,i,k));
                 }
              }
           }

           // Order the singular values.

           while (k < pp) {
              REAL t;
              if (s[k] >= s[k+1]) {
                 break;
              }
              t = s[k];
              s[k] = s[k+1];
              s[k+1] = t;
              if (wantv && (k < n-1)) {
                 for (i = 0; i < n; i++) {
                    t = MatGet0(V,i,k+1); MatSet0(V,i,k+1, MatGet0(V,i,k)); MatSet0(V,i,k, t);
                 }
              }
              if (wantu && (k < m-1)) {
                 for (i = 0; i < m; i++) {
                    t = MatGet0(U,i,k+1); MatSet0(U,i,k+1, MatGet0(U,i,k)); MatSet0(U,i,k, t);
                 }
              }
              k++;
           }
           iter = 0;
           p--;
        }
        break;
     }
  }

MatUnDim(A);
MatUnDim(work_);
MatUnDim(e_);

return 0;

};

