/* smleig.c: Eigenvalues and eigenvectors of a real (non-complex) matrix. */
/**********************************************************************/

#include <math.h>
#include "smlcomp.h"

/**********************************************************************/
/**

    If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is
    diagonal and the eigenvector matrix V is orthogonal. That is,
    the diagonal values of D are the eigenvalues, and
    V*V' = I, where I is the identity matrix.  The columns of V
    represent the eigenvectors in the sense that A*V = V*D.

<P>
    If A is not symmetric, then the eigenvalue matrix D is block diagonal
    with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
    a + i*b, in 2-by-2 blocks, [a, b; -b, a].  That is, if the complex
    eigenvalues look like
<pre>

          u + iv     .        .          .      .    .
            .      u - iv     .          .      .    .
            .        .      a + ib       .      .    .
            .        .        .        a - ib   .    .
            .        .        .          .      x    .
            .        .        .          .      .    y
</pre>
        then D looks like
<pre>

            u        v        .          .      .    .
           -v        u        .          .      .    .
            .        .        a          b      .    .
            .        .       -b          a      .    .
            .        .        .          .      x    .
            .        .        .          .      .    y
</pre>
    This keeps V a real matrix in both symmetric and non-symmetric
    cases, and A*V = V*D.



    <p>
    The matrix V may be badly
    conditioned, or even singular, so the validity of the equation
    A = V*D*inverse(V) depends upon the condition number of V.

   <p>
 (Adapted from JAMA, a Java Matrix Library, developed by jointly
 by the Mathworks and NIST; see  http://math.nist.gov/javanumerics/jama).
**/
/**********************************************************************/


// Symmetric Householder reduction to tridiagonal form.
static void tred2 (MATRIX V, MATRIX d_, MATRIX e_)
{

//  This is derived from the Algol procedures tred2 by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

   int   n = MatRows(V);
   REAL* d = MatData(d_);
   REAL* e = MatData(e_);
   int i, j;


   for (j = 0; j < n; j++) {
      d[j] = MatGet0(V,n-1,j);
   }

   // Householder reduction to tridiagonal form.

   for (i = n-1; i > 0; i--) {

      // Scale to avoid under/overflow.

      REAL   scale = 0.0;
      REAL   h = 0.0;
      int k;
      for (k = 0; k < i; k++) {
         scale = scale + sml_fabs(d[k]);
      }
      if (scale == 0.0) {
         int j;
         e[i] = d[i-1];
         for (j = 0; j < i; j++) {
            d[j] = MatGet0(V,i-1,j);
            MatSet0(V,i,j, 0.0);
            MatSet0(V,j,i, 0.0);
         }
      } else {

         // Generate Householder vector.

         int k, j;
         REAL f, g, hh;
         for (k = 0; k < i; k++) {
            d[k] /= scale;
            h += d[k] * d[k];
         }
         f = d[i-1];
         g = sml_sqrt(h);
         if (f > 0) {
            g = -g;
         }
         e[i] = scale * g;
         h = h - f * g;
         d[i-1] = f - g;
         for (j = 0; j < i; j++) {
            e[j] = 0.0;
         }

         // Apply similarity transformation to remaining columns.

         for (j = 0; j < i; j++) {
            int k;
            f = d[j];
            MatSet0(V,j,i, f);
            g = e[j] + MatGet0(V,j,j) * f;
            for (k = j+1; k <= i-1; k++) {
               g += MatGet0(V,k,j) * d[k];
               e[k] += MatGet0(V,k,j) * f;
            }
            e[j] = g;
         }
         f = 0.0;
         for (j = 0; j < i; j++) {
            e[j] /= h;
            f += e[j] * d[j];
         }
         hh = f / (h + h);
         for (j = 0; j < i; j++) {
            e[j] -= hh * d[j];
         }
         for (j = 0; j < i; j++) {
            int k;
            f = d[j];
            g = e[j];
            for (k = j; k <= i-1; k++) {
               MatSetME0(V,k,j,  (f * e[k] + g * d[k]));
            }
            d[j] = MatGet0(V,i-1,j);
            MatSet0(V,i,j,  0.0);
         }
      }
      d[i] = h;
   }

   // Accumulate transformations.

   for (i = 0; i < n-1; i++) {
      REAL h;
      int k;
      MatSet0(V,n-1,i, MatGet0(V,i,i));
      MatSet0(V,i,i,  1.0);
      h = d[i+1];
      if (h != 0.0) {
         int k, j;
         for (k = 0; k <= i; k++) {
            d[k] = MatGet0(V,k,i+1) / h;
         }
         for (j = 0; j <= i; j++) {
            int k;
            REAL g = 0.0;
            for (k = 0; k <= i; k++) {
               g += MatGet0(V,k,i+1) * MatGet0(V,k,j);
            }
            for (k = 0; k <= i; k++) {
               MatSetME0(V,k,j,  g * d[k]);
            }
         }
      }
      for (k = 0; k <= i; k++) {
         MatSet0(V,k,i+1,  0.0);
      }
   }
   for (j = 0; j < n; j++) {
      d[j] = MatGet0(V,n-1,j);
      MatSet0(V,n-1,j, 0.0);
   }
   MatSet0(V,n-1,n-1,  1.0);
   e[0] = 0.0;
}

// Symmetric tridiagonal QL algorithm.

static void tql2 (MATRIX V, MATRIX d_, MATRIX e_)
{
//  This is derived from the Algol procedures tql2, by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

   int   n = MatRows(V);
   REAL* d = MatData(d_);
   REAL* e = MatData(e_);
   REAL f    = 0.0;
   REAL tst1 = 0.0;
   REAL eps  = REAL_EPSILON;
   int i, j, l;


   for (i = 1; i < n; i++) {
      e[i-1] = e[i];
   }
   e[n-1] = 0.0;

   for (l = 0; l < n; l++) {

      // Find small subdiagonal element

      int m = l;
      tst1 = Matmaxr(tst1,sml_fabs(d[l]) + sml_fabs(e[l]));

     // Original while-loop from Java code
      while (m < n) {
         if (sml_fabs(e[m]) <= eps*tst1) {
            break;
         }
         m++;
      }


      // If m == l, d[l] is an eigenvalue,
      // otherwise, iterate.

      if (m > l) {
         int iter = 0;
         do {
            REAL g, p, r,  dl1, h;
            REAL c=1.0, c2=1.0, c3=1.0, el1, s=0, s2=0;
            int i;
            iter++;           // (Could check iteration count here.)

            // Compute implicit shift

            g = d[l];
            p = (d[l+1] - g) / (2.0 * e[l]);
            r = sml_hypot(p,1.0);
            if (p < 0) {
               r = -r;
            }
            d[l] = e[l] / (p + r);
            d[l+1] = e[l] * (p + r);
            dl1 = d[l+1];
            h = g - d[l];
            for (i = l+2; i < n; i++) {
               d[i] -= h;
            }
            f = f + h;

            // Implicit QL transformation.

            p = d[m];
            el1 = e[l+1];
            for (i = m-1; i >= l; i--) {
               int k;
               c3 = c2;
               c2 = c;
               s2 = s;
               g = c * e[i];
               h = c * p;
               r = sml_hypot(p,e[i]);
               e[i+1] = s * r;
               s = e[i] / r;
               c = p / r;
               p = c * d[i] - s * g;
               d[i+1] = h + s * (c * g + s * d[i]);

               // Accumulate transformation.

               for (k = 0; k < n; k++) {
                  h = MatGet0(V,k,i+1);
                  MatSet0(V,k,i+1, s * MatGet0(V,k,i) + c * h);
                  MatSet0(V,k,i,  c * MatGet0(V,k,i) - s * h);
               }
            }
            p = -s * s2 * c3 * el1 * e[l] / dl1;
            e[l] = s * p;
            d[l] = c * p;

            // Check for convergence.

         } while (sml_fabs(e[l]) > eps*tst1);
      }
      d[l] = d[l] + f;
      e[l] = 0.0;
   }

   // Sort eigenvalues and corresponding vectors.

   for (i = 0; i < n-1; i++) {
      int k = i;
      int j;
      REAL p = d[i];
      for (j = i+1; j < n; j++) {
         if (d[j] < p) {
            k = j;
            p = d[j];
         }
      }
      if (k != i) {
         int j;
         d[k] = d[i];
         d[i] = p;
         for (j = 0; j < n; j++) {
            p = MatGet0(V,j,i);
            MatSet0(V,j,i, MatGet0(V,j,k));
            MatSet0(V,j,k, p);
         }
      }
   }
}

// Nonsymmetric reduction to Hessenberg form.

void orthes (MATRIX H, MATRIX V, MATRIX ort_)
{

   //  This is derived from the Algol procedures orthes and ortran,
   //  by Martin and Wilkinson, Handbook for Auto. Comp.,
   //  Vol.ii-Linear Algebra, and the corresponding
   //  Fortran subroutines in EISPACK.

   int  n    = MatRows(V);
   int  low  = 0;
   int  high = n-1;
   REAL* ort = MatData(ort_);
   int m, i;


   for (m = low+1; m <= high-1; m++) {

      // Scale column.

      REAL scale = 0.0;
      int i;
      for (i = m; i <= high; i++) {
         scale = scale + sml_fabs(MatGet0(H,i,m-1));
      }
      if (scale != 0.0) {

         // Compute Householder transformation.

         REAL h = 0.0, g;
         int i, j;
         for (i = high; i >= m; i--) {
            ort[i] = MatGet0(H,i,m-1)/scale;
            h += ort[i] * ort[i];
         }
         g = sml_sqrt(h);
         if (ort[m] > 0) {
            g = -g;
         }
         h = h - ort[m] * g;
         ort[m] = ort[m] - g;

         // Apply Householder similarity transformation
         // H = (I-u*u'/h)*H*(I-u*u'/h)

         for (j = m; j < n; j++) {
            REAL f = 0.0;
            int i;
            for (i = high; i >= m; i--) {
               f += ort[i]*MatGet0(H,i,j);
            }
            f = f/h;
            for (i = m; i <= high; i++) {
               MatSetME0(H,i,j, f*ort[i]);
            }
        }

        for (i = 0; i <= high; i++) {
            REAL f = 0.0;
            int j;
            for (j = high; j >= m; j--) {
               f += ort[j]*MatGet0(H,i,j);
            }
            f = f/h;
            for (j = m; j <= high; j++) {
               MatSetME0(H,i,j, f*ort[j]);
            }
         }
         ort[m] = scale*ort[m];
         MatSet0(H,m,m-1,  scale*g);
      }
   }

   // Accumulate transformations (Algol's ortran).

   for (i = 0; i < n; i++) {
      int j;
      for (j = 0; j < n; j++) {
         MatSet0(V,i,j,  (i == j ? 1.0 : 0.0));
      }
   }

   for (m = high-1; m >= low+1; m--) {
      if (MatGet0(H,m,m-1) != 0.0) {
         int i, j;
         for (i = m+1; i <= high; i++) {
            ort[i] = MatGet0(H,i,m-1);
         }
         for (j = m; j <= high; j++) {
            REAL g = 0.0;
            int i;
            for (i = m; i <= high; i++) {
               g += ort[i] * MatGet0(V,i,j);
            }
            // Double division avoids possible underflow
            g = (g / ort[m]) / MatGet0(H,m,m-1);
            for (i = m; i <= high; i++) {
               MatSetPE0(V,i,j, g * ort[i]);
            }
         }
      }
   }
}


// Complex scalar division.

static REAL cdivr, cdivi;     // result of division

static void cdiv(REAL xr, REAL xi, REAL yr, REAL yi) {
   REAL r,d;
   if (sml_fabs(yr) > sml_fabs(yi)) {
      r = yi/yr;
      d = yr + r*yi;
      cdivr = (xr + r*xi)/d;
      cdivi = (xi - r*xr)/d;
   } else {
      r = yr/yi;
      d = yi + r*yr;
      cdivr = (r*xr + xi)/d;
      cdivi = (r*xi - xr)/d;
   }
}


// Nonsymmetric reduction from Hessenberg to real Schur form.

static void hqr2 (MATRIX H, MATRIX V, MATRIX d_, MATRIX e_)
{

   //  This is derived from the Algol procedure hqr2,
   //  by Martin and Wilkinson, Handbook for Auto. Comp.,
   //  Vol.ii-Linear Algebra, and the corresponding
   //  Fortran subroutine in EISPACK.

   // Initialize

   int nn = MatRows(H);
   REAL* d = MatData(d_);
   REAL* e = MatData(e_);
   int n = nn-1;
   int low = 0;
   int high = nn-1;
   REAL eps  = REAL_EPSILON;
   REAL exshift = 0.0;
   REAL p=0,q=0,r=0,s=0,z=0,t,w,x,y;
   int i, j;
   int iter = 0;

   // Store roots isolated by balanc and compute matrix norm

   REAL norm = 0.0;
   for (i = 0; i < nn; i++) {
      int j;
      if ((i < low) || (i > high)) {
         d[i] = MatGet0(H,i,i);
         e[i] = 0.0;
      }
      for (j = Matmaxi(i-1,0); j < nn; j++) {
         norm = norm + sml_fabs(MatGet0(H,i,j));
      }
   }

   // Outer loop over eigenvalue index

   while (n >= low) {

      // Look for single small sub-diagonal element

      int l = n;
      while (l > low) {
         s = sml_fabs(MatGet0(H,l-1,l-1)) + sml_fabs(MatGet0(H,l,l));
         if (s == 0.0) {
            s = norm;
         }
         if (sml_fabs(MatGet0(H,l,l-1)) < eps * s) {
            break;
         }
         l--;
      }

      // Check for convergence
      // One root found

      if (l == n) {
         MatSet0(H,n,n, MatGet0(H,n,n) + exshift);
         d[n] = MatGet0(H,n,n);
         e[n] = 0.0;
         n--;
         iter = 0;

      // Two roots found

      } else if (l == n-1) {
         w = MatGet0(H,n,n-1) * MatGet0(H,n-1,n);
         p = (MatGet0(H,n-1,n-1) - MatGet0(H,n,n)) / 2.0;
         q = p * p + w;
         z = sml_sqrt(sml_fabs(q));
         MatSet0(H,n,n, MatGet0(H,n,n) + exshift);
         MatSet0(H,n-1,n-1, MatGet0(H,n-1,n-1) + exshift);
         x = MatGet0(H,n,n);

         // REAL pair

         if (q >= 0) {
            int j, i;
            if (p >= 0) {
               z = p + z;
            } else {
               z = p - z;
            }
            d[n-1] = x + z;
            d[n] = d[n-1];
            if (z != 0.0) {
               d[n] = x - w / z;
            }
            e[n-1] = 0.0;
            e[n] = 0.0;
            x = MatGet0(H,n,n-1);
            s = sml_fabs(x) + sml_fabs(z);
            p = x / s;
            q = z / s;
            r = sml_sqrt(p * p+q * q);
            p = p / r;
            q = q / r;

            // Row modification

            for (j = n-1; j < nn; j++) {
               z = MatGet0(H,n-1,j);
               MatSet0(H,n-1,j,  q * z + p * MatGet0(H,n,j));
               MatSet0(H,n,j,    q * MatGet0(H,n,j) - p * z);
            }

            // Column modification

            for (i = 0; i <= n; i++) {
               z = MatGet0(H,i,n-1);
               MatSet0(H,i,n-1,  q * z + p * MatGet0(H,i,n));
               MatSet0(H,i,n,    q * MatGet0(H,i,n) - p * z);
            }

            // Accumulate transformations

            for (i = low; i <= high; i++) {
               z = MatGet0(V,i,n-1);
               MatSet0(V,i,n-1, q * z + p * MatGet0(V,i,n));
               MatSet0(V,i,n,   q * MatGet0(V,i,n) - p * z);
            }

         // Complex pair

         } else {
            d[n-1] = x + p;
            d[n] = x + p;
            e[n-1] = z;
            e[n] = -z;
         }
         n = n - 2;
         iter = 0;

      // No convergence yet

      } else {

         int i, k, m;
         // Form shift

         x = MatGet0(H,n,n);
         y = 0.0;
         w = 0.0;
         if (l < n) {
            y = MatGet0(H,n-1,n-1);
            w = MatGet0(H,n,n-1) * MatGet0(H,n-1,n);
         }

         // Wilkinson's original ad hoc shift

         if (iter == 10) {
            int i;
            exshift += x;
            for (i = low; i <= n; i++) {
               MatSetME0(H,i,i, x);
            }
            s = sml_fabs(MatGet0(H,n,n-1)) + sml_fabs(MatGet0(H,n-1,n-2));
            x = y = 0.75 * s;
            w = -0.4375 * s * s;
         }

         // MATLAB's new ad hoc shift

         if (iter == 30) {
             s = (y - x) / 2.0;
             s = s * s + w;
             if (s > 0) {
                 int i;
                 s = sml_sqrt(s);
                 if (y < x) {
                    s = -s;
                 }
                 s = x - w / ((y - x) / 2.0 + s);
                 for (i = low; i <= n; i++) {
                    MatSetME0(H,i,i,s);
                 }
                 exshift += s;
                 x = y = w = 0.964;
             }
         }

         iter = iter + 1;   // (Could check iteration count here.)

         // Look for two consecutive small sub-diagonal elements

         m = n-2;
         while (m >= l) {
            z = MatGet0(H,m,m);
            r = x - z;
            s = y - z;
            p = (r * s - w) / MatGet0(H,m+1,m) + MatGet0(H,m,m+1);
            q = MatGet0(H,m+1,m+1) - z - r - s;
            r = MatGet0(H,m+2,m+1);
            s = sml_fabs(p) + sml_fabs(q) + sml_fabs(r);
            p = p / s;
            q = q / s;
            r = r / s;
            if (m == l) {
               break;
            }
            if (sml_fabs(MatGet0(H,m,m-1)) * (sml_fabs(q) + sml_fabs(r)) <
               eps * (sml_fabs(p) * (sml_fabs(MatGet0(H,m-1,m-1)) + sml_fabs(z) +
               sml_fabs(MatGet0(H,m+1,m+1))))) {
                  break;
            }
            m--;
         }

         for (i = m+2; i <= n; i++) {
            MatSet0(H,i,i-2, 0.0);
            if (i > m+2) {
               MatSet0(H,i,i-3,  0.0);
            }
         }

         // Double QR step involving rows l:n and columns m:n

         for (k = m; k <= n-1; k++) {
            int notlast = (k != n-1);
            if (k != m) {
               p = MatGet0(H,k,k-1);
               q = MatGet0(H,k+1,k-1);
               r = (notlast ? MatGet0(H,k+2,k-1) : 0.0);
               x = sml_fabs(p) + sml_fabs(q) + sml_fabs(r);
               if (x != 0.0) {
                  p = p / x;
                  q = q / x;
                  r = r / x;
               }
            }
            if (x == 0.0) {
               break;
            }
            s = sml_sqrt(p * p + q * q + r * r);
            if (p < 0) {
               s = -s;
            }
            if (s != 0) {
               int j, i;
               if (k != m) {
                  MatSet0(H,k,k-1,  -s * x);
               } else if (l != m) {
                  MatSet0(H,k,k-1,  -MatGet0(H,k,k-1));
               }
               p = p + s;
               x = p / s;
               y = q / s;
               z = r / s;
               q = q / p;
               r = r / p;

               // Row modification

               for (j = k; j < nn; j++) {
                  p = MatGet0(H,k,j) + q * MatGet0(H,k+1,j);
                  if (notlast) {
                     p = p + r * MatGet0(H,k+2,j);
                     MatSetME0(H,k+2,j,  p * z);
                  }
                  MatSetME0(H,k,j, p * x);
                  MatSetME0(H,k+1,j, p * y);
               }

               // Column modification

               for (i = 0; i <= Matmini(n,k+3); i++) {
                  p = x * MatGet0(H,i,k) + y * MatGet0(H,i,k+1);
                  if (notlast) {
                     p = p + z * MatGet0(H,i,k+2);
                     MatSetME0(H,i,k+2, p * r);
                  }
                  MatSetME0(H,i,k,   p);
                  MatSetME0(H,i,k+1, p * q);
               }

               // Accumulate transformations

               for (i = low; i <= high; i++) {
                  p = x * MatGet0(V,i,k) + y * MatGet0(V,i,k+1);
                  if (notlast) {
                     p = p + z * MatGet0(V,i,k+2);
                     MatSetME0(V,i,k+2, p * r);
                  }
                  MatSetME0(V,i,k  , p);
                  MatSetME0(V,i,k+1, p * q);
               }
            }  // (s != 0)
         }  // k loop
      }  // check convergence
   }  // while (n >= low)

   // Backsubstitute to find vectors of upper triangular form

   if (norm == 0.0) {
      return;
   }

   for (n = nn-1; n >= 0; n--) {
      p = d[n];
      q = e[n];

      // REAL vector

      if (q == 0) {
         int l = n;
         int i;
         MatSet0(H,n,n, 1.0);
         for (i = n-1; i >= 0; i--) {
            int j;
            w = MatGet0(H,i,i) - p;
            r = 0.0;
            for (j = l; j <= n; j++) {
               r = r + MatGet0(H,i,j) * MatGet0(H,j,n);
            }
            if (e[i] < 0.0) {
               z = w;
               s = r;
            } else {
               l = i;
               if (e[i] == 0.0) {
                  if (w != 0.0) {
                     MatSet0(H,i,n,  -r / w);
                  } else {
                     MatSet0(H,i,n,  -r / (eps * norm));
                  }

               // Solve real equations

               } else {
                  x = MatGet0(H,i,i+1);
                  y = MatGet0(H,i+1,i);
                  q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
                  t = (x * s - z * r) / q;
                  MatSet0(H,i,n,  t);
                  if (sml_fabs(x) > sml_fabs(z)) {
                     MatSet0(H,i+1,n, (-r - w * t) / x);
                  } else {
                     MatSet0(H,i+1,n, (-s - y * t) / z);
                  }
               }

               // Overflow control

               t = sml_fabs(MatGet0(H,i,n));
               if ((eps * t) * t > 1) {
                  int j;
                  for (j = i; j <= n; j++) {
                     MatSet0(H,j,n,  MatGet0(H,j,n) / t);
                  }
               }
            }
         }

      // Complex vector

      } else if (q < 0) {
         int l = n-1;
         int i;

         // Last vector component imaginary so matrix is triangular

         if (sml_fabs(MatGet0(H,n,n-1)) > sml_fabs(MatGet0(H,n-1,n))) {
            MatSet0(H,n-1,n-1,  q / MatGet0(H,n,n-1));
            MatSet0(H,n-1,n,   -(MatGet0(H,n,n) - p) / MatGet0(H,n,n-1));
         } else {
            cdiv(0.0,-MatGet0(H,n-1,n),MatGet0(H,n-1,n-1)-p,q);
            MatSet0(H,n-1,n-1,   cdivr);
            MatSet0(H,n-1,n,   cdivi);
         }
         MatSet0(H,n,n-1,   0.0);
         MatSet0(H,n,n,   1.0);
         for (i = n-2; i >= 0; i--) {
            REAL ra,sa,vr,vi;
            int j;
            ra = 0.0;
            sa = 0.0;
            for (j = l; j <= n; j++) {
               ra = ra + MatGet0(H,i,j) * MatGet0(H,j,n-1);
               sa = sa + MatGet0(H,i,j) * MatGet0(H,j,n);
            }
            w = MatGet0(H,i,i) - p;

            if (e[i] < 0.0) {
               z = w;
               r = ra;
               s = sa;
            } else {
               l = i;
               if (e[i] == 0) {
                  cdiv(-ra,-sa,w,q);
                  MatSet0(H,i,n-1, cdivr);
                  MatSet0(H,i,n,   cdivi);
               } else {

                  // Solve complex equations

                  x = MatGet0(H,i,i+1);
                  y = MatGet0(H,i+1,i);
                  vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
                  vi = (d[i] - p) * 2.0 * q;
                  if ((vr == 0.0) && (vi == 0.0)) {
                     vr = eps * norm * (sml_fabs(w) + sml_fabs(q) +
                     sml_fabs(x) + sml_fabs(y) + sml_fabs(z));
                  }
                  cdiv(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi);
                  MatSet0(H,i,n-1, cdivr);
                  MatSet0(H,i,n,   cdivi);
                  if (sml_fabs(x) > (sml_fabs(z) + sml_fabs(q))) {
                     MatSet0(H,i+1,n-1, (-ra - w * MatGet0(H,i,n-1) + q * MatGet0(H,i,n)) / x);
                     MatSet0(H,i+1,n,   (-sa - w * MatGet0(H,i,n) - q * MatGet0(H,i,n-1)) / x);
                  } else {
                     cdiv(-r-y*MatGet0(H,i,n-1),-s-y*MatGet0(H,i,n),z,q);
                     MatSet0(H,i+1,n-1, cdivr);
                     MatSet0(H,i+1,n,   cdivi);
                  }
               }

               // Overflow control

               t = Matmaxr(sml_fabs(MatGet0(H,i,n-1)),sml_fabs(MatGet0(H,i,n)));
               if ((eps * t) * t > 1) {
                  int j;
                  for (j = i; j <= n; j++) {
                     MatSetDE0(H,j,n-1, t);
                     MatSetDE0(H,j,n,   t);
                  }
               }
            }
         }
      }
   }

   // Vectors of isolated roots

   for (i = 0; i < nn; i++) {
      if (i < low || i > high) {
         int j;
         for (j = i; j < nn; j++) {
            MatSet0(V,i,j,  MatGet0(H,i,j));
         }
      }
   }

   // Back transformation to get eigenvectors of original matrix

   for (j = nn-1; j >= low; j--) {
      int i;
      for (i = low; i <= high; i++) {
         int k;
         z = 0.0;
         for (k = low; k <= Matmini(j,high); k++) {
            z = z + MatGet0(V,i,k) * MatGet0(H,k,j);
         }
         MatSet0(V,i,j,  z);
      }
   }
}

// public:

int  MatEigen (MATRIX V, MATRIX d_, MATRIX e_, MATRIX A)
/*
Input: A is a  Square real (non-complex) matrix
V  <- the eigenvector matrix
d_ <- the real parts of the eigenvalues
e_ <- the imaginary parts of the eigenvalues
Returns 0 if no errors.

Check for symmetry, then construct the eigenvalue decomposition
*/
{

   int n  = MatCols(A); // Row and column dimension (square matrix)
   int issymmetric;     // boolean

   if (!MatIsSquare(A)) {
     MatErrThrow("MatEigen: Not a square matrix.");
     return 1;
  }

  if (V == d_ || V == e_ || d_== e_)  {
    MatErrThrow("MatEigen: arguments must be distinct.");
    return 2;
  }

  MatReDim(V,  n, n);  // eigenvectors.
  MatReDim(d_, n, 1);  // real parts of the eigenvalues
  MatReDim(e_, n, 1);  // imaginary parts of the eigenvalues

  issymmetric = MatIsSymmetric(A);
  if (issymmetric)  {
     MatCopy(V, A);

     // Tridiagonalize.
     tred2(V, d_, e_);

     // Diagonalize.
     tql2(V, d_, e_);

  }
  else {
     MATRIX H = MatDim(0, 0);    // storage of nonsymmetric Hessenberg form.
     MATRIX ort_ = MatDim(n, 1); // storage for nonsymmetric algorithm.
     MatCopy(H, A);

     // Reduce to Hessenberg form.
     orthes(H, V, ort_);

     // Reduce Hessenberg to real Schur form.
     hqr2(H, V, d_, e_);
  }

  return 0;
}


// don't need getD for now, may be activated if needed
#if 0

/**
 Computes the block diagonal eigenvalue matrix.
    If the original matrix A is not symmetric, then the eigenvalue
 matrix D is block diagonal with the real eigenvalues in 1-by-1
 blocks and any complex eigenvalues,
    a + i*b, in 2-by-2 blocks, [a, b; -b, a].  That is, if the complex
    eigenvalues look like
<pre>

          u + iv     .        .          .      .    .
            .      u - iv     .          .      .    .
            .        .      a + ib       .      .    .
            .        .        .        a - ib   .    .
            .        .        .          .      x    .
            .        .        .          .      .    y
</pre>
        then D looks like
<pre>

            u        v        .          .      .    .
           -v        u        .          .      .    .
            .        .        a          b      .    .
            .        .       -b          a      .    .
            .        .        .          .      x    .
            .        .        .          .      .    y
</pre>
    This keeps V a real matrix in both symmetric and non-symmetric
    cases, and A*V = V*D.

 @param D: upon return, the matrix is filled with the block diagonal
 eigenvalue matrix.

*/
void getD (MATRIX D, MATRIX d_, MATRIX e_)
{
   REAL*  d  = MatData(d_);
   REAL*  e  = MatData(e_);
   int n = MatRows(d_);
   int i, j;

   MatReDim(D, n, n);
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++)   MatSet0(D,i,j, 0.0);
      MatGet0(D,i,i) = d[i];
      if (e[i] > 0)   MatSet0(D,i,i+1,  e[i]);
        else if (e[i] < 0)   MatSet0(D,i,i-1,  e[i]);
   }
}

#endif
