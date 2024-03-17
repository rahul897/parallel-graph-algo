# include <cstdlib>
#include <fstream>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
#include<sys/time.h>
# include "cal.cpp"
using namespace std;

int dsvdc ( double a[], int lda, int m, int n, double s[], double e[],
  double u[], int ldu, double v[], int ldv, double work[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    DSVDC computes the singular value decomposition of a real rectangular matrix.
//
//  Discussion:
//
//    This routine reduces an M by N matrix A to diagonal form by orthogonal
//    transformations U and V.  The diagonal elements S(I) are the singular
//    values of A.  The columns of U are the corresponding left singular
//    vectors, and the columns of V the right singular vectors.
//
//    The form of the singular value decomposition is then
//
//      A(MxN) = U(MxM) * S(MxN) * V(NxN)'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 May 2007
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*N].  On input, the M by N matrix whose
//    singular value decomposition is to be computed.  On output, the matrix
//    has been destroyed.  Depending on the user's requests, the matrix may
//    contain other useful information.
//
//    Input, int LDA, the leading dimension of the array A.
//    LDA must be at least M.
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix A.
//
//    Output, double S[MM], where MM = min(M+1,N).  The first
//    min(M,N) entries of S contain the singular values of A arranged in
//    descending order of magnitude.
//
//    Output, double E[MM], where MM = min(M+1,N), ordinarily contains zeros.
//    However see the discussion of INFO for exceptions.
//
//    Output, double U[LDU*K].  If JOBA = 1 then K = M;
//    if 2 <= JOBA, then K = min(M,N).  U contains the M by M matrix of left singular
//    vectors.  U is not referenced if JOBA = 0.  If M <= N or if JOBA = 2, then
//    U may be identified with A in the subroutine call.
//
//    Input, int LDU, the leading dimension of the array U.
//    LDU must be at least M.
//
//    Output, double V[LDV*N], the N by N matrix of right singular vectors.
//    V is not referenced if JOB is 0.  If N <= M, then V may be identified
//    with A in the subroutine call.
//
//    Input, int LDV, the leading dimension of the array V.
//    LDV must be at least N.
//
//    Workspace, double WORK[M].
//
//    Input, int JOB, controls the computation of the singular
//    vectors.  It has the decimal expansion AB with the following meaning:
//      A =  0, do not compute the left singular vectors.
//      A =  1, return the M left singular vectors in U.
//      A >= 2, return the first min(M,N) singular vectors in U.
//      B =  0, do not compute the right singular vectors.
//      B =  1, return the right singular vectors in V.
//
//    Output, int *DSVDC, status indicator INFO.
//    The singular values (and their corresponding singular vectors)
//    S(*INFO+1), S(*INFO+2),...,S(MN) are correct.  Here MN = min ( M, N ).
//    Thus if *INFO is 0, all the singular values and their vectors are
//    correct.  In any event, the matrix B = U' * A * V is the bidiagonal
//    matrix with the elements of S on its diagonal and the elements of E on
//    its superdiagonal.  Thus the singular values of A and B are the same.
//
{
  double b;
  double c;
  double cs;
  double el;
  double emm1;
  double f;
  double g;
  int i;
  int info;
  int iter;
  int j;
  int jobu;
  int k;
  int kase;
  int kk;
  int l;
  int ll;
  int lls;
  int ls;
  int lu;
  int maxit = 30;
  int mm;
  int mm1;
  int mn;
  int mp1;
  int nct;
  int nctp1;
  int ncu;
  int nrt;
  int nrtp1;
  double scale;
  double shift;
  double sl;
  double sm;
  double smm1;
  double sn;
  double t;
  double t1;
  double test;
  bool wantu;
  bool wantv;
  double ztest;
//
//  Determine what is to be computed.
//
  info = 0;
  wantu = false;
  wantv = false;
  jobu = ( job % 100 ) / 10;

  if ( 1 < jobu )
  {
    ncu = i4_min ( m, n );
  }
  else
  {
    ncu = m;
  }

  if ( jobu != 0 )
  {
    wantu = true;
  }

  if ( ( job % 10 ) != 0 )
  {
    wantv = true;
  }
//
//  Reduce A to bidiagonal form, storing the diagonal elements
//  in S and the super-diagonal elements in E.
//
  nct = i4_min ( m-1, n );
  nrt = i4_max ( 0, i4_min ( m, n-2 ) );
  lu = i4_max ( nct, nrt );

  for ( l = 1; l <= lu; l++ )
  {
//
//  Compute the transformation for the L-th column and
//  place the L-th diagonal in S(L).
//
    if ( l <= nct )
    {
      s[l-1] = dnrm2 ( m-l+1, a+l-1+(l-1)*lda, 1 );

      if ( s[l-1] != 0.0 )
      {
        if ( a[l-1+(l-1)*lda] != 0.0 )
        {
          s[l-1] = r8_sign ( a[l-1+(l-1)*lda] ) * r8_abs ( s[l-1] );
        }
        dscal ( m-l+1, 1.0 / s[l-1], a+l-1+(l-1)*lda, 1 );
        a[l-1+(l-1)*lda] = 1.0 + a[l-1+(l-1)*lda];
      }
      s[l-1] = -s[l-1];
    }

    for ( j = l+1; j <= n; j++ )
    {
//
//  Apply the transformation.
//
      if ( l <= nct && s[l-1] != 0.0 )
      {
        t = - ddot ( m-l+1, a+l-1+(l-1)*lda, 1, a+l-1+(j-1)*lda, 1 )
          / a[l-1+(l-1)*lda];
        daxpy ( m-l+1, t, a+l-1+(l-1)*lda, 1, a+l-1+(j-1)*lda, 1 );
      }
//
//  Place the L-th row of A into E for the
//  subsequent calculation of the row transformation.
//
      e[j-1] = a[l-1+(j-1)*lda];
    }
//
//  Place the transformation in U for subsequent back multiplication.
//
    if ( wantu && l <= nct )
    {
      for ( i = l; i <= m; i++ )
      {
        u[i-1+(l-1)*ldu] = a[i-1+(l-1)*lda];
      }
    }

    if ( l <= nrt )
    {
//
//  Compute the L-th row transformation and place the
//  L-th superdiagonal in E(L).
//
      e[l-1] = dnrm2 ( n-l, e+l, 1 );

      if ( e[l-1] != 0.0 )
      {
        if ( e[l] != 0.0 )
        {
          e[l-1] = r8_sign ( e[l] ) * r8_abs ( e[l-1] );
        }
        dscal ( n-l, 1.0 / e[l-1], e+l, 1 );
        e[l] = 1.0 + e[l];
      }

      e[l-1] = -e[l-1];
//
//  Apply the transformation.
//
      if ( l+1 <= m && e[l-1] != 0.0 )
      {
        for ( j = l+1; j <= m; j++ )
        {
          work[j-1] = 0.0;
        }

        for ( j = l+1; j <= n; j++ )
        {
          daxpy ( m-l, e[j-1], a+l+(j-1)*lda, 1, work+l, 1 );
        }

        for ( j = l+1; j <= n; j++ )
        {
          daxpy ( m-l, -e[j-1]/e[l], work+l, 1, a+l+(j-1)*lda, 1 );
        }
      }
//
//  Place the transformation in V for subsequent back multiplication.
//
      if ( wantv )
      {
        for ( j = l+1; j <= n; j++ )
        {
          v[j-1+(l-1)*ldv] = e[j-1];
        }
      }
    }
  }
//
//  Set up the final bidiagonal matrix of order MN.
//
  mn = i4_min ( m + 1, n );
  nctp1 = nct + 1;
  nrtp1 = nrt + 1;

  if ( nct < n )
  {
    s[nctp1-1] = a[nctp1-1+(nctp1-1)*lda];
  }

  if ( m < mn )
  {
    s[mn-1] = 0.0;
  }

  if ( nrtp1 < mn )
  {
    e[nrtp1-1] = a[nrtp1-1+(mn-1)*lda];
  }

  e[mn-1] = 0.0;
//
//  If required, generate U.
//
  if ( wantu )
  {
    for ( i = 1; i <= m; i++ )
    {
      for ( j = nctp1; j <= ncu; j++ )
      {
        u[(i-1)+(j-1)*ldu] = 0.0;
      }
    }

    for ( j = nctp1; j <= ncu; j++ )
    {
      u[j-1+(j-1)*ldu] = 1.0;
    }

    for ( ll = 1; ll <= nct; ll++ )
    {
      l = nct - ll + 1;

      if ( s[l-1] != 0.0 )
      {
        for ( j = l+1; j <= ncu; j++ )
        {
          t = - ddot ( m-l+1, u+(l-1)+(l-1)*ldu, 1, u+(l-1)+(j-1)*ldu, 1 )
            / u[l-1+(l-1)*ldu];
          daxpy ( m-l+1, t, u+(l-1)+(l-1)*ldu, 1, u+(l-1)+(j-1)*ldu, 1 );
        }

        dscal ( m-l+1, -1.0, u+(l-1)+(l-1)*ldu, 1 );
        u[l-1+(l-1)*ldu] = 1.0 + u[l-1+(l-1)*ldu];
        for ( i = 1; i <= l-1; i++ )
        {
          u[i-1+(l-1)*ldu] = 0.0;
        }
      }
      else
      {
        for ( i = 1; i <= m; i++ )
        {
          u[i-1+(l-1)*ldu] = 0.0;
        }
        u[l-1+(l-1)*ldu] = 1.0;
      }
    }
  }
//
//  If it is required, generate V.
//
  if ( wantv )
  {
    for ( ll = 1; ll <= n; ll++ )
    {
      l = n - ll + 1;

      if ( l <= nrt && e[l-1] != 0.0 )
      {
        for ( j = l+1; j <= n; j++ )
        {
          t = - ddot ( n-l, v+l+(l-1)*ldv, 1, v+l+(j-1)*ldv, 1 )
            / v[l+(l-1)*ldv];
          daxpy ( n-l, t, v+l+(l-1)*ldv, 1, v+l+(j-1)*ldv, 1 );
        }

      }
      for ( i = 1; i <= n; i++ )
      {
        v[i-1+(l-1)*ldv] = 0.0;
      }
      v[l-1+(l-1)*ldv] = 1.0;
    }
  }
//
//  Main iteration loop for the singular values.
//
  mm = mn;
  iter = 0;

  while ( 0 < mn )
  {
//
//  If too many iterations have been performed, set flag and return.
//
    if ( maxit <= iter )
    {
      info = mn;
      return info;
    }
//
//  This section of the program inspects for
//  negligible elements in the S and E arrays.
//
//  On completion the variables KASE and L are set as follows:
//
//  KASE = 1     if S(MN) and E(L-1) are negligible and L < MN
//  KASE = 2     if S(L) is negligible and L < MN
//  KASE = 3     if E(L-1) is negligible, L < MN, and
//               S(L), ..., S(MN) are not negligible (QR step).
//  KASE = 4     if E(MN-1) is negligible (convergence).
//
    for ( ll = 1; ll <= mn; ll++ )
    {
      l = mn - ll;

      if ( l == 0 )
      {
        break;
      }

      test = r8_abs ( s[l-1] ) + r8_abs ( s[l] );
      ztest = test + r8_abs ( e[l-1] );

      if ( ztest == test )
      {
        e[l-1] = 0.0;
        break;
      }
    }

    if ( l == mn - 1 )
    {
      kase = 4;
    }
    else
    {
      mp1 = mn + 1;

      for ( lls = l+1; lls <= mn+1; lls++ )
      {
        ls = mn - lls + l + 1;

        if ( ls == l )
        {
          break;
        }

        test = 0.0;
        if ( ls != mn )
        {
          test = test + r8_abs ( e[ls-1] );
        }

        if ( ls != l + 1 )
        {
          test = test + r8_abs ( e[ls-2] );
        }

        ztest = test + r8_abs ( s[ls-1] );

        if ( ztest == test )
        {
          s[ls-1] = 0.0;
          break;
        }

      }

      if ( ls == l )
      {
        kase = 3;
      }
      else if ( ls == mn )
      {
        kase = 1;
      }
      else
      {
        kase = 2;
        l = ls;
      }
    }

    l = l + 1;
//
//  Deflate negligible S(MN).
//
    if ( kase == 1 )
    {
      mm1 = mn - 1;
      f = e[mn-2];
      e[mn-2] = 0.0;

      for ( kk = 1; kk <= mm1; kk++ )
      {
        k = mm1 - kk + l;
        t1 = s[k-1];
        drotg ( &t1, &f, &cs, &sn );
        s[k-1] = t1;

        if ( k != l )
        {
          f = -sn * e[k-2];
          e[k-2] = cs * e[k-2];
        }

        if ( wantv )
        {
          drot ( n, v+0+(k-1)*ldv, 1, v+0+(mn-1)*ldv, 1, cs, sn );
        }
      }
    }
//
//  Split at negligible S(L).
//
    else if ( kase == 2 )
    {
      f = e[l-2];
      e[l-2] = 0.0;

      for ( k = l; k <= mn; k++ )
      {
        t1 = s[k-1];
        drotg ( &t1, &f, &cs, &sn );
        s[k-1] = t1;
        f = - sn * e[k-1];
        e[k-1] = cs * e[k-1];
        if ( wantu )
        {
          drot ( m, u+0+(k-1)*ldu, 1, u+0+(l-2)*ldu, 1, cs, sn );
        }
      }
    }
//
//  Perform one QR step.
//
    else if ( kase == 3 )
    {
//
//  Calculate the shift.
//
      scale = r8_max ( r8_abs ( s[mn-1] ),
              r8_max ( r8_abs ( s[mn-2] ),
              r8_max ( r8_abs ( e[mn-2] ),
              r8_max ( r8_abs ( s[l-1] ), r8_abs ( e[l-1] ) ) ) ) );

      sm = s[mn-1] / scale;
      smm1 = s[mn-2] / scale;
      emm1 = e[mn-2] / scale;
      sl = s[l-1] / scale;
      el = e[l-1] / scale;
      b = ( ( smm1 + sm ) * ( smm1 - sm ) + emm1 * emm1 ) / 2.0;
      c = ( sm * emm1 ) * ( sm * emm1 );
      shift = 0.0;

      if ( b != 0.0 || c != 0.0 )
      {
        shift = sqrt ( b * b + c );
        if ( b < 0.0 )
        {
          shift = -shift;
        }
        shift = c / ( b + shift );
      }

      f = ( sl + sm ) * ( sl - sm ) - shift;
      g = sl * el;
//
//  Chase zeros.
//
      mm1 = mn - 1;

      for ( k = l; k <= mm1; k++ )
      {
        drotg ( &f, &g, &cs, &sn );

        if ( k != l )
        {
          e[k-2] = f;
        }

        f = cs * s[k-1] + sn * e[k-1];
        e[k-1] = cs * e[k-1] - sn * s[k-1];
        g = sn * s[k];
        s[k] = cs * s[k];

        if ( wantv )
        {
          drot ( n, v+0+(k-1)*ldv, 1, v+0+k*ldv, 1, cs, sn );
        }

        drotg ( &f, &g, &cs, &sn );
        s[k-1] = f;
        f = cs * e[k-1] + sn * s[k];
        s[k] = -sn * e[k-1] + cs * s[k];
        g = sn * e[k];
        e[k] = cs * e[k];

        if ( wantu && k < m )
        {
          drot ( m, u+0+(k-1)*ldu, 1, u+0+k*ldu, 1, cs, sn );
        }
      }
      e[mn-2] = f;
      iter = iter + 1;
    }
//
//  Convergence.
//
    else if ( kase == 4 )
    {
//
//  Make the singular value nonnegative.
//
      if ( s[l-1] < 0.0 )
      {
        s[l-1] = -s[l-1];
        if ( wantv )
        {
          dscal ( n, -1.0, v+0+(l-1)*ldv, 1 );
        }
      }
//
//  Order the singular value.
//
      for ( ; ; )
      {
        if ( l == mm )
        {
          break;
        }

        if ( s[l] <= s[l-1] )
        {
          break;
        }

        t = s[l-1];
        s[l-1] = s[l];
        s[l] = t;

        if ( wantv && l < n )
        {
          dswap ( n, v+0+(l-1)*ldv, 1, v+0+l*ldv, 1 );
        }

        if ( wantu && l < m )
        {
          dswap ( m, u+0+(l-1)*ldu, 1, u+0+l*ldu, 1 );
        }

        l = l + 1;
      }
      iter = 0;
      mn = mn - 1;
    }
  }

  return info;
}

void r8mat_svd_linpack ( int m, int n, double a[], double u[], double s[],
  double v[] )
{
  double *a_copy;
  double *e;
  int i;
  int info;
  int j;
  int lda;
  int ldu;
  int ldv;
  int job;
  int lwork;
  double *sdiag;
  double *work;
//
//  The correct size of E and SDIAG is min ( m+1, n).
//
  a_copy = new double[m*n];
  e = new double[m+n];
  sdiag = new double[m+n];
  work = new double[m];
//
//  Compute the eigenvalues and eigenvectors.
//
  job = 11;
  lda = m;
  ldu = m;
  ldv = n;
//
//  The input matrix is destroyed by the routine.  Since we need to keep
//  it around, we only pass a copy to the routine.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a_copy[i+j*m] = a[i+j*m];
    }
  }
  info = dsvdc ( a_copy, lda, m, n, sdiag, e, u, ldu, v, ldv, work, job );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "R8MAT_SVD_LINPACK - Failure!\n";
    cout << "  The SVD could not be calculated.\n";
    cout << "  LINPACK routine DSVDC returned a nonzero\n";
    cout << "  value of the error flag, INFO = " << info << "\n";
    return;
  }
//
//  Make the MxN matrix S from the diagonal values in SDIAG.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( i == j )
      {
        s[i+j*m] = sdiag[i];
      }
      else
      {
        s[i+j*m] = 0.0;
      }
    }
  }
//
//  Note that we do NOT need to transpose the V that comes out of LINPACK!
//
  delete [] a_copy;
  delete [] e;
  delete [] sdiag;
  delete [] work;

  return;
}

int getsvd(int m,int n)
{
  int i,j;
  double *u,*s,*a,*v;
  //cin>>m>>n;
 struct timeval it;
  gettimeofday(&it,NULL);
  unsigned int l=it.tv_usec;
  srand(l);
  fstream fout;
  fout.open("svdf.txt",ios::in|ios::out|ios::trunc);
  u = new double[m*m];
  s = new double[m*n];
  a = new double[m*n];
  v = new double[n*n];
  //
  //  Generate the matrix.
  //
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+m*j] = ((double)rand() / (double)(RAND_MAX));
    }
  }
  cout<<"-----------H-------------\n";
  for (int i = 0; i < m; ++i,cout<<endl)
  {
    for (int j = 0; j < n; j++)
    {
        cout<<a[i+m*j]<<" ";
    }
  }


  r8mat_svd_linpack ( m, n, a, u, s, v );
 // cout<<"-----------U-------------\n";
  for (int i = 0; i < m; ++i,fout<<endl)
  {
    for (int j = 0; j < m; j++)
    {
        fout<<-1*u[i+m*j]<<" ";
    }
  }
  fout<<"\n";
  //cout<<"-----------V-------------\n";
  for (int i = 0; i < n; ++i,fout<<endl)
  {
    for (int j = 0; j < n; j++)
    {
        fout<<-1*v[i+m*j]<<" ";
    }
  }
 // cout<<"-----------sigma-------------\n";
  fout<<"\n";

  for (int i = 0; i < m; ++i,fout<<endl)
  {
    for (int j = 0; j < n; j++)
    {
        fout<<s[i+m*j]<<" ";
    }
  }
  fout.close();
  return 0;
}
