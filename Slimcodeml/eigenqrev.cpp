//CHHS eigenqrev.cpp
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defines.h"
#include "global.h"
#include "functions.h"
#include "lapack.h"

#define USE_LAPACK

#ifndef USE_LAPACK
//CHHS removed old code
#else

#include "lapack.h"

//
//  Using LAPACK DSYEVR driver routine compute the eigenvalues and eigenvector of the symmetric input matrix U[n*n]
//  Reorders the output values so they are ordered as the ones computed by the original eigenRealSym() routine.
//
static void inline
eigenRealSym (double *U, int n, double *R, double ignored[])
{
  const int N = 61;
  const int I0 = 0;
  const double D0 = 0.0;
  int m;
  int info;
  int isuppz[2 * N];
  const int lwork = 26 * N;	//CHHS size should be determined 
  double work[lwork];
  const int liwork = 10 * N;	//CHHS size should be determined
  int iwork[liwork];
  double U1[N * N];

  // Compute eigenvalues and eigenvectors for the full symmetric matrix
  dsyevr_ ("V", "A", "U", &n, U, &n, &D0, &D0, &I0, &I0, &D0, &m, R, U1, &n,
	   isuppz, work, &lwork, iwork, &liwork, &info);

  // Reorder eigenvalues
  int i;
  double t;
  int mid = n / 2;
  for (i = 0; i < mid; ++i)
    {
      t = R[i];
      R[i] = R[n - 1 - i];
      R[n - 1 - i] = t;
    }

  // Reorder eigenvectors
  int c, r;
  for (c = 0; c < n; ++c)
    {
      for (r = 0; r < n; ++r)
	{
	  U[r * n + c] = U1[(n - 1 - c) * n + r];
	}
    }
}
#endif

/* eigen solution for real symmetric matrix */
void
eigenQREV (double Q[], double pi[], int n,
	   double Root[], double U[], double V[], double spacesqrtpi[])
{
  /*
     This finds the eigen solution of the rate matrix Q for a time-reversible
     Markov process, using the algorithm for a real symmetric matrix.
     Rate matrix Q = S * diag{pi} = U * diag{Root} * V,
     where S is symmetrical, all elements of pi are positive, and U*V = I.
     space[n] is for storing sqrt(pi).

     [U 0] [Q_0 0] [U^-1 0]    [Root  0]
     [0 I] [0   0] [0    I]  = [0     0]

     Ziheng Yang, 25 December 2001 (ref is CME/eigenQ.pdf)
   */
  int i, j, inew, jnew, nnew, status;
  double *pi_sqrt = spacesqrtpi, small = 1e-100;
  int ONE = 1;			//CHHS For BLAS / LAPACK
  for (j = 0, nnew = 0; j < n; j++)
    if (pi[j] > small)
      {
	pi_sqrt[nnew++] = sqrt (pi[j]);
      }

  /* store in U the symmetrical matrix S = sqrt(D) * Q * sqrt(-D) */

  //CHHS VERY important: we are calling BLAS (not CBLAS), therefore routines expect all matrices in column major order, although in C we use row major order! This refers to all 6 steps!
  if (nnew == n)		//CHHS Step 1
    {				//CHHS Here we compute PI^(1/2) * Q* PI^(-1/2); not rewritten in BLAS because i) here performance does not matter, ii) BLAS solution needs to dscal from both sides which is less efficient
      for (i = 0; i < n; i++)
	for (j = 0, U[i * n + i] = Q[i * n + i]; j < i; j++)
	  U[i * n + j] = U[j * n + i] =
	    (Q[i * n + j] * pi_sqrt[i] / pi_sqrt[j]);
      eigenRealSym (U, n, Root, V);	//CHHS Step 2, matrix decomposition to follow
    }
  else
    {				//CHHS Only necessary if one or more pi values 0
      printf ("Zero base frequencies are currently not supported, exiting!\n");	//CHHS else-part currently not working (CodonFreq=0 works)
      exit (1);
      for (i = 0, inew = 0; i < n; i++)
	{
	  if (pi[i] > small)
	    {
	      for (j = 0, jnew = 0; j < i; j++)
		if (pi[j] > small)
		  {
		    U[inew * nnew + jnew] = U[jnew * nnew + inew]
		      = Q[i * n + j] * pi_sqrt[inew] / pi_sqrt[jnew];
		    jnew++;
		  }
	      U[inew * nnew + inew] = Q[i * n + i];
	      inew++;
	    }
	}
      eigenRealSym (U, nnew, Root, V);

      for (i = n - 1, inew = nnew - 1; i >= 0; i--)	/* construct Root */
	Root[i] = (pi[i] > small ? Root[inew--] : 0);
      for (i = n - 1, inew = nnew - 1; i >= 0; i--)
	{			/* construct V */
	  if (pi[i] > small)
	    {
	      for (j = n - 1, jnew = nnew - 1; j >= 0; j--)
		if (pi[j] > small)
		  {
		    V[i * n + j] = U[jnew * nnew + inew] * pi_sqrt[jnew];
		    jnew--;
		  }
		else
		  V[i * n + j] = (i == j);
	      inew--;
	    }
	  else
	    for (j = 0; j < n; j++)
	      V[i * n + j] = (i == j);
	}
      for (i = n - 1, inew = nnew - 1; i >= 0; i--)
	{			/* construct U */
	  if (pi[i] > small)
	    {
	      for (j = n - 1, jnew = nnew - 1; j >= 0; j--)
		if (pi[j] > small)
		  {
		    U[i * n + j] = U[inew * nnew + jnew] / pi_sqrt[inew];
		    jnew--;
		  }
		else
		  U[i * n + j] = (i == j);
	      inew--;
	    }
	  else
	    for (j = 0; j < n; j++)
	      U[i * n + j] = (i == j);
	}
    }

  /*   This routine works on P(t) as well as Q. */
  /*
     if(fabs(Root[0])>1e-10 && noisy) printf("Root[0] = %.5e\n",Root[0]);
     Root[0]=0;
   */
}
