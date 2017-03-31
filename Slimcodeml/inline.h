//CHHS inline.h contains templates and function headers which may be inlined
#ifndef INLINE_H
#define INLINE_H

#include "blas.h"

inline void
zero (double x[], int n)
{
  memset (x, 0, n * sizeof (double));
}

inline void
xtoy (double x[], double y[], int n)
{
  memcpy (y, x, n * sizeof (double));
}

inline double
sum (double x[], int n)
{
  double t = 0.;

  for (int i = 0; i < n; i++)
    t += x[i];

  return t;
}

template < class T > inline T min2 (T a, T b)
{
  return (a < b) ? a : b;
}

template < class T > inline T max2 (T a, T b)
{
  return (a > b) ? a : b;
}

template < class T > inline T square (T a)
{
  return a * a;
}

inline unsigned int
spaceming2 (unsigned int n)
{
  return n * (n * 2 + 9 + 2) * sizeof (double);
}

inline void
identity (double x[], int n)
{
  memset (x, 0, n * n * sizeof (double));
  for (int i = 0; i < n * (n + 1); i += (n + 1))
    x[i] = 1.;
}

inline void
axtoy (double a, double x[], double y[], int n)
{
  for (int i = 0; i < n; ++i)
    y[i] = a * x[i];
}

inline void
abyx (double a, double x[], int n)
{
  dscal_ (&n, &a, x, &I1);	//CHHS BLAS function
}

inline void
fillxc (double x[], double c, int n)
{
  for (int i = 0; i < n; ++i)
    x[i] = c;
}

inline double
norm (double x[], int n)
{
  // Tried dnrm2_(&n, x, &I1);
  // but it changes results too much (weird!?)
  // also using ddot, same problem
  int i;
  double t = 0.;

  for (i = 0; i < n; ++i)
    t += x[i] * x[i];

  return sqrt (t);
}

inline double
distance (double x[], double y[], int n)
{
  int i;
  double t = 0.;

  for (i = 0; i < n; ++i)
    t += square (x[i] - y[i]);

  return sqrt (t);
}

inline double
innerp (double x[], double y[], int n)
{
  return ddot_ (&n, x, &I1, y, &I1);	//CHHS BLAS routine
}

/// Transpose a matrix.
///       x[n][m] --> y[m][n]
///
inline void
mattransp2 (double x[], double y[], int n, int m)
{
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j)
      y[j * n + i] = x[i * m + j];
}
#endif
