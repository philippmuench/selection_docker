//CHHS blas.h contains headers for BLAS functions
#ifndef BLAS_H
#define BLAS_H

#ifdef __cplusplus
extern "C"
{
#endif
  void dgemv_ (const char *trans,
	       const int *m,
	       const int *n,
	       const double *alpha,
	       double *a,
	       const int *lda,
	       const double *x,
	       const int *incx,
	       const double *beta, double *y, const int *incy);

  void dsyrk_ (const char *uplo,
	       const char *trans,
	       const int *n,
	       const int *k,
	       const double *alpha,
	       double *a,
	       const int *lda, const double *beta, double *c, const int *ldc);

  void dgemm_ (const char *transa,
	       const char *transb,
	       const int *m,
	       const int *n,
	       const int *k,
	       const double *alpha,
	       double *a,
	       const int *lda,
	       double *b,
	       const int *ldb, const double *beta, double *c, const int *ldc);

  double ddot_ (const int *n,
		const double *dx,
		const int *incx, const double *dy, const int *incy);

  double dnrm2_ (const int *n, const double *dx, const int *incx);

  void dscal_ (const int *n,
	       const double *alpha, double *dx, const int *incx);
  void daxpy_ (const int *n,
	       const double *alpha,
	       const double *x, const int *incx, double *y, const int *incy);

  void dsymv_ (const char *uplo, const int *n, const double *alpha, double *a,
	       const int *lda, const double *x, const int *incx,
	       const double *beta, double *y, const int *incy);

#ifdef __cplusplus
}
#endif

const int I1 = 1;		//CHHS these constants may be useful to call BLAS routines (e.g., INCX, INCY, ALPHA, ...)
const double D1 = 1.;
const double D0 = 0.;

#endif
