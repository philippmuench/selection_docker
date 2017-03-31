//CHHS lapack.h contains LAPACK function signatures
#ifndef LAPACK_H
#define LAPACK_H

extern "C" void dsyevr_ (const char *jobz,
			 const char *range,
			 const char *uplo,
			 const int *n,
			 double *a,
			 const int *lda,
			 const double *vl,
			 const double *vu,
			 const int *il,
			 const int *iu,
			 const double *abstol,
			 int *m,
			 double *w,
			 double *z,
			 const int *ldz,
			 int *isuppz,
			 double *work,
			 const int *lwork,
			 int *iwork, const int *liwork, int *info);

extern "C" double dlange_ (const char *norm, const int *m, const int *n,
			   double *a, const int *lda, double *work);

#endif
