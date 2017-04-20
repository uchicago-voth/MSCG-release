//
//  external_matrix_routines.h
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#ifndef _external_matrix_routines_h
#define _external_matrix_routines_h

# ifdef __cplusplus
extern "C" {
# endif

# if _mkl_flag == 1
	#include "mkl.h"
# else
	#include <gsl/gsl_blas.h>
#  undef _mkl_flag
#  define _mkl_flag 0
# endif

extern double cblas_ddot(const int n, const double* dx, const int incx, const double* dy, const int incy);

	
# if _mkl_flag == 0
// Exclude these function definitions when compiling with MKL

extern void dgesvd_(char* jobu, char* jobvt, int* m, int* n, double* a, int* lda,
                    double* s, double* u, int* ldu, double* vt, int* ldvt, double* lapack_temp_workspace,
                    int* lapack_setup_flag, int* info);

extern void dgetrf_(const int* m, const int* n, double* a, const int* lda, int* ipiv, const int* info);

extern void dgeqrf_(int* m, int* n, double* a, int* lda, double* lapack_tau, double* lapack_temp_workspace,
                    int* lapack_setup_flag, int* info);

extern void dgelsd_(int* m, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb,
                    double* s, double* rcond, int* rank, double* lapack_temp_workspace, int* lapack_setup_flag, int* iwork, int* info);

extern void dgelss_(int* m, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, double* s,
                    double* rcond, int* rank, double* lapack_temp_workspace, int* lapack_setup_flag, int* info);

extern void dgetri_(const int* n, double* a, const int* lda, int* ipiv, double* work, const int* lwork, int *info);

# endif
					
#ifdef __cplusplus
}
#endif

#endif
