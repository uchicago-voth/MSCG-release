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
# else
#  undef _mkl_flag
#  define _mkl_flag 0
# endif

extern void dsyrk_(const char* uplo, const char* trans, const int* n, const int* k,
                   const double* alpha, const double* a, const int* lda, const double* beta,
                   double* c, const int* ldc);

extern void dgemv_(const char* trans, const int* m, const int* n, const double* alpha,
                   const double* a, const int* lda, const double* x, const int* incx,
                   const double* beta, double* y, const int* incy);
	
# if _mkl_flag == 0
// Exclude these function definitions when compiling with MKL
extern void daxpy( const int* lda, const double* a, const double* beta, const int* ldb, const double* b, const int* ldc);

extern void dgesvd_(char* jobu, char* jobvt, int* m, int* n, double* a, int* lda,
                    double* s, double* u, int* ldu, double* vt, int* ldvt, double* lapack_temp_workspace,
                    int* lapack_setup_flag, int* info);

extern void dgeqrf_(int* m, int* n, double* a, int* lda, double* lapack_tau, double* lapack_temp_workspace,
                    int* lapack_setup_flag, int* info);

extern void dgelsd_(int* m, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb,
                    double* s, double* rcond, int* rank, double* lapack_temp_workspace, int* lapack_setup_flag, int* iwork, int* info);

extern void dgelss_(int* m, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, double* s,
                    double* rcond, int* rank, double* lapack_temp_workspace, int* lapack_setup_flag, int* info);
# endif
					
#ifdef __cplusplus
}
#endif

#endif
