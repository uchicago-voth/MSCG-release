//standard libs
#include "stdio.h"
#include "stdlib.h"
#include "malloc/malloc.h"
#include "math.h"
#include "string.h"
#include "time.h"

//Gromacs trajectory libs
#include "xdrfile_xtc.h"
#include "xdrfile_trr.h"

//math libs
//#include "mkl.h" //Intel mkl
#include "gsl/gsl_bspline.h" //GNU math lib for b-splines

//Lapack functions
extern void dsyrk_(const char *uplo, const char *trans, const int *n, const int *k,
                    const double *alpha, const double *a, const int *lda, const double *beta,
                    double *c, const int *ldc);

extern void dgemv_(const char *trans, const int *m, const int *n, const double *alpha,
                   const double *a, const int *lda, const double *x, const int *incx,
                   const double *beta, double *y, const int *incy);

extern void dgesvd_(char* jobu, char* jobvt, int* m, int* n, double* a, int* lda,
                    double* s, double* u, int* ldu, double* vt, int* ldvt, double* work,
                    int* lwork, int* info );

extern void dgeqrf_(int * m, int* n, double* a, int* lda, double* tau, double* work,
                    int* lwork, int* info );

extern void dgelsd_(int* m, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb,
                    double* s, double* rcond, int* rank, double* work, int* lwork, int* iwork, int* info );

extern void dgelss_(int* m, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, double* s, 
                    double* rcond, int* rank, double* work, int* lwork, int* info );


//the solver; google "least square QR" for details
extern void pda_lsqr_(const int *m, const int *n, const double *damp,
                      const int *leniw, const int *lenjw, const int *lenrw,
                      const int *iw, const int *jw, const double *rw,
                      double *u, double *v, double *w, double *x, double *se,
                      const double *atol, const double *btol, const double *conlim,
                      const int *itnlim, int *istop, int *itn, double *anorm,
                      double *acond, double *rnorm, double *arnorm, double *xnorm);

//User defined functions
#include "types.h" //structures and typedefs
#include "control.h" //read input file control.in
#include "top.h" //read input file top.in
#include "trj.h" //read trajectory and build cell linked-list
#include "matrix.h" //matrix computation and equation solving
#include "compute.h" //FM calculation
#include "output.h" //write output files

