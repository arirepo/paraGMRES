#ifndef _SPARSE_BLAS_H_
#define _SPARSE_BLAS_H_
#include <mpi.h>

//performs sparse matrix-vector multiplication
int sparse_matmult(double *A, double *q, double *u, int *ia, int *ja,int nrowsI, MPI_Comm *commun);

int dotproduct(double *a, double *b, double *c, int nrowsI, MPI_Comm *commun);

#endif
