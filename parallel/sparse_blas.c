#include <stdio.h>
#include "sparse_blas.h"
#include <mpi.h>

//performs sparse matrix-vector multiplication
int sparse_matmult(double *A, double *q, double *u, int *ia, int *ja,int nrowsI, MPI_Comm *commun)
{

  //locals
  int r, j;
  // reset u
  for( r = 0 ; r < nrowsI; r++)
    u[r] = 0.;

  for( r = 0 ; r < nrowsI; r++)
      for ( j = ia[r]; j < ia[r+1]; j++)
	{
	  //printf("(u[r=%d]=%f) += ((A[j=%d]=%f) * (q[ja[j]=%d]=%f))\n", r, u[r], j, A[j], ja[j], q[ja[j]]); 
	  u[r] += (A[j] * q[ja[j]]);
	}

  //explicit barrier
  //MPI_Barrier( *commun );

  //completed successfully!
  return 0;

} 

//performs dot product of the given vectors "a" and "b" (assuming chunks stored COMPLETELY per process) and stored the single scalar value obtaing by collective calls among process in "c".
int dotproduct(double *a, double *b, double *c, int nrowsI, MPI_Comm *commun)
{
  //locals
  int i = 0;
  double local_dot = 0.;
  //*c = 0.; //hard reset
  
  //perform local dot product
  for( i = 0; i < nrowsI; i++)
    local_dot += (a[i] * b[i]);

  //reduceall to single values
  /* int MPI_Allreduce( */
  /*   void *sendbuf, */
  /*   void *recvbuf, */
  /*   int count, */
  /*   MPI_Datatype datatype, */
  /*   MPI_Op op, */
  /*   MPI_Comm comm */
  /* ); */

  MPI_Allreduce(
		&local_dot,
		c,
		1,
		MPI_DOUBLE,
		MPI_SUM,
		(*commun));

  //completed successfully!
  return 0;

}

