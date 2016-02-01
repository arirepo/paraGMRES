#ifndef _LOCALIZE_H_
#define _LOCALIZE_H_
#include <mpi.h>

int localize(int my_rank, int np, int nrows, int *ia, int *ja, MPI_Comm *commun, int ***receive, int **cntowner, int ***to_be_sent, int **cntsent);

int echo_localized_forms(int my_rank, int np, int nrows, int *ia, int *ja, int **receive, int *cntowner, int **to_be_sent, int *cntsent);

int gather_sync(int my_rank, int np, int nrows, int *ia, int *ja, int **receive, int *cntowner, int **to_be_sent, int *cntsent, double *qk, MPI_Comm *commun);

#endif
