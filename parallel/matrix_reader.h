#ifndef _MATRIX_READER_H
#define _MATRIX_READER_H_

int matrix_reader(char *filename, int np, int my_rank, int *nnz, int *my_nnz, int *nrows, int **ia, int **ja, double **A);

//prints local sparse matrix
int echo_matrix(int my_nnz, int nrowsI, int *ia, int *ja, double *A);

int max_int_array(int *input, int n);

//reads the rhs vector from input file.
int read_rhs(char *filename, int np, int my_rank, int nrows, int DI, double *b);


#endif
     
