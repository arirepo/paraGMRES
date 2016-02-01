#include "matrix_reader.h"
#include <stdio.h>
#include <stdlib.h>

int matrix_reader(char *filename, int np, int my_rank, int *nnz, int *my_nnz, int *nrows, int **ia, int **ja, double **A)
{
     //locals
     int i, i_start, i_end;
     FILE *fid = fopen(filename, "r");
     int irow, icol;
     double val = 0.;
     int DI = 0;

     //resetting (to be safe!)
     *my_nnz = 0;
     *nrows = 0;
     *nnz = 0;
     
     //reading the header
     fscanf(fid, "%d %d %d\n", nrows, nrows, nnz);

     // check the divisibility of the number of rows (nrows) to the number of processes (np)
     if ( ((*nrows) % np) != 0 ) //not divisible!
     {
	  if( my_rank == 0)
	  {
	       printf("\n number of rows of the input matrix is not divisible by the number of processes. please choose another 'p' and try again.\n");
	       exit(0);
	  }
	  else
	  {
	       exit(0);
	  }
     }

     DI = (int)(*nrows / np);
     int *col_at_row = (int *)calloc( DI , sizeof(int) );
          
     //setting the range which I should read the data
     i_start = my_rank * DI;
     i_end = (my_rank+1) * DI;
     
     //scanning through the file and finding the number of items     
     for( i = 0; i < (*nnz); i++) //HINT: it can be done with EOF checking method
     {
	  fscanf(fid, "%d %d %lg\n", &irow, &icol, &val);
	  irow--; icol--; //converting to zero-based
	  if ( (i_start <= irow) && (irow < i_end) )
	       col_at_row[(irow-i_start)]++; //that row has this additional column!
     }
     //allocate (*ia)
     (*ia) = (int *)malloc( (DI+1) * sizeof(int) );

     //find my total number of non-zero items (my_nnz) and also update ia[]
     for(i= 0; i < DI; i++)
     {
	  (*ia)[i] = (*my_nnz); 	  
	  (*my_nnz) += col_at_row[i];
     }
     (*ia)[i] = (*my_nnz);
     
     //allocate (*A) one big chunck of doubles for local CRS array for storing data
     (*A) = (double *)malloc( (*my_nnz) * sizeof(double) );

     //allocate (*ja)
     (*ja) = (int *)malloc( (*my_nnz) * sizeof(int) );

     //rescanning through the file to fill data matrix A[] and ia[] and ja[]
     rewind(fid);
     //reading the header
     fscanf(fid, "%d %d %d\n", nrows, nrows, nnz);     
     //reset col_at_row[]
     for(i= 0; i < DI; i++)
	  col_at_row[i] = 0;
     
     for( i = 0; i < (*nnz); i++)
     {
	  fscanf(fid, "%d %d %lg\n", &irow, &icol, &val);
	  irow--; icol--; //converting to zero-based
	  if ( (i_start <= irow) && (irow < i_end) )
	  {
	       irow -= i_start;
	       (*A)[(*ia)[irow] + col_at_row[irow]] = val;
	       (*ja)[(*ia)[irow] + col_at_row[irow]] = icol;	       
	       col_at_row[irow]++; //GENERAL FORM :) !!!
	  }
     }
     
     //closing the input file
     fclose(fid);

     //clean-up
     free(col_at_row);
     
     //completed successfully!
     return 0;

}

//prints local sparse matrix
int echo_matrix(int my_nnz, int nrowsI, int *ia, int *ja, double *A)
{
     int i = 0, j = 0;
     int i_start, i_end;
     
     printf("my total number of entries is: %d \n", my_nnz);
     printf("my total number of rows is: %d \n", nrowsI);
     for( i = 0; i < nrowsI; i++)
     {
	  i_start = ia[i];
	  i_end = ia[i+1]-1;
	  
	  printf("ia[%d] = %d \n", i, i_start);
	  printf("ia[%d+1]-1 = %d \n", i, i_end);
	  for( j = i_start; j <= i_end; j++)
	  {
	       printf("ja[%d] = %d , val = %19.19f \n", j, ja[j], A[j]);
     
	  
	  }
     }
     
     //completed successfully!
     return 0;
     
}

//finds the maximum of the given integer array with size "n"
int max_int_array(int *input, int n)
{
  //locals
  int i = 0;
  int max = input[0];

  for( i = 0; i < n ; i++)
    if( input[i] > max)
      max = input[i];

  //returning the maximum value
  return max;

}

//reads the rhs vector from input file.
int read_rhs(char *filename, int np, int my_rank, int nrows, int DI, double *b)
{
     //locals
     int i;
     FILE *fid = fopen(filename, "r");
     int i_start, i_end, irow;
     double val = 0.;
     int nrows_file = 0, ncols_file = 0;
     
     //reading the header
     fscanf(fid, "%d %d\n", &nrows_file, &ncols_file);

     //checking  ...
     if ( nrows_file != nrows )
     {
	  printf("FATAL error: number of rows in the rhs file doesn't matches with the number of rows of the input matrix. Try different rhs file.\n");
	  exit(0);
     }
     else if ( ncols_file != 1)
     {
	  printf("FATAL error: number of columns in the rhs file must be one.\n");
	  exit(0);
     }
     else
     {
	  //OK, we are all set!
     }
          
     //setting the range which I should read the data
     i_start = my_rank * DI;
     i_end = (my_rank+1) * DI;
     
     //scanning through the file
     irow = 0;
     for( i = 0; i < nrows; i++) 
     {
	  fscanf(fid, "%lg\n", &val);
	  
	  if ( (i_start <= irow) && (irow < i_end) )
	       b[(irow-i_start)] = val; //put it in local [b]
	  
	  irow++;
     }
     
     //closing the input file
     fclose(fid);

     
     //completed successfully!
     return 0;

}
