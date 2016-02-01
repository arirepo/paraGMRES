#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "matrix_reader.h"
#include "localize.h"
#include "sparse_blas.h"
#include "util.h"

int main(int argc, char *argv[])
{
  /* local vars */
  int my_rank = 0, p = 0, root = 0;

  MPI_Comm *commun = (MPI_Comm *)malloc(1*sizeof(MPI_Comm));
  commun[0] = MPI_COMM_WORLD;

  MPI_Init(&argc,&argv); //initializing MPI environment
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank); //getting the current rank
  MPI_Comm_size(MPI_COMM_WORLD,&p); //getting number of processes 

  int i = 0, j = 0, k, ss,  k_end = atoi(argv[3]);      
     
  //setting appropriate output environment
  char stdoutname[128];
  FILE *newstdout;
  sprintf(stdoutname,"out.%d",my_rank);
  newstdout = freopen(stdoutname,"w",stdout);
  if (newstdout == NULL)
    {
      printf("Failed to redirect stdout to %s.\n",stdoutname);
    }

  ////////////////// CODE BODY /////////////////////////  
  //reading the portion of
  int nnz, my_nnz, nrows;
  int *ia = NULL, *ja = NULL;
  double *A = NULL;
  int nrowsI = 0;
  matrix_reader(argv[1], p, my_rank, &nnz, &my_nnz, &nrows, &ia, &ja, &A);
  nrowsI = (int)(nrows/p);

  //echo_matrix(my_nnz, nrowsI, ia, ja, A);
	  
  // reading the rhs vector from the input file
  double *b = (double *)calloc(  nrowsI , sizeof(double) );
  for ( i = 0 ; i < nrowsI; i++ )
    b[i] = 1.;
     
  //read_rhs(argv[2], p, my_rank , nrows, nrowsI , b);     
  //print_matrix_double("b", b, nrowsI , 1);

  //localize
  int **receive = NULL;
  int *cntowner= NULL;
  int **to_be_sent= NULL;
  int *cntsent = NULL;

  localize(my_rank, p, nrows, ia, ja, commun, &receive, &cntowner, &to_be_sent, &cntsent);
     
  //prints localized contents
  //     echo_localized_forms(my_rank, p, nrows, ia, ja, receive, cntowner, to_be_sent, cntsent);

  double t_start = 0., t_end = 0.;  
  t_start = MPI_Wtime(); //put a time mark!
  
  int precond_flag = atoi(argv[4]);
  double gmres_res = 1.e-10; //gmres residual
  int n_total = max_int_array(ja, my_nnz)+1; //ntotal is equal to nlocal+n phantom points
  double norm_r0 = 0.;

  double *u = (double *)calloc( nrowsI , sizeof(double) ); //this memory would be enough for all! 
  double **v = (double **)malloc( 1 * sizeof(double *) ); //we'll add appropriate columns inside gmres loop
  double **h = (double **)malloc( 1 * sizeof(double *) ); //we'll add appropriate columns inside gmres loop
  double *x0 = (double *)calloc( n_total , sizeof(double) );      
  double *r0 = (double *)calloc( nrowsI , sizeof(double) );      
  double *g = (double *)calloc( 1 , sizeof(double) ); //we'll add additional allocations inside gmres loop
  double *c = NULL, *s = NULL, *resi = NULL;
  double delta = 0., gamma = 0.;
  double *X = (double *)calloc(nrowsI , sizeof(double) );

  //pick the diagonal entries for preconditioning (if it is preconditioned)
  double *diag = (double *) calloc( nrowsI , sizeof(double) );
  double *yk = (double *) calloc( n_total , sizeof(double) );          
  if( precond_flag ) //should it be preconditioned?
    {
      for ( i = 0 ; i < nrowsI; i++ )
	for ( j = ia[i] ; j < ia[i+1]; j++ )
	  if( ja[j] == i ) //then it's on the main diagonal
	    diag[i] = A[j];
      //print_matrix_double("diag", diag, nrowsI , 1);
    }
     
  //init x0
  for( i = 0; i < nrowsI; i++)
    x0[i] = 1.0;

  //restarting vars
  int rst = 0, mm = atoi(argv[5]), completed = 0, total_itr = 1;
     
  //do restarting
  for( rst = 0; rst < mm; rst++)
    {

      //gather/sync 
      gather_sync(my_rank, p, nrows, ia, ja, receive, cntowner, to_be_sent, cntsent, x0, commun);
      /* //matrix vector product */
      sparse_matmult(A, x0, r0, ia, ja, nrowsI, commun);
      //updating the final r0
      for( i = 0; i < nrowsI; i++)
	r0[i] = b[i] - r0[i];

      // normalization ...
      dotproduct(r0, r0, &norm_r0, nrowsI, commun); //dot product
      norm_r0 = sqrt(norm_r0); //root
      for( i = 0; i < nrowsI; i++)
	r0[i] /= norm_r0; //normalizing
          
      // initial allocation of v
      h[0] = (double *)calloc( 1 , sizeof(double) );           
      v[0] = (double *)calloc( n_total , sizeof(double) );
      for( i = 0; i < nrowsI; i++)
	v[0][i] = r0[i]; //initializing ...

      g[0] = norm_r0; //initial g
      //report the residual on the root process
      /* if ( my_rank == root ) */
      /* 	printf("%d %e\n", 1, fabs(g[0])); */
	       
      //main gmres loop
      for( k = 0; k < k_end; k++) //will be fancier soon! 
	{
	  //reallocating vars to have enough space for the next items
	  g = (double *) realloc (g,  (k+2) * sizeof(double) );
	  g[k+1] = 0.;
	  c = (double *) realloc (c,  (k+1) * sizeof(double) ); 
	  s = (double *) realloc (s,  (k+1) * sizeof(double) );
	  resi = (double *) realloc (resi,  (k+1) * sizeof(double) );
	  //reallocating v and h
	  v = (double **)realloc (v , (k+2) * sizeof(double *) );
	  v[k+1] = (double *)calloc( n_total , sizeof(double) );
	  h = (double **)realloc (h , (k+2) * sizeof(double *) );
	  h[k+1] = (double *)calloc( 1 , sizeof(double) );     
	  
	  for ( j = 0 ; j <= (k+1); j++)
	    h[j] = (double *)realloc( h[j], (k+1) *sizeof(double) );

	  //NOTE : THE FOLLOWING IF STATEMENT INCREASES PERFORMANCE WHEN IT IS
	  //NOT SUPPOSED TO BE PRECONDITIONED BECAUSE IT JUST SKIPS THE INITIALIZATION
	  // OF VECTOR [yk].
	  if( precond_flag ) //should it be preconditioned? 
	    {
	      for ( i = 0 ; i < nrowsI; i++ )
		if( diag[i] )
		  yk[i] = v[k][i] / diag[i];
		else
		  yk[i] = v[k][i];
	      //gather/sync
	      gather_sync(my_rank, p, nrows, ia, ja, receive, cntowner, to_be_sent, cntsent, yk, commun);
	      /* //matrix vector product */
	      sparse_matmult(A, yk, u, ia, ja, nrowsI, commun);
	    }
	  else
	    {
	      //gather/sync
	      gather_sync(my_rank, p, nrows, ia, ja, receive, cntowner, to_be_sent, cntsent, v[k], commun);
	      /* //matrix vector product */
	      sparse_matmult(A, v[k], u, ia, ja, nrowsI, commun);
	    }
	  
	  for ( j = 0 ; j <= k; j++)
	    {
	      dotproduct(v[j], u, &h[j][k], nrowsI, commun); //dot product
	      for( ss = 0; ss < nrowsI; ss++)
		u[ss] -= h[j][k] * v[j][ss];
	    }
	  dotproduct(u, u, &h[k+1][k], nrowsI, commun); //dot product
	  h[k+1][k] = sqrt(h[k+1][k]); //norm
	  //updating v[k+1]
	  for( ss = 0; ss < nrowsI; ss++)
	    v[k+1][ss] = u[ss] / h[k+1][k];

	  for ( j = 0 ; j < k; j++)
	    {
	      delta = h[j][k];
	      h[j][k] = c[j] * delta + s[j] * h[j+1][k];
	      h[j+1][k] = -s[j] * delta + c[j] * h[j+1][k];
	    }
	  gamma = sqrt(h[k][k] * h[k][k] + h[k+1][k]*h[k+1][k]);
	  c[k] = h[k][k] / gamma;
	  s[k] = h[k+1][k] / gamma;
	  h[k][k] = gamma;
	  h[k+1][k] = 0.;
	  delta = g[k];
	  g[k] = c[k] * delta + s[k] * g[k+1];
	  g[k+1] = -s[k] * delta + c[k] * g[k+1];
	  resi[k] = fabs(g[k+1]);
	  //report the residual on the root process
	  /* if ( my_rank == root ) */
	  /*   printf("%d %e\n", total_itr++, resi[k]); */
	  if ( resi[k] <= gmres_res )
	  {
	       completed = 1; //set the completed flag
	       break; //terminate gmres loop ASAP!
	  }
	} //end of the main gmres loop
      if(completed)
	   break; //terminate restart loop ASAP!
      
      k--; //reducing k to emulate the inside the loop effect
	  //compute alpha
	  //allocate alpha
      double *alpha = (double *)calloc(k+1 , sizeof(double) );
      //solve backward
      for ( j = k ; j >= 0; j--)
	{
	  alpha[j] = g[j]/h[j][j];
	  for ( ss = (j+1) ; ss <= k; ss++)
	    alpha[j] -= (h[j][ss]/h[j][j] * alpha[ss]);
	}
      //compute zk
      double *zk = (double *)calloc(nrowsI , sizeof(double) );
      for ( j = 0 ; j <= k; j++)
	for ( ss = 0 ; ss < nrowsI; ss++)
	  zk[ss] += alpha[j] * v[j][ss];

      if( precond_flag ) //should it be preconditioned?
	for ( i = 0 ; i < nrowsI; i++ )
	  if ( diag[i] )
	    zk[i] /= diag[i];
	       
      //compute solution
      for ( ss = 0 ; ss < nrowsI; ss++)
	X[ss] = x0[ss] + zk[ss];
      //report the solution

      //preparinf for restart operation ...
      //put x into x0
      for( i = 0; i < nrowsI; i++)
	x0[i] = X[i]; //put last [X] into x0 
	  
      //clean ups
      free(alpha);
      free(zk);

      for ( j = 0 ; j <= (k+1); j++)
	{
	  free(v[j]);
	  free(h[j]);	       
	}
	  
    } //end of restarting loop

  t_end = MPI_Wtime(); //put a time mark!

  printf("Entire program run time including reading the input files and clean-up is : %f seconds.\n", (t_end - t_start) );

  //int_matrix_double("x", X, nrowsI, 1);
  

  ///////////////END OF CODE BODY //////////////////////
     
  //closing outputs
  fclose(newstdout);     
  MPI_Finalize(); //finishing MPI session.


  /* completed successfully! */
  return 0;

}

