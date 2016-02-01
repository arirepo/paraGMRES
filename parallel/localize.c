#include "localize.h"
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

int localize(int my_rank, int np, int nrows, int *ia, int *ja, MPI_Comm *commun, int ***receive, int **cntowner, int ***to_be_sent, int **cntsent)
{

  // locals
  int i, j, indx, owner;
  int nlocal = (int)(nrows/np);
  int phn = nlocal;
  int *tag = (int *)calloc( nrows , sizeof(int) );
  (*receive) = (int **)malloc(np * sizeof(int *));
  (*to_be_sent) = (int **)malloc(np * sizeof(int *));
  int **expect = (int **)malloc(np * sizeof(int *));
  (*cntowner) = (int *)calloc(np , sizeof(int));
  (*cntsent) = (int *)calloc(np , sizeof(int));

  int *tmp_receive = NULL, *tmp_expect = NULL;
  int msg_tag1 = 40, msg_tag2 = 60;
  MPI_Status status1[np];
  MPI_Request request1[np];
  MPI_Status status2[np];
  MPI_Request request2[np];


  //initializing to NULL (safe)
  for ( i = 0; i < np; i++) 
    {
      (*receive)[i] = NULL;
      expect[i] = NULL;
      (*to_be_sent)[i] = NULL;
    }

  // main loop for hacking ja[]
  for ( i = 0; i < nlocal; i++) // i-local always starts from zero
    for ( indx = ia[i]; indx < ia[i+1]; indx++) //ia[] always starts from zero for this process
      {
	j=ja[indx]; //find the global column (zero-based)
	owner = (int)(j/nlocal); //and the owner rank
	if (owner == my_rank)
	  ja[indx] -= (my_rank*nlocal); // NOTE: my_rank*nlocal=vect_start_indx for this process 
	else
	  {
	    if( !tag[j] ) //is it already localized??? if no do it!
	      {
		//updating receive list
		tmp_receive = (int *) realloc ((*receive)[owner],  ((*cntowner)[owner]+1)*sizeof(int));
		if(tmp_receive == NULL)
		  {
		    printf("\ncan't realloc receive list!\n");
		    fflush(stdout);
		    exit(0);
		  }
		else
		    (*receive)[owner] = tmp_receive;

		(*receive)[owner][(*cntowner)[owner]] = phn;
		//updating expect list
		tmp_expect = (int *) realloc (expect[owner],  ((*cntowner)[owner]+1) * sizeof(int));
		if(tmp_expect == NULL)
		  {
		    printf("\ncan't realloc expect list!\n");
		    fflush(stdout);
		    exit(0);
		  }
		else
		    expect[owner] = tmp_expect;

		expect[owner][(*cntowner)[owner]] = j;
		// updating countors
		tag[j] = phn;
		((*cntowner)[owner])++;
		phn++;
	      }
	    ja[indx] = tag[j];
	  }
      }

  // sending expect list to other processes
  int rq_inc1 = 0, rq_inc2 = 0; //increaments for requests
  for( owner = 0; owner < np; owner++) //loop over process which we want to send expect list
    if ( owner == my_rank) //dont send to myself
  	continue;
    else
      {
  	/*              <<<  This is just template for MPI_Isend >>> */
  	/* int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, */
  	/*               MPI_Comm comm, MPI_Request *request) */

  	//sending the numbers of expected values to that process
  	MPI_Isend(&(*cntowner)[owner], 1, MPI_INT, owner, msg_tag1, (*commun), request1+rq_inc1++);

  	//sending expected values to that process
  	MPI_Isend(expect[owner], (*cntowner)[owner], MPI_INT, owner, msg_tag2, (*commun), request2+ rq_inc2++);
      }

  // receiving to_be_sent list from other processes
  int st_inc1 = 0, st_inc2 = 0; //increaments for status  
  for( owner = 0; owner < np; owner++) //loop over process which we want to receive from them.
    if ( owner == my_rank)
      continue;
    else
      {
  	/*                <<<  This is just template for MPI_Recv >>>
             int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, */
  	/*              MPI_Comm comm, MPI_Status *status) */

  	//receiving the numbers of expected values for that process
  	MPI_Recv(&(*cntsent)[owner], 1, MPI_INT, owner, msg_tag1, (*commun), status1+st_inc1++ );
  	//once get the size of data, reallocating enough space
  	tmp_expect = (int *)realloc ((*to_be_sent)[owner], (*cntsent)[owner]*sizeof(int) );
  	if(tmp_expect == NULL)
  	  {
  	    printf("\ncan't realloc to_be_sent list!\n");
  	    exit(0);
  	  }
  	else
  	  {
  	    (*to_be_sent)[owner] = tmp_expect;
  	  }
  	// reciving data and putting them in place
  	MPI_Recv((*to_be_sent)[owner],(*cntsent)[owner], MPI_INT, owner, msg_tag2, (*commun), status2+st_inc2++ );
  	for ( i = 0; i < (*cntsent)[owner]; i++) //localizing data
  	  (*to_be_sent)[owner][i] -= (my_rank*nlocal);
      }
  
  // wait until all send and receives ar complete
  // THIS WAS OK -> MPI_Barrier( *commun );
  /* int MPI_Waitall(int count, MPI_Request array_of_requests[],  */
  /*              MPI_Status array_of_statuses[]) */
  MPI_Waitall(np-1, request1, status1);
  MPI_Waitall(np-1, request2, status2);


  // done!
  // clean-up
  for( i = 0; i < np; i++)
       if( i != my_rank)
	    free(expect[i]);
  free(expect);

  free(tag);


  //completed successfully!
  return 0;

}


//prints localized contents
int echo_localized_forms(int my_rank, int np, int nrows, int *ia, int *ja, int **receive, int *cntowner, int **to_be_sent, int *cntsent)
{

  //locals
  int i = 0, j = 0;
  int i_start, i_end;
  int nlocal = (int)(nrows / np);

  printf("my_rank = %d, np= %d, nrows = %d\n", my_rank, np, nrows);

  //showing contents of ia[] and ja[]
  for( i = 0; i < nlocal; i++)
    {
      i_start = ia[i];
      i_end = ia[i+1]-1;
      printf("ia[i=%d] = %d, ia[i+1]-1 = %d\n", i, i_start, i_end);
      for( j = i_start; j <= i_end; j++)
  	printf("ja[j=%d] = %d, ", j, ja[j]);

      printf("\n");
    }

  //showing contents of receive list
  for( i = 0; i < np; i++)
    if(my_rank == i)
      continue;
    else
      {
  	printf("from process %d, I need to get %d values which come to sit in local positions: ",i, cntowner[i]);
  	for( j = 0; j < cntowner[i]; j++)
  	  printf("%d, ", receive[i][j]);

  	printf("\n");
      }

  //showing contents of to be sent list
  for( i = 0; i < np; i++)
    if(my_rank == i)
      continue;
    else
      {
  	printf("to process %d, I need to send %d number of my local values: ",i, cntsent[i]);
  	for( j = 0; j < cntsent[i]; j++)
  	  printf("%d, ", to_be_sent[i][j]);

  	printf("\n");
      }



  /* completed successfully! */
  return 0;

}

// for gathering data in one location in the local process and synchronization among processes
// Important : assumes qk is already allocated with appropriate size.
int gather_sync(int my_rank, int np, int nrows, int *ia, int *ja, int **receive, int *cntowner, int **to_be_sent, int *cntsent, double *qk, MPI_Comm *commun)
{

  int i = 0;
  int owner = 0; //is exactly the rank starting from 0.
  int msg_tag = 40;
  MPI_Status status[np];
  MPI_Request request[np];
  double **buff_snd = (double **)malloc(np * sizeof(double *));
  double **buff_rcv = (double **)malloc(np * sizeof(double *));
 
  // sending doubles to other processes
  int rq_inc = 0; //increament for requests  
  for( owner = 0; owner < np; owner++) //loop over process which we want to send
    if ( owner == my_rank) //dont send to myself
      continue;
    else
      {
  	/*              <<<  This is just template for MPI_Isend >>> */
  	/* int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, */
  	/*               MPI_Comm comm, MPI_Request *request) */

	// first pack doubles into the send buffer
	buff_snd[owner] = (double *)malloc(cntsent[owner] * sizeof(double));
	for( i = 0; i < cntsent[owner]; i++)
	  buff_snd[owner][i] = qk[to_be_sent[owner][i]];

  	//sending expected values to that process
  	MPI_Isend(buff_snd[owner], cntsent[owner], MPI_DOUBLE, owner, msg_tag, (*commun), request+rq_inc++);
      }

  // receiving doubles from other processes and putting them in the right places
  int st_inc = 0; //increament for status    
  for( owner = 0; owner < np; owner++) //loop over process which we want to receive from them.
    if ( owner == my_rank)
      continue;
    else
      {
  	/*                <<<  This is just template for MPI_Recv >>>
			  int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, */
  	/*              MPI_Comm comm, MPI_Status *status) */

	// first allocate receive buffer
	buff_rcv[owner] = (double *)malloc(cntowner[owner] * sizeof(double));

  	// reciving data and putting them in receive buffer
  	MPI_Recv(buff_rcv[owner], cntowner[owner], MPI_DOUBLE, owner, msg_tag, (*commun), status+st_inc++ );
  	for ( i = 0; i < cntowner[owner]; i++) //putting in right place indicated by ja[]
  	  qk[receive[owner][i]] = buff_rcv[owner][i];
      }
  
  // wait until all send and receives ar complete
  // important: each process after finishing synchronization can independently goes to do the matrix product, but for this moment we just keep the bare-bone mode to get everything working!
  // THIS WAS GOOD! -> MPI_Barrier( *commun );
  /* int MPI_Waitall(int count, MPI_Request array_of_requests[],  */
  /*              MPI_Status array_of_statuses[]) */
  MPI_Waitall(np-1, request, status);
  

  //clean - ups
  for( i = 0; i < np ; i++)
    if( i != my_rank)
      {
	free(buff_snd[i]);
	free(buff_rcv[i]);
      }

  free(buff_rcv);
  free(buff_snd);

  /* completed successfully! */
  return 0;

}
