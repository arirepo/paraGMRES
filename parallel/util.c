#include <stdio.h>
#include <stdlib.h>
#include "util.h"

void print_prog_syntax(void)
{

  printf("\n not enough or wrong input arguments. \n syntax: mpirun -np p ./arnoldi [input matrix file] \n"); 
}

void print_matrix_double(const char *name, double *inmat, int n1, int n2)
{
     int i,j;
     printf("\n contents of %s : \n", name);
     for( i = 0; i < n1; i++)
     {
	  for( j = 0; j < n2; j++)
	       printf("%14.14f, ", inmat[i*n2+j]);
	  printf("\n"); 

     }
     //done!
}

void print_matrix_int(const char *name, int *inmat, int n1, int n2)
{
     int i,j;
     printf("\n contents of %s : \n", name);
     for( i = 0; i < n1; i++)
     {
	  for( j = 0; j < n2; j++)
	       printf("%d, ", inmat[i*n2+j]);
	  printf("\n"); 

     }
     //done!
}

void print_array_double(const char *name, double *inmat, int n1)
{
     int i;
     printf("\n contents of %s : \n", name);
     for( i = 0; i < n1; i++)
       {
	  printf("%e \n", inmat[i]);
	  fflush(stdout);
       }
     //done!
}

void print_array_int(const char *name, int *inmat, int n1)
{
     int i;
     printf("\n contents of %s : \n", name);
     for( i = 0; i < n1; i++)
	  printf("%d \n", inmat[i]);
     //done!
}


int sum_int(int *array, int len)
{

  int i = 0, tmp = 0;
  for( i = 0; i < len; i++)
    tmp += array[i];

  return tmp;

}

