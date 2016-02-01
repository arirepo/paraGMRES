# Makefile to build parallel arnoldi project 
# --- macros
CC=mpiCC
CFLAGS=  -Wall -O3
OBJECTS= mpi_gmres.o matrix_reader.o localize.o sparse_blas.o util.o
LIBS = -lm

# --- targets
all: gmres
gmres: $(OBJECTS) 
	$(CC)  -o gmres  $(OBJECTS) $(LIBS)
mpi_gmres.o : mpi_gmres.c
	$(CC) $(CFLAGS) -c mpi_gmres.c
matrix_reader.o : matrix_reader.c
	$(CC) $(CFLAGS) -c matrix_reader.c
localize.o : localize.c
	$(CC) $(CFLAGS) -c localize.c
sparse_blas.o : sparse_blas.c
	$(CC) $(CFLAGS) -c sparse_blas.c
util.o : util.c
	$(CC) $(CFLAGS) -c util.c

# --- remove object and executable files
clean:
	rm -f gmres $(OBJECTS)
