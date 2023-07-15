#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>
#include <time.h>
#define DIM 20000
#define ITER 1000
#define MULTI 1
int main( int argc, char **argv ) {
    MPI_Init(&argc, &argv);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    char** matrix_multi;
    char* matrix_plain;
    matrix_plain = calloc(DIM*DIM, sizeof(char));
    matrix_multi = malloc(DIM*sizeof(char*));
    double start, end;
    start = MPI_Wtime();

    #if (MULTI == 1) 
        for (int i = 0; i < DIM; i++) {
            matrix_multi[i] = &matrix_plain[i*DIM];
        }
    #endif
    
    for (int k = 0; k < ITER; k++)
        for (int i = 0; i < DIM; i++) 
            for (int j = 0; j < DIM; j++) {
                #if (MULTI == 1) 
                    matrix_multi[i][j] = 1;
                #else 
                    matrix_plain[i*DIM+j] = 1;
                #endif
            }
    end = MPI_Wtime();
    printf("Time: %f\n", end-start);
    MPI_Finalize();
    return 0;
}