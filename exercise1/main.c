#if defined(__STDC__)
    #if (__STDC_VERSION__ >= 199901L)
        #define _XOPEN_SOURCE 700
    #endif
#endif
#if !defined(_OPENMP)
    #error "OpenMP support needed for this code"
#endif

#define PROJECT_RULES 0
#define EXTENSION ".pgm"
#define IMAGETEST 0 //This macro enable a mode that prints a glider matrix

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <stdbool.h>
#include <string.h>
#include "read_write_pgm_image.h"



//-----------TEST RESULTS------------
//NOT OPTIMIZED----------------------
//srun project -i -k 15 -e 1 -f fava -n 70 -s 0     WORKS
//srun project -i -k 16 -e 1 -f fava -n 70 -s 0     WORKS
//srun project -i -k 17 -e 1 -f fava -n 70 -s 0     WORKS
//srun project -i -k 15 -e 1 -f fava -n 70 -s 1     WORKS  
//srun project -i -k 16 -e 1 -f fava -n 70 -s 1     WORKS 
//srun project -i -k 16 -e 0 -f fava -n 70 -s 1     WORKS
//srun project -r -k 16 -e 1 -f fava -n 70 -s 1     WORKS
//srun project -r -k 16 -e 0 -f fava -n 70 -s 1     WORKS
char* snapshotName(int n);
void orderedUpdate(unsigned char* field, long int k);
unsigned char* staticUpdate(unsigned char* block, int block_rows, long int k, unsigned char* blockB);
void updateMargins(unsigned char* block,long int k, int block_rows, int world_rank, int world_size, MPI_Datatype row);

int main( int argc, char **argv )
{
    //Initialize MPI environment
    MPI_Init(&argc, &argv);
    //get number of processes
    double start, end;
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    MPI_Barrier(MPI_COMM_WORLD);
    //printf("Version 4.1\n");
    start = MPI_Wtime();


    int maxval = MAXVAL;
    bool initPlayground = 1; //Flag to state if is a new game or a loaded one. Default is new
    long int k = 0; //Size of playground, must be specified
    int e = 1; //Update mode, default is 1 (static)
    char* file = NULL; //String containing file name, must be specified
    int n = 0; //Number of steps of the computation, must be specified else the program will do nothing
    int s = 0; //Number of steps to take a snapshot of the system. Default is 0 (no snapshot)

    for (int i = 0; i < argc; i++)
        if (argv[i][0] == '-') {
            if (argv[i][1] == 'i')
                initPlayground = 1;
            else if (argv[i][1] == 'r')
                initPlayground = 0;
            else if (argv[i][1] == 'k')
                k = atoi(argv[i+1]);
            else if (argv[i][1] == 'e')
                e = atoi(argv[i+1]);
            else if (argv[i][1] == 'f')
                file = argv[i+1];
            else if (argv[i][1] == 'n')
                n = atoi(argv[i+1]);
            else if (argv[i][1] == 's')
                s = atoi(argv[i+1]);
            else {
                printf("ERROR: invalid argument");
                return -1;
            }
        }
    unsigned char* block = NULL;

    if (e) {
        int* blocksizes = NULL; //An array to store the size of every block of data passed to the process with scatter
        int* displacements = NULL; //An array to store the displacements of the blocks 
        if((blocksizes = malloc(world_size*sizeof(int))) == NULL) return -1;
        if((displacements = malloc(world_size*sizeof(int))) == NULL) return -1;
        int block_rows;
        MPI_Datatype row;
        if (initPlayground) {
            MPI_Type_contiguous(k, MPI_UNSIGNED_CHAR, &row);
            MPI_Type_commit(&row);
            makeBlockSizes(k, world_size, displacements, blocksizes); //These are useful to not calculate every time the displacements and the blocksizes
            block_rows = blocksizes[world_rank];
            if ((block = malloc((block_rows * k)*sizeof( unsigned char)))==NULL) return -1; //Allocate the subgrid
        } 
        else {
            #if IMAGETEST != 1
                if (!initPlayground) {
                    block = parallel_read(&maxval, &k, file, world_rank, world_size, &block_rows, displacements, blocksizes, &row);
                }
            #else
                MPI_Type_contiguous(k, MPI_UNSIGNED_CHAR, &row);
                MPI_Type_commit(&row);
                makeBlockSizes(k, world_size, displacements, blocksizes);
                block_rows = blocksizes[world_rank];    
                if ((block = malloc((block_rows * k)*sizeof( unsigned char)))==NULL) return -1;    
            #endif
        }

        #if IMAGETEST == 1 
        unsigned char* workBlocks = NULL;
        if (!initPlayground) {
            if (world_rank == 0) {     
                //This area of the memory will contain all the workblocks being passed to the processes with scatter
                if((workBlocks = malloc(((k/world_size + 2) * k)*world_size*sizeof(unsigned char) + (k % world_size) * k * sizeof(unsigned char)))==NULL) return -1;
                unsigned char* image = NULL;
                if((image = calloc(k * k, sizeof(unsigned char*))) == NULL) return -1;
                for (int i = 0; i < k-4; i+=7) {
                    image[(i+0) * k + 2] = maxval;
                    image[(i+1) * k + 2] = maxval;
                    image[(i+2) * k + 2] = maxval;
                    image[(i+2) * k + 1] = maxval;
                    image[(i+1) * k + 0] = maxval;
                }
                printf("2--------------\n");
                makeBlockSizes(k, world_size, displacements, blocksizes);
                printf("2.1------------\n");
                size_t workBlocksOffset;
                int imageCut;
                for (int blockCount = 0; blockCount < world_size; blockCount++) {//We fill the workblocks vector block by block
                    #pragma omp parallel
                    {
                        workBlocksOffset = displacements[blockCount]*k;
                        imageCut = displacements[blockCount] -1 -(2*blockCount);
                        #pragma omp for collapse(2) 
                            for (int i = imageCut; i < imageCut + blocksizes[blockCount]; i++) //We are attaching the external rows contiguosly to the subgrid, so we start 1 row before and finish 1 row after 
                                for (int j = 0; j < k; j++)                     
                                    workBlocks[workBlocksOffset + ((i - imageCut) * k + j)] = image[((i >= 0)? (i%k) : (k -1)) * k + j]; //Here the ternary operator defines the toroid behauviour                        
                    }
                }
                printf("2.5------------\n");
                
            }

            //DEBUG
            //for (int i= 0; i < world_size; i++) printf("%d, ", displacements[i]); printf("\n");
            //for (int i= 0; i < world_size; i++) printf("%d, ", blocksizes[i]); printf("\n");
            //if (world_rank == 0) for (int i= 0; i < k*k + 2*k*world_size; i++) printf("%d, ", workBlocks[i]); printf("\n");

            printf("\n%d\n", block_rows);
            //DEBUG
            MPI_Scatterv(workBlocks, blocksizes, displacements, row, block, block_rows, row, 0, MPI_COMM_WORLD);
            printf("2.6------------\n");
        }
        #endif



        if (initPlayground)  //Fill the playground with random values
            for (int i = 1; i < block_rows-1; i++) 
                for (int j= 0; j < k; j++) 
                    block[i*k + j] = rand() % 2 == 0? MAXVAL : 0;

        updateMargins(block, k, block_rows, world_rank, world_size, row);   //We update the margins because every process is blind about the others

        unsigned char* tmp = NULL;
        unsigned char* blockB = malloc(((block_rows) * k)*sizeof(unsigned char));


        if (blockB == NULL) return -1;

        for (int i = 0; i < n; i++) {
            tmp = block;
            block = staticUpdate(block, block_rows, k, blockB);//not the most elegant way but is efficient
            blockB = tmp; //swap
            updateMargins(block, k, block_rows, world_rank, world_size, row); //Update the external rows
            if (s > 0 && (i+1) % s == 0) 
                efficient_write2(block, MAXVAL, k, snapshotName(i), world_rank, world_size, block_rows, displacements, blocksizes, row); //Different routine for snapshot saving
        }
        
        MPI_Barrier(MPI_COMM_WORLD); //If you take out this the file writing goes wrong
        efficient_write(block, MAXVAL, k, "end.pgm", world_rank, world_size, block_rows, displacements, blocksizes, row);
    }
    else {
        if (world_rank == 0) {
            block = malloc(k*k*sizeof(unsigned char));
            if (initPlayground) 
                for (int i = 0; i < k*k; i++) 
                    block[i] = rand() % 2 == 0? MAXVAL : 0;
            else 
                block = read_pgm_image(&maxval, &k, &k, file);
                
            for (int i = 0; i < n; i++) {
                orderedUpdate(block, k);
                if (s > 0 && (i+1) % s == 0)
                    write_pgm_image(block, MAXVAL, k, snapshotName(i)); 
            }
            
            write_pgm_image(block, MAXVAL, k, "end.pgm"); 
        }
    }
    
    

    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    if (world_rank == 0) {
        printf("%f", end-start);
    }
    MPI_Finalize();
    
    
    return 0;
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//--------------------------------------------------END OF MAIN-----------------------------------------------------------------
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

void updateMargins(unsigned char* block,long int k, int block_rows, int world_rank, int world_size, const MPI_Datatype row) {
    MPI_Request req_lower_send, req_upper_send, req_lower_recv, req_upper_recv;
    MPI_Status stat_lower_send, stat_upper_send, stat_lower_recv, stat_upper_recv;
    //Non blocking receiving calls, using module and ternary operator to respect the toroid behaviour
    MPI_Irecv(&block[(block_rows -1) * k] , 1, row, (world_rank + 1) % world_size, 1, MPI_COMM_WORLD, &req_lower_recv);
    MPI_Irecv(block, 1, row, world_rank == 0? world_size -1 : world_rank - 1, 0, MPI_COMM_WORLD, &req_upper_recv);
    //Non blocking sending calls, using module and ternary operator to respect the toroid behaviour
    MPI_Isend(&block[(block_rows - 2) * k], 1, row, (world_rank + 1) % world_size, 0, MPI_COMM_WORLD, &req_lower_send);
    MPI_Isend(&block[k], 1, row, world_rank == 0? world_size -1 : world_rank - 1, 1, MPI_COMM_WORLD, &req_upper_send);
    //Waiting for completion of the request 
    MPI_Wait(&req_upper_recv, &stat_upper_recv);
    MPI_Wait(&req_lower_recv, &stat_lower_recv);
    MPI_Wait(&req_upper_send, &stat_upper_send);
    MPI_Wait(&req_lower_send, &stat_lower_send);    
}

void orderedUpdate(unsigned char* field, long int k) {
    int neighbours = 0;
    for ( int i = 0; i < k; i++ )
        for ( int j = 0; j < k; j++ ) {
            neighbours = ( field[((i > 0 ? i : k)-1) * k + (j > 0? j : k)-1] + field[((i > 0 ? i : k)-1) * k + j] + 
                        field[((i > 0 ? i : k)-1) * k + ((j + 1) % k)] +  field[i * k + ((j + 1) % k)] +
                        field[((i+1)%k) * k + ((j + 1) % k)] + field[((i + 1) % k) * k + j] + 
                        field[((i + 1) % k) * k + (j > 0? j : k)-1] + field[i * k + (j > 0? j : k)-1] )/MAXVAL;
            //Here the ternary operator and module defines the toroid behauviour
  
            //The project rules are different from the original rules, so i inserted this directive to enable to switch between official rules and project rules
            #if PROJECT_RULES == 1
            if (neighbours == 2 || neighbours ==3)
                field[i][j] = MAXVAL;
            else
                field[i][j] = 0;
            #else
            if (((neighbours == 2 || neighbours == 3) && field[i * k + j] > 0 ) || (neighbours == 3 && field[i * k + j] == 0)) //We apply the three rules in one single check
                field[i * k + j] = MAXVAL;
            else
                field[i * k + j] = 0;
            #endif 
        }
    return;
}

unsigned char* staticUpdate(unsigned char* block, int block_rows, long int k, unsigned char* blockB) {
    #pragma omp parallel 
    {
        int neighbours = 0;
        #pragma omp for collapse(2) 
            for ( int i = 1; i < block_rows - 1; i++ )
                for ( int j = 0; j < k; j++ ) {
                    neighbours = (block[(i - 1) * k + (j > 0? j : k)-1] + block[(i - 1) * k + j] + 
                                block[(i - 1) * k + (j+1)%k] +  block[i * k + (j+1)%k] +
                                block[(i + 1) * k + (j+1)%k] + block[(i + 1) * k + j] + 
                                block[(i + 1) * k + (j > 0? j : k)-1] + block[i * k + (j > 0? j : k)-1])/MAXVAL;
                    //Here the ternary operator defines the toroid behauviour

                    //The project rules are different from the original rules, so i inserted this directive to enable to switch between official rules and project rules
                    #if PROJECT_RULES == 1 
                    if (neighbours == 2 || neighbours ==3)
                        blockB[i * k + j] = MAXVAL;
                    else
                        blockB[i * k + j] = 0;
                    #else
                    //We apply the three rules in one single check
                    if (((neighbours == 2 || neighbours == 3) && block[i * k + j] > 0 ) || (neighbours == 3 && block[i * k + j] == 0)) 
                        blockB[i * k + j] = MAXVAL;
                    else 
                        blockB[i * k + j] = 0;
                    #endif
                }
    }
    return blockB;
}

char* snapshotName(int n) {
    char* s = malloc(19*sizeof(char));

    if (n < 0)
        return NULL;
    else if (n < 10) {
        sprintf(s, "snapshot_0000%d", n);
    }
    else if (n < 100) {
        sprintf(s, "snapshot_000%d", n);
    }
    else if (n < 1000) {
        sprintf(s, "snapshot_00%d", n);
    }
    else if (n < 10000) {
        sprintf(s, "snapshot_0%d", n);
    }
    else if (n < 100000) {
        sprintf(s, "snapshot_%d", n);
    }

    strcat(s, EXTENSION);
    return s;
}

