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
#define IMAGETEST 1

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
//-i -k 15 -e 1 -f fava -n 70 -s 0     WORKS
//-i -k 16 -e 1 -f fava -n 70 -s 0     WORKS
//-i -k 17 -e 1 -f fava -n 70 -s 0     WORKS
//-i -k 15 -e 1 -f fava -n 70 -s 1     WORKS  
//-i -k 16 -e 1 -f fava -n 70 -s 1     WORKS 
//-i -k 16 -e 0 -f fava -n 70 -s 1     WORKS
//-r -k 16 -e 1 -f fava -n 70 -s 1     WORKS
//-r -k 16 -e 0 -f fava -n 70 -s 1     WORKS
void collectivePrint(unsigned char* block, int world_size, int world_rank, int block_rows, int k);
void orderedUpdate(unsigned char** field, int k_i, int k);
unsigned char* staticUpdate(unsigned char* block, int block_rows, long int k, unsigned char* blockB);
char* snapshotName(int n);
char* intToString(int value, char* buffer);
char* reverse(char *buffer, int i, int j);
void swapc(char *x, char *y);
void closestFactors(int nOfProcess, int* n_i, int* n_j);
void updateMargins(unsigned char* block, long int k, int block_rows, int world_rank, int world_size);
void mergeMatrix(unsigned char** image, char* workblocks, long int k, int rows_per_block, int* displacements, int* blocksizes);

int main( int argc, char **argv )
{
    //Initialize MPI environment
    MPI_Init(&argc, &argv);
    //get number of processes
    double start, end;
    double fileread_start, fileread_end;
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    MPI_Barrier(MPI_COMM_WORLD);
    printf("Version 2.1\n");
    start = MPI_Wtime();


    //printf("0--------------\n");
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

    #if IMAGETEST != 1
    if (!initPlayground && world_rank != 0) {
        read_pgm_image( &maxval, &k, &k, file);
    }
    #endif
    //printf("1--------------\n");
    unsigned char* workBlocks = NULL; //This area of the memory will contain all the workblocks being passed to the processes with scatter
    unsigned char* image = NULL;
    int* blocksizes = NULL; //An array to store the size of every block of data passed to the process with scatter
    int* displacements = NULL; //An array to store the displacements of the blocks 
    MPI_Datatype row;
    MPI_Type_contiguous(k, MPI_UNSIGNED_CHAR, &row);
    MPI_Type_commit(&row);
    int block_rows = (world_rank >= k % world_size) ? (k/world_size + 2) : (k/world_size + 3); //Every process will have it's block to store the information
    unsigned char* block = malloc((block_rows * k)*sizeof( unsigned char));
    if (block == NULL) return -1;
    int sum = 0;
    fileread_start = MPI_Wtime();  
    if(!initPlayground) {
        if((workBlocks = malloc(((k/world_size + 2) * k)*world_size*sizeof(unsigned char) + (k % world_size) * k * sizeof(unsigned char)))==NULL) return -1;
    }
    if((blocksizes = malloc(world_size*sizeof(int))) == NULL) return -1;
    if((displacements = malloc(world_size*sizeof(int))) == NULL) return -1;
    for (int i = 0; i < world_size; i++) {
        //If the dimension of the grid is not divisible by the number of processes, we have to distribute the rest
        blocksizes[i] = (i >= k % world_size) ? (k/world_size + 2) : (k/world_size + 3); //In this way we spread the rest of the division around the first processes
        displacements[i] = sum; //We set the correct value for the displacement
        sum += blocksizes[i];
    }
    if (world_rank == 0 ) {
        if (!initPlayground) {
            #if IMAGETEST == 1 
                if((image = calloc(k * k, sizeof(unsigned char*))) == NULL) return -1;
                for (int i = 0; i < k-4; i+=7) {
                    image[(i+0) * k + 2] = MAXVAL;
                    image[(i+1) * k + 2] = MAXVAL;
                    image[(i+2) * k + 2] = MAXVAL;
                    image[(i+2) * k + 1] = MAXVAL;
                    image[(i+1) * k + 0] = MAXVAL;
                }
            #else
                image = read_pgm_image( &maxval, &k, &k, file);
            #endif
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
        }
    }
    
    if (initPlayground) {
        for (int i = 1; i < block_rows-1; i++) 
            for (int j= 0; j < k; j++) 
                block[i*k + j] = rand() % 2 == 0? MAXVAL : 0;
        updateMargins(block, k, block_rows, world_rank, world_size);              
    }
    else {
        MPI_Scatterv(workBlocks, blocksizes, displacements, row, block, block_rows, row, 0, MPI_COMM_WORLD);
    }
    //printf("3\n");
    

    //printf("I am process %d, with corners: %d, %d, %d, %d and with external margins: upper(%d,%d), lower(%d,%d)\n", world_rank, block[k], block[k*2 - 1], block[(block_rows - 1) * k - 1], block[(block_rows - 2) * k], block[0], block[k-1], block[block_rows*k-1], block[(block_rows - 1)*k]); //DEBUG
    fileread_end = MPI_Wtime();  
    unsigned char** writebuffer;

    if (world_rank== 0) 
            if((writebuffer = (unsigned char**) malloc(k *sizeof(unsigned char*)))==NULL)return -1;
    if (e) {

        unsigned char* tmp = NULL;
        unsigned char* blockB = malloc(((block_rows) * k)*sizeof(unsigned char));
        //printf("3.1\n");
        if (blockB == NULL) return -1;
        //printf("3.2\n");
        for (int i = 0; i < n; i++) {
            tmp = block;
            //printf("3.3\n");
            block = staticUpdate(block, block_rows, k, blockB);//not the most elegant way but is efficient
            blockB = tmp;
            //printf("3.4\n");
            updateMargins(block, k, block_rows, world_rank, world_size);
            //printf("3.5\n");
            if ( s > 0 && i % s == 0) {
                MPI_Gatherv(block, block_rows, row, workBlocks, blocksizes, displacements, row, 0, MPI_COMM_WORLD); //We need to recollect the original grid in order to save the file
                if (world_rank== 0) {
                    mergeMatrix(writebuffer, workBlocks, k, world_size, displacements, blocksizes);
                    write_pgm_image(writebuffer, MAXVAL, k, k, snapshotName(i));
                }
                MPI_Scatterv(workBlocks, blocksizes, displacements,  row, block, block_rows, row, 0, MPI_COMM_WORLD);
            }
        }
    }
    else {
        if (world_rank == 0) {
            for (int i = 0; i < k; i++)
                if ((writebuffer[i] = (unsigned char*)malloc(k * sizeof(unsigned char)))==NULL) return -1;
            for (int i = 0; i < k*k; i++) 
                writebuffer[i/k][i%k] = image[i];
                
            for (int i = 0; i < n; i++) {
                orderedUpdate(writebuffer, k, k);
                if (i%s==0)
                    write_pgm_image(writebuffer, MAXVAL, k, k, snapshotName(i)); 
            }
        }
    }
    double writetime_start;
    double writetime_end;
    writetime_start = MPI_Wtime();  
    /*for (int i = 0; i < world_size; i++) {
        if (world_rank == i) write_your_block(block, MAXVAL, k, "end.pgm", world_rank, block_rows);
        MPI_Barrier(MPI_COMM_WORLD);
    }*/
    //parallel_write(block, MAXVAL, k, "end.pgm", world_rank, block_rows, displacements, row);
    efficient_write(block, MAXVAL, k, "end.pgm", world_rank, world_size, block_rows, displacements, blocksizes, row);
    writetime_end = MPI_Wtime();
    //printf("4-----------------\n");
    /*MPI_Gatherv(block, block_rows, row, workBlocks, blocksizes, displacements, row, 0, MPI_COMM_WORLD);

    if (world_rank==0) {    
        mergetime_start = MPI_Wtime();    
        mergeMatrix(writebuffer, workBlocks, k, world_size, displacements, blocksizes);
        mergetime_end = MPI_Wtime();
        write_pgm_image(writebuffer, MAXVAL, k, k, "end.pgm");
    }
    //printf("5-------------\n");*/
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    if (world_rank == 0) {
        printf("Execution time on process %d: %f |", world_rank, end-start);
        printf(" Execution time for fileread making: %f |", fileread_end-fileread_start);
        printf(" Execution time for file writing (version 2): %f |", writetime_end-writetime_start);
    }
    MPI_Finalize();
    
    
    return 0;
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//--------------------------------------------------END OF MAIN-----------------------------------------------------------------
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


void mergeMatrix(unsigned char** image, char* workblocks, long int k, int world_size, int* displacements, int*blockSizes) {
    int imageIndex = 0;
    for (int i = 0; i<world_size; i++) {
        for (int j = displacements[i]+1; j < displacements[i] + blockSizes[i] - 1; j++) 
            image[imageIndex + (j - displacements[i]-1)]  = &workblocks[j*k];
        imageIndex += blockSizes[i]-2;
    } 
}

void updateMargins(unsigned char* block,long int k, int block_rows, int world_rank, int world_size) {
    //Here we preapre the rows to send and the buffer for the ones to receive
    unsigned char* lower_row_to_send = &block[(block_rows - 2) * k];
    unsigned char* upper_row_to_send = &block[k];
    MPI_Request req_lower_send;
    MPI_Status stat_lower_send;
    MPI_Request req_upper_send;
    MPI_Status stat_upper_send;
    MPI_Request req_lower_recv;
    MPI_Status stat_lower_recv;
    MPI_Request req_upper_recv;
    MPI_Status stat_upper_recv;

    MPI_Irecv(&block[(block_rows -1) * k] , k, MPI_UNSIGNED_CHAR, (world_rank + 1) % world_size, 1, MPI_COMM_WORLD, &req_lower_recv);
    MPI_Irecv(block, k, MPI_UNSIGNED_CHAR, world_rank == 0? world_size -1 : world_rank - 1, 0, MPI_COMM_WORLD, &req_upper_recv);
    MPI_Isend(lower_row_to_send, k, MPI_UNSIGNED_CHAR, (world_rank + 1) % world_size, 0, MPI_COMM_WORLD, &req_lower_send);
    MPI_Isend(upper_row_to_send, k, MPI_UNSIGNED_CHAR, world_rank == 0? world_size -1 : world_rank - 1, 1, MPI_COMM_WORLD, &req_upper_send);

    MPI_Wait(&req_upper_recv, &stat_upper_recv);
    MPI_Wait(&req_lower_recv, &stat_lower_recv);
    MPI_Wait(&req_upper_send, &stat_upper_send);
    MPI_Wait(&req_lower_send, &stat_lower_send);
}

void orderedUpdate(unsigned char** field, int k_i, int k) {
    int neighbours = 0;
    for ( int i = 0; i < k_i; i++ )
        for ( int j = 0; j < k; j++ ) {
            neighbours = ( field[(i > 0 ? i : k_i)-1][(j > 0? j : k)-1] + field[(i > 0 ? i : k_i)-1][j] + 
                        field[(i > 0 ? i : k_i)-1][(j+1)%k] +  field[i][(j+1)%k] +
                        field[(i+1)%k_i][(j+1)%k] + field[(i+1)%k_i][j] + 
                        field[(i+1)%k_i][(j > 0? j : k)-1] + field[i][(j > 0? j : k)-1] )/MAXVAL;
            //Here the ternary operator defines the toroid behauviour
            
            //The project rules are different from the original rules, so i inserted this directive to enable to switch between official rules and project rules
            #if PROJECT_RULES == 1
            if (neighbours == 2 || neighbours ==3)
                field[i][j] = MAXVAL;
            else
                field[i][j] = 0;
            #else
            if (((neighbours == 2 || neighbours == 3) && field[i][j] > 0 ) || (neighbours == 3 && field[i][j] == 0)) //We apply the three rules in one single check
                field[i][j] = MAXVAL;
            else
                field[i][j] = 0;
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
                    if (((neighbours == 2 || neighbours == 3) && block[i * k + j] > 0 ) || (neighbours == 3 && block[i * k + j] == 0)) //We apply the three rules in one single check
                        blockB[i * k + j] = MAXVAL;
                    else 
                        blockB[i * k + j] = 0;
                    #endif
                }
    }
    return blockB;
}

char* snapshotName(int n) {
    char* buffer;
    char* s = malloc(19*sizeof(char));

    if (n < 0)
        return NULL;
    else if (n < 10) {
        buffer = (char* )malloc(1*sizeof(char));
        strcpy(s, "snapshot_0000");
        strcat(s, intToString(n, buffer));
    }
    else if (n < 100) {
        buffer = (char* )malloc(2*sizeof(char));
        strcpy(s, "snapshot_000");
        strcat(s, intToString(n, buffer));
    }
    else if (n < 1000) {
        buffer = (char* )malloc(3*sizeof(char));
        strcpy(s, "snapshot_00");
        strcat(s, intToString(n, buffer));
    }
    else if (n < 10000) {
        buffer = (char* )malloc(4*sizeof(char));
        strcpy(s, "snapshot_0");
        strcat(s, intToString(n, buffer));
    }

    strcat(s, EXTENSION);

    return s;
}

char* intToString(int value, char* buffer) {
    int base = 10;
    
    int i = 0;
    while (value)
    {
        int r = value % base;
 
        if (r >= 10) {
            buffer[i++] = 65 + (r - 10);
        }
        else {
            buffer[i++] = 48 + r;
        }
 
        value = value / base;
    }
    // se il numero è 0
    if (i == 0) {
        buffer[i++] = '0';
    }
 
    // Se la base è 10 e il valore è negativo, la stringa risultante
    // è preceduto da un segno meno (-)
    // Con qualsiasi altra base, il valore è sempre considerato senza segno
    if (value < 0 && base == 10) {
        buffer[i++] = '-';
    }
    buffer[i] = '\0'; // stringa di terminazione nulla
    // inverte la stringa e la restituisce
    return reverse(buffer, 0, i - 1);
}

char* reverse(char *buffer, int i, int j)
{
    while (i < j) {
        swapc(&buffer[i++], &buffer[j--]);
    }
 
    return buffer;
}

void swapc(char *x, char *y) {
    char t = *x; *x = *y; *y = t;
}

/*void collectivePrint(unsigned char* block, int world_size, int world_rank, int block_rows, int k) {
    for (int k = 0; k < world_size; k++) {
        
        if (world_rank == k) {
            printf("Processo: %d\n", k);
            for ( int i = 1; i < block_rows - 1; i++ ) {
                for ( int j = 0; j < k; j++ )
                    if (block[i*k + j] > 0)
                        printf("1 ");
                    else
                        printf("0 ");
                printf("\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}*/