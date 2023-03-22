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
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <stdbool.h>
#include <string.h>
#include <mpi.h>
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

void orderedUpdate(unsigned char** field, int k_i, int k_j);
unsigned char* staticUpdate(unsigned char* block, int k_i, int k_j, unsigned char* blockB);
char* snapshotName(int n);
char* intToString(int value, char* buffer);
char* reverse(char *buffer, int i, int j);
void swapc(char *x, char *y);
void closestFactors(int nOfProcess, int* n_i, int* n_j);
void updateMargins(unsigned char* block, int k_j, int block_rows, int world_rank, int world_size);
void mergeMatrix(unsigned char** image, char* workblocks, int k_i, int k_j, int block_rows);

int main( int argc, char **argv )
{

    //printf("AAAAAAAAA 0\n"); //DEBUG
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
    start = MPI_Wtime();

    //printf("AAAAAAAAA 1\n"); //DEBUG
    
    int maxval = MAXVAL;
    bool initPlayground = 1; //Flag to state if is a new game or a loaded one. Default is new
    int k = 0; //Size of playground, must be specified
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


    if (!initPlayground && world_rank != 0) {
        read_pgm_image( &maxval, &k, &k, file);
    }

    unsigned char* workBlocks = NULL; //This area of the memory will contain all the workblocks being passed to the processes with scatter
    unsigned char* image = NULL;
    int* blocksizes = NULL; //An array to store the size of every block of data passed to the process with scatter
    int* displacements = NULL; //An array to store the displacements of the blocks 
    if (world_rank == 0) {     

              
        if (!initPlayground) { 
            image = read_pgm_image( &maxval, &k, &k, file);
        }
        else {
            if((image = calloc(k * k, sizeof(unsigned char*))) == NULL) return -1;
        }
    
        
        if (initPlayground) {
            for (int i = 0; i < k-4; i+=7) {
                image[(i+0) * k + 2] = MAXVAL;
                image[(i+1) * k + 2] = MAXVAL;
                image[(i+2) * k + 2] = MAXVAL;
                image[(i+2) * k + 1] = MAXVAL;
                image[(i+1) * k + 0] = MAXVAL;
            }
        } //I draw a column of gliders in case of a new playground, to be able to check if the algorithm is working
        
        
        int sum = 0;
        if((workBlocks = malloc(((k/world_size + 2) * k)*world_size*sizeof(unsigned char) + (k % world_size) * k * sizeof(unsigned char)))==NULL) return -1;
        if((blocksizes = malloc(world_size*sizeof(int))) == NULL) return -1;
        if((displacements = malloc(world_size*sizeof(int))) == NULL) return -1;
        for (int i = 0; i < world_size; i++) {
            //If the dimension of the grid is not divisible by the number of processes, we have to distribute the rest
            blocksizes[i] = (i >= k % world_size) ? ((k/world_size + 2) * k) : ((k/world_size + 3) * k); //In this way we spread the rest of the division around the first processes
            displacements[i] = sum; //We set the correct value for the displacement
            sum += blocksizes[i];
        }

        int workBlocksIndex = 0;

        for (int blockCount = 0; blockCount < world_size; blockCount++) //We fill the workblocks vector block by block
            for (int i = displacements[blockCount]/k -1 -(2*blockCount); i < displacements[blockCount]/k -1 -(2*blockCount) + blocksizes[blockCount]/k; i++) //We are attaching the external rows contiguosly to the subgrid, so we start 1 row before and finish 1 row after
                for (int j = 0; j < k; j++) {
                    workBlocks[workBlocksIndex] = image[((i >= 0)? (i%k) : (k -1)) * k + j]; //Here the ternary operator defines the toroid behauviour
                    workBlocksIndex++;
                }
            
        

    }
    int block_rows= (world_rank >= k % world_size) ? (k/world_size + 2) : (k/world_size + 3); //Every process will have it's block to store the information 
    unsigned char* block = malloc((block_rows * k)*sizeof( unsigned char));
    MPI_Scatterv(workBlocks, blocksizes, displacements, MPI_UNSIGNED_CHAR, block, block_rows* k, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
    

    //printf("I am process %d, with corners: %d, %d, %d, %d and with external margins: upper(%d,%d), lower(%d,%d)\n", world_rank, block[k_j], block[k_j*2 - 1], block[(block_rows - 1) * k_j - 1], block[(block_rows - 2) * k_j], block[0], block[k_j-1], block[block_rows*k_j-1], block[(block_rows - 1)*k_j]); //DEBUG
    
    unsigned char** writebuffer;
    if (world_rank== 0) 
            if((writebuffer = (unsigned char**) malloc(k *sizeof(unsigned char*)))==NULL)return -1;
    if (e) {

        unsigned char* tmp = NULL;
        unsigned char* blockB = malloc(((block_rows) * k)*sizeof(unsigned char));
        if (blockB == NULL) return -1;
        for (int i = 0; i < n; i++) {
            tmp = block;
            block = staticUpdate(block, block_rows, k, blockB);//not the most elegant way but is efficient
            blockB = tmp;
            updateMargins(block, k, block_rows, world_rank, world_size);
            if ( s > 0 && i % s == 0) {
                MPI_Gatherv(block, block_rows* k, MPI_UNSIGNED_CHAR, workBlocks, blocksizes, displacements, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD); //We need to recollect the original grid in order to save the file
                if (world_rank== 0) {
                    mergeMatrix(writebuffer, workBlocks, k, k, block_rows);
                    write_pgm_image(writebuffer, MAXVAL, k, k, snapshotName(i));
                }
                MPI_Scatterv(workBlocks, blocksizes, displacements, MPI_UNSIGNED_CHAR, block, block_rows* k, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
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
    //getchar();
    //printf("---END---"); //DEBUG
    //MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(block, block_rows* k, MPI_UNSIGNED_CHAR, workBlocks, blocksizes, displacements, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
    
    if (world_rank==0) {
        mergeMatrix(writebuffer, workBlocks, k, k, block_rows);
        write_pgm_image(writebuffer, MAXVAL, k, k, "end.pgm");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    printf(" Execution time on process %d: %f ", world_rank, end-start);
    MPI_Finalize();
    
    
    return 0;
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//--------------------------------------------------END OF MAIN-----------------------------------------------------------------
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


void mergeMatrix(unsigned char** image, char* workblocks, int k_i, int k_j, int block_rows) {
    int j = 1;
    for(int i = 0; i < k_i; i++) {
        if (j%block_rows == block_rows -1)
            j += 2;
        image[i] = &workblocks[j * k_j]; //That's the most efficient way i found to merge the matrix, simply point at the rows jumping the ones that are repeated due to message passing mechanism
        j++;
    }
    printf("j: %d", j);
}

void updateMargins(unsigned char* block, int k_j, int block_rows, int world_rank, int world_size) {
    //Here we preapre the rows to send and the buffer for the ones to receive
    unsigned char* lower_row_to_send = &block[(block_rows - 2) * k_j];
    unsigned char* upper_row_to_send = &block[k_j];
    unsigned char* lower_row_to_update = malloc(k_j * sizeof(unsigned char));
    unsigned char* upper_row_to_update = malloc(k_j * sizeof(unsigned char));


    if (world_rank == 0) {
        MPI_Send(lower_row_to_send, k_j, MPI_UNSIGNED_CHAR, 1, 0, MPI_COMM_WORLD); //Process 0 does the send first to avoid the deadlock
        MPI_Recv(upper_row_to_update, k_j, MPI_UNSIGNED_CHAR, world_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } 
    else {
        MPI_Recv(upper_row_to_update, k_j, MPI_UNSIGNED_CHAR, world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //receive
        MPI_Send(lower_row_to_send, k_j, MPI_UNSIGNED_CHAR, (world_rank + 1) % world_size, 0, MPI_COMM_WORLD);//and send
    }
    for (int i = 0; i < k_j; i++) //Update values
        block[i] = upper_row_to_update[i];
    free(upper_row_to_update);


    if (world_rank == 0) {
        MPI_Send(upper_row_to_send, k_j, MPI_UNSIGNED_CHAR, world_size - 1, 1, MPI_COMM_WORLD); //Process 0 does the send first to avoid the deadlock
        MPI_Recv(lower_row_to_update, k_j, MPI_UNSIGNED_CHAR, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   
    } 
    else {
        MPI_Recv(lower_row_to_update, k_j, MPI_UNSIGNED_CHAR, (world_rank + 1) % world_size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//receive
        MPI_Send(upper_row_to_send, k_j, MPI_UNSIGNED_CHAR, world_rank - 1, 1, MPI_COMM_WORLD);//and send
    }
    for (int i = 0; i < k_j; i++)  //Update values
        block[(block_rows -1) * k_j + i] = lower_row_to_update[i];
    free(lower_row_to_update);
}

void orderedUpdate(unsigned char** field, int k_i, int k_j) {
    int neighbours = 0;
    for ( int i = 0; i < k_i; i++ )
        for ( int j = 0; j < k_j; j++ ) {
            neighbours = ( field[(i > 0 ? i : k_i)-1][(j > 0? j : k_j)-1] + field[(i > 0 ? i : k_i)-1][j] + 
                        field[(i > 0 ? i : k_i)-1][(j+1)%k_j] +  field[i][(j+1)%k_j] +
                        field[(i+1)%k_i][(j+1)%k_j] + field[(i+1)%k_i][j] + 
                        field[(i+1)%k_i][(j > 0? j : k_j)-1] + field[i][(j > 0? j : k_j)-1] )/MAXVAL;
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
    printf("Debug 4\n");
    return;
}
unsigned char* staticUpdate(unsigned char* block, int block_rows, int k_j, unsigned char* blockB) {
    #pragma omp parallel 
    {
        int neighbours = 0;
        #pragma omp for collapse(2) 
            for ( int i = 1; i < block_rows - 1; i++ )
                for ( int j = 0; j < k_j; j++ ) {
                    neighbours = (block[(i - 1) * k_j + (j > 0? j : k_j)-1] + block[(i - 1) * k_j + j] + 
                                block[(i - 1) * k_j + (j+1)%k_j] +  block[i * k_j + (j+1)%k_j] +
                                block[(i + 1) * k_j + (j+1)%k_j] + block[(i + 1) * k_j + j] + 
                                block[(i + 1) * k_j + (j > 0? j : k_j)-1] + block[i * k_j + (j > 0? j : k_j)-1])/MAXVAL;
                    //Here the ternary operator defines the toroid behauviour

                    //The project rules are different from the original rules, so i inserted this directive to enable to switch between official rules and project rules
                    #if PROJECT_RULES == 1 
                    if (neighbours == 2 || neighbours ==3)
                        blockB[i * k_j + j] = MAXVAL;
                    else
                        blockB[i * k_j + j] = 0;
                    #else
                    if (((neighbours == 2 || neighbours == 3) && block[i * k_j + j] > 0 ) || (neighbours == 3 && block[i * k_j + j] == 0)) //We apply the three rules in one single check
                        blockB[i * k_j + j] = MAXVAL;
                    else 
                        blockB[i * k_j + j] = 0;
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
    // if the number is 0
    if (i == 0) {
        buffer[i++] = '0';
    }
 
    // if the base is 10 and the value negative, the string
    // is preceded by a (-) sign
    // With any other base, the value is considered without sign
    if (value < 0 && base == 10) {
        buffer[i++] = '-';
    }
    buffer[i] = '\0'; // Null termination string
    // Reversing the string and returning it
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
 
