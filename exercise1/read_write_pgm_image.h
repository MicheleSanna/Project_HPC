#include <stdlib.h>
#include <stdio.h> 
#include <mpi.h>
#include <string.h>

#define XWIDTH 256
#define YWIDTH 256
#define MAXVAL 255


#if ((0x100 & 0xf) == 0x0)
#define I_M_LITTLE_ENDIAN 1
#define swap(mem) (( (mem) & (short int)0xff00) >> 8) +	\
  ( ((mem) & (short int)0x00ff) << 8)
#else
#define I_M_LITTLE_ENDIAN 0
#define swap(mem) (mem)
#endif

void efficient_write(unsigned char *block, int maxval, long int k, const char* image_name, int world_rank, int world_size, int block_rows,const int *displacements,const int* blocksizes,const MPI_Datatype row);
void efficient_write2(const unsigned char *block, int maxval, long int k, const char* image_name, int world_rank, int world_size, int block_rows,const int *displacements,const int* blocksizes,const MPI_Datatype row);
void parallel_write(const unsigned char *block, int maxval, long int k, const char* image_name, int world_rank, int block_rows,const int *displacements,const MPI_Datatype row);
void write_pgm_image(const unsigned char *image, int maxval,long int k, const char *image_name);
void write_your_block(const unsigned char *block, int maxval, long int x, const char* image_name, int world_rank, int block_rows);

unsigned char* read_pgm_image(int *maxval,long int *xsize,long int *ysize, const char *image_name);
unsigned char* parallel_read(int* maxval, long int* k, const char* image_name, int world_rank, int world_size, int* block_rows, int* displacements, int* blocksizes, MPI_Datatype* row);
unsigned char* read_your_segment(int* maxval, long int* x, const char* image_name, int world_size, int world_rank);

void makeBlockSizes(int k, int world_size, int *displacements, int *blocksizes);
void swap_image( void *image, int xsize, int ysize, int maxval );
void * generate_gradient( int maxval, int xsize, int ysize );