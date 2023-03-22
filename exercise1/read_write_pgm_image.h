#include <stdlib.h>
#include <stdio.h> 


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

void write_pgm_image( unsigned char **image, int maxval, int xsize, int ysize, const char *image_name);
unsigned char* read_pgm_image(int *maxval, int *xsize, int *ysize, const char *image_name);
void swap_image( void *image, int xsize, int ysize, int maxval );
void * generate_gradient( int maxval, int xsize, int ysize );