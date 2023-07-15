#include "read_write_pgm_image.h"

// =============================================================
//  utilities for managinf pgm files
//
//  * write_pgm_image
//  * read_pgm_image
//  * swap_image
//
// =============================================================


void efficient_write(unsigned char *block, int maxval, long int k, const char* image_name, int world_rank, int world_size, int block_rows, const int *displacements, const int* blocksizes, const MPI_Datatype row) {
  MPI_Status status;
  MPI_Request req;
  MPI_File fh;
  int offset;
  int written = 1;
  char* str = malloc(70*sizeof(char));

  sprintf(str, "P5\n# generated by\n# Michele\n%d %d\n%d\n", k, k, maxval);
  offset = strlen(str) + 1;
  if (world_rank == 0) {
    MPI_File_open(MPI_COMM_SELF, image_name, MPI_MODE_WRONLY | MPI_MODE_CREATE ,MPI_INFO_NULL, &fh);
    MPI_File_set_size(fh, (k*k+offset)*sizeof(char));
    MPI_File_write_at(fh, 0, str, offset, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_write_at(fh, offset-1, &block[k], block_rows - 2, row, MPI_STATUS_IGNORE);
    
    while (written < world_size) {
      MPI_Recv(block,  (block_rows-2), row, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
      MPI_File_write_at(fh, offset-1+((displacements[status.MPI_SOURCE]-(2*status.MPI_SOURCE))*k), block, blocksizes[status.MPI_SOURCE] -2, row, MPI_STATUS_IGNORE);
      written++;
    }

    MPI_File_close(&fh);
  }
  else {
    MPI_Isend(&block[k],  (block_rows-2), row, 0, 0, MPI_COMM_WORLD, &req);
    MPI_Request_free(&req);
  }
}

void efficient_write2(const unsigned char *block, int maxval, long int k, const char* image_name, int world_rank, int world_size, int block_rows, const int *displacements, const int* blocksizes, const MPI_Datatype row) {
  MPI_Status status;
  MPI_Request req;
  MPI_File fh;
  int offset;
  int written = 1;
  char* str = malloc(70*sizeof(char));

  sprintf(str, "P5\n# generated by\n# Michele\n%d %d\n%d\n", k, k, maxval);
  offset = strlen(str) + 1;
  if (world_rank == 0) {
    unsigned char* buffer;
    if ((buffer =  malloc(k*block_rows*sizeof(unsigned char)))==NULL) return;
    for (long int i = 0; i < k*block_rows; i++) buffer[i] = block[i];

    MPI_File_open(MPI_COMM_SELF, image_name, MPI_MODE_WRONLY | MPI_MODE_CREATE ,MPI_INFO_NULL, &fh);
    MPI_File_set_size(fh, (k*k+offset)*sizeof(char));
    MPI_File_write_at(fh, 0, str, offset, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_write_at(fh, offset-1, &buffer[k], block_rows - 2, row, MPI_STATUS_IGNORE);

    while (written < world_size) {
      MPI_Recv(buffer,  (block_rows-2), row, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
      MPI_File_write_at(fh, offset-1+((displacements[status.MPI_SOURCE]-(2*status.MPI_SOURCE))*k), buffer, blocksizes[status.MPI_SOURCE] -2, row, MPI_STATUS_IGNORE);
      written++;
    }

    MPI_File_close(&fh);
  }
  else {
    MPI_Isend(&block[k],  (block_rows-2), row, 0, 0, MPI_COMM_WORLD, &req);
    MPI_Request_free(&req);
  }
}

void parallel_write(const unsigned char *block, int maxval, long int k, const char* image_name, int world_rank, int block_rows, const int *displacements, const MPI_Datatype row) {
  MPI_File fh;
  int offset;
  char* str = malloc(70*sizeof(char));
  sprintf(str, "P5\n# generated by\n# Michele\n%d %d\n%d\n", k, k, maxval);
  offset = strlen(str) + 1;

  MPI_File_open(MPI_COMM_WORLD, image_name, MPI_MODE_WRONLY | MPI_MODE_CREATE ,MPI_INFO_NULL, &fh);

  if (world_rank == 0) MPI_File_write_at(fh, 0, str, offset, MPI_CHAR, MPI_STATUS_IGNORE);

  MPI_File_write_at(fh, offset+((displacements[world_rank]-(2*world_rank))*k), &block[k], block_rows - 2, row, MPI_STATUS_IGNORE);

  MPI_File_close(&fh);
}

void write_your_block(const unsigned char *block, int maxval, long int k, const char* image_name, int world_rank, int block_rows) {
  int color_depth = 1 + ( maxval > 255 );
  FILE* image_file;
  if (world_rank == 0) { 
    image_file = fopen(image_name, "w"); 
    fprintf(image_file, "P5\n# generated by\n# Michele\n%d %d\n%d\n", k, k, maxval);  
    fclose(image_file);
  }
  if (k < 30000) {  
      image_file = fopen(image_name, "a"); 
      for (int i = 1; i < block_rows-1; i++) {
        fwrite(&block[i*k], 1, k*color_depth, image_file);  
      }
      fclose(image_file);
  } 
  else {
    for(int i = 1; i < block_rows-1; i++) {
      image_file = fopen(image_name, "a"); 
      fwrite(&block[i*k], 1, k*color_depth, image_file); 
      fclose(image_file);
    }
  }
  return ;


}

void write_pgm_image( const unsigned char *image, int maxval, long int k, const char *image_name) {

  int color_depth = 1 + ( maxval > 255 );
  FILE* image_file;
  image_file = fopen(image_name, "w"); 
  fprintf(image_file, "P5\n# generated by\n# Michele\n%d %d\n%d\n", k, k, maxval); 
  if (k < 30000) {  
      // Writing file
      for (int i = 0; i < k; i++) {
        fwrite(&image[i*k], 1, k*color_depth, image_file);  
      }
      fclose(image_file);
  } 
  else {
    fclose(image_file);
    for(int i = 0; i < k; i++) {
      image_file = fopen(image_name, "a"); 
      fwrite(&image[i*k], 1, k*color_depth, image_file); 
      fclose(image_file);
    }
  }
  return ;

  /* ---------------------------------------------------------------
     TYPE    MAGIC NUM     EXTENSION   COLOR RANGE
           ASCII  BINARY
     PBM   P1     P4       .pbm        [0-1]
     PGM   P2     P5       .pgm        [0-255]
     PPM   P3     P6       .ppm        [0-2^16[
  
  ------------------------------------------------------------------ */
}

unsigned char* parallel_read(int* maxval, long int* k, const char* image_name, int world_rank, int world_size, int* block_rows, int* displacements, int* blocksizes, MPI_Datatype* row) {
  MPI_File fh;
  FILE* image_file; 
  unsigned char* block;
  *k = *maxval = 0;
  int sum;
  char MagicN[2];
  char *line = NULL;
  size_t  m, n = 0;

  image_file = fopen(image_name, "r"); 
  if (image_file == NULL) return NULL;
  
  // get the Magic Number
  m = fscanf(image_file, "%2s%*c", MagicN );

  // skip all the comments
  m = getline( &line, &n, image_file);
  while ( (m > 0) && (line[0]=='#') )
    m = getline( &line, &n, image_file);
  if (m > 0) {
    m = sscanf(line, "%d%*c%d%*c%d%*c", k, k, maxval);
    if ( m < 3 )
	    fscanf(image_file, "%d%*c", maxval);
    MPI_Type_contiguous(*k, MPI_UNSIGNED_CHAR, row);
    MPI_Type_commit(row);
  }
  else return NULL;
  
  free( line );
  
  int color_depth = 1 + ( *maxval > 255 );
   
  size_t offset = ftell(image_file);

  fclose(image_file);

  makeBlockSizes(*k, world_size, displacements, blocksizes);
  *block_rows = blocksizes[world_rank];

  if ( (block = malloc((*block_rows) * (*k) * sizeof(unsigned char))) == NULL ) return NULL;

  MPI_File_open(MPI_COMM_WORLD, image_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  MPI_File_read_at(fh, offset + ((displacements[world_rank]-(2*world_rank))*(*k)), &block[*k], (*block_rows)-2,*row, MPI_STATUS_IGNORE);
  MPI_File_close(&fh);
  return block;
}

unsigned char* read_pgm_image( int *maxval,long int *xsize,long int *ysize, const char *image_name)
/*
 * image        : a pointer to the pointer that will contain the image
 * maxval       : a pointer to the int that will store the maximum intensity in the image
 * xsize, ysize : pointers to the x and y sizes
 * image_name   : the name of the file to be read
 *
 */
{
  FILE* image_file; 
  unsigned char* image;
  image_file = fopen(image_name, "r"); 
  *xsize = *ysize = *maxval = 0;
  if (image_file == NULL)
    return NULL;
  char MagicN[2];
  char *line = NULL;
  size_t  k, n = 0;
  // get the Magic Number
  k = fscanf(image_file, "%2s%*c", MagicN );

  // skip all the comments
  k = getline( &line, &n, image_file);
  while ( (k > 0) && (line[0]=='#') )
    k = getline( &line, &n, image_file);
  if (k > 0)
  {
    k = sscanf(line, "%d%*c%d%*c%d%*c", xsize, ysize, maxval);
    if ( k < 3 )
	    fscanf(image_file, "%d%*c", maxval);
  }
  else
  {
    *maxval = -1;         // this is the signal that there was an I/O error
			    // while reading the image header
    free( line );
    return NULL;
  }
  free( line );
  
  int color_depth = 1 + ( *maxval > 255 );
  long int size = *xsize * *ysize * color_depth;
  
  if ( (image = malloc( size )) == NULL )
    {
      fclose(image_file);
      *maxval = -2;         // this is the signal that memory was insufficient
      *xsize  = 0;
      *ysize  = 0;
      return NULL;
    }
  for (int i=0;i< *ysize;i++) {
    if ( fread(&image[i * *xsize], 1, *xsize, image_file) != *xsize)
      {
        free( image );
        image   = NULL;
        *maxval = -3;         // this is the signal that there was an i/o error
        *xsize  = 0;
        *ysize  = 0;
      } 
  } 
  fclose(image_file);
  return image;
}


void swap_image( void *image, int xsize, int ysize, int maxval )
/*
 * This routine swaps the endianism of the memory area pointed
 * to by ptr, by blocks of 2 bytes
 *
 */
{
  if ( maxval > 255 )
    {
      // pgm files has the short int written in
      // big endian;
      // here we swap the content of the image from
      // one to another
      //
      unsigned int size = xsize * ysize;
      for ( int i = 0; i < size; i++ )
  	((unsigned short int*)image)[i] = swap(((unsigned short int*)image)[i]);
    }
  return;
}



// =============================================================
//

void * generate_gradient( int maxval, int xsize, int ysize )
/*
 * just and example about how to generate a vertical gradient
 * maxval is either 255 or 65536, xsize and ysize are the
 * x and y dimensions of the image to be generated.
 * The memory region that will contain the image is returned
 * by the function as a void *
 */
{
  char      *cImage;   // the image when a single byte is used for each pixel
  short int *sImage;   // the image when a two bytes are used for each pixel
  void      *ptr;
  
  int minval      = 0; 
  int delta       = (maxval - minval) / ysize;
  
  if(delta < 1 )
    delta = 1;
  
  if( maxval < 256 )
    // generate a gradient with 1 byte of color depth
    {
      cImage = (char*)calloc( xsize*ysize, sizeof(char) );
      unsigned char _maxval = (char)maxval;
      int idx = 0;
      for ( int yy = 0; yy < ysize; yy++ )
	{
	  unsigned char value = minval + yy*delta;
	  for( int xx = 0; xx < xsize; xx++ )
	    cImage[idx++] = (value > _maxval)?_maxval:value;
	}
      ptr = (void*)cImage;
    }
  else
    // generate a gradient with 2 bytes of color depth
    {
      sImage = (unsigned short int*)calloc( xsize*ysize, sizeof(short int) );
      unsigned short int _maxval = swap((unsigned short int)maxval);
      int idx = 0;
      for ( int yy = 0; yy < ysize; yy++ )
	{
	  unsigned short int value  = (short int) (minval+ yy*delta);
	  unsigned short int _value = swap( value );    // swap high and low bytes, the format expect big-endianism
	  
	  for( int xx = 0; xx < xsize; xx++ )
	    sImage[idx++] = (value > maxval)?_maxval:_value;
	}
      ptr = (void*)sImage;	
    }

  return ptr;
}


void makeBlockSizes(int k, int world_size, int* displacements, int* blocksizes) {
    int sum = 0;

    for (int i = 0; i < world_size; i++) {
        //If the dimension of the grid is not divisible by the number of processes, we have to distribute the rest
        blocksizes[i] = (i >= k % world_size) ? (k/world_size + 2) : (k/world_size + 3); //In this way we spread the rest of the division around the first processes
        displacements[i] = sum; //We set the correct value for the displacement
        sum += blocksizes[i];
    }
}


