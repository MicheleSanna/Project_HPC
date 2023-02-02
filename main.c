#if defined(__STDC__)
#  if (__STDC_VERSION__ >= 199901L)
#     define _XOPEN_SOURCE 700
#  endif
#endif
#if !defined(_OPENMP)
//#error "OpenMP support needed for this code"
#endif
#define PROJECT_RULES 0
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <stdbool.h>
#include "read_write_pgm_image.h"
#include <string.h>

void orderedUpdate(char** field, int k_i, int k_j, int s);
char** staticUpdate(char** field, int k_i, int k_j, int s, char** fieldB);
char* snapshotName(int n);
char* itoa(int value, char* buffer);
char* reverse(char *buffer, int i, int j);
void swapc(char *x, char *y);

int main( int argc, char **argv )
{
    int maxval = MAXVAL;
    bool initPlayground = 1; //Flag to state if is a new game or a loaded one. Default is new
    int k = 0; //Size of playground, must be specified
    int e = 1; //Update mode, default is 1 (static)
    char* file = NULL; //String containing file name, must be specified
    int n = 0; //Number of steps of the computation, must be specified else the program will do nothing
    int s = 0; //Number of steps to take a snapshot of the system. Default is 0 (no snapshot)
    printf("diocane1");
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
    printf("diocane2");
    int k_i = k, k_j = k; //For now those will be the same. maybe later i will differentiate them
    void** image = calloc(k_i * k_j, sizeof(bool)); 
    if (!initPlayground)
        read_pgm_image( image, &maxval, &k, &k, file);
    char** fieldA = (char**) image;
    printf("diocane3");
    //-------------------TEST CODE-----------------------
    fieldA[3][3] = 1;
    fieldA[3][4] = 1;
    fieldA[3][5] = 1;
    fieldA[2][5] = 1;
    fieldA[1][4] = 1;
    //----------------END OF TEST CODE--------------------
    printf("diocane4");
    int i;
    if (e) {
        char** fieldB = malloc(k_i * k_j * sizeof(bool)); 
        char** tmp = NULL;
        int i = 0;
        for (i = 0; i < n; i++) {
            tmp = fieldA;
            fieldA = staticUpdate(fieldA, k_i, k_j, s,fieldB);//Done this trick with pointers to avoid everytime the copy of a matrix, not the most elegant way but is efficient
            fieldB = tmp;
        }
    }
    else {
        for (i = 0; i < n; i++)
            orderedUpdate(fieldA, k_i, k_j, s);
    }
    printf("diocane5");
    write_pgm_image(fieldA, MAXVAL, k_i, k_j, file);
    printf("diocane6");
    getchar();
    return 0;
}




void orderedUpdate(char** field, int k_i, int k_j, int s) {
    int neighbours = 0;
    for ( int i = 0; i < k_i; i++ )
        for ( int j = 0; j < k_j; j++ ) {
            neighbours = field[i-1][j-1] + field[i-1][j] + field[i-1][j+1] + field[i][j+1] + field[i+1][j+1] + field[i+1][j] + field[i+1][j-1] + field[i][j-1];
            #if PROJECT_RULES == 1
            if (neighbours == 2 || neighbours ==3)
                field[i][j] = MAXVAL;
            else
                field[i][j] = 0;
            #else
            if (((neighbours == 2 || neighbours ==3) && field[i][j] == 1 ) || (neighbours == 3 && field[i][j] == 0))
                field[i][j] = MAXVAL;
            else
                field[i][j] = 0;
            #endif
            if (i%s==0)
                write_pgm_image(field, MAXVAL, k_i, k_j, snapshotName(i));

        }
    return;
}
char** staticUpdate(char** field, int k_i, int k_j, int s, char** fieldB) {
    int neighbours = 0;
    for ( int i = 0; i < k_i; i++ )
        for ( int j = 0; j < k_j; j++ ) {
            neighbours = field[i-1][j-1] + field[i-1][j] + field[i-1][j+1] + field[i][j+1] + field[i+1][j+1] + field[i+1][j] + field[i+1][j-1] + field[i][j-1];
            #if PROJECT_RULES == 1
            if (neighbours == 2 || neighbours ==3)
                fieldB[i][j] = MAXVAL;
            else
                fieldB[i][j] = 0;
            #else
            if (((neighbours == 2 || neighbours ==3) && field[i][j] == 1 ) || (neighbours == 3 && field[i][j] == 0))
                fieldB[i][j] = MAXVAL;
            else
                fieldB[i][j] = 0;
            #endif
            if (i%s==0)
                write_pgm_image(field, MAXVAL, k_i, k_j, snapshotName(i));
        }
    return fieldB;
}

char* snapshotName(int n) {
    char* buffer;

    if (n < 0)
        return NULL;
    if (n < 10) {
        buffer = (char* )malloc(1*sizeof(char));
        return strcat("snapshot_0000", itoa(n, buffer));
    }
    if (n < 100) {
        buffer = (char* )malloc(2*sizeof(char));
        return strcat("snapshot_000", itoa(n, buffer));
    }
    if (n < 1000) {
        buffer = (char* )malloc(3*sizeof(char));
        return strcat("snapshot_00", itoa(n, buffer));
    }
    if (n < 10000) {
        buffer = (char* )malloc(4*sizeof(char));
        return strcat("snapshot_0", itoa(n, buffer));
    }
}

char* itoa(int value, char* buffer) {
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
 