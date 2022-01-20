#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

int** buildMatrix(int N) {
    /* allocate the n*m contiguous items */
    int *p = (int *)malloc(N * N *sizeof(int));

    /* allocate the row pointers into the memory */
    int** mat = (int **)malloc(N * sizeof(int *));

    /* set up the pointers into the contiguous memory */
    for (int i=0; i<N; i++)
       mat[i] = &(p[i * N]);

    return mat;
}

void destroyMatrix(int** mat) {
    /* free the memory - the first element of the array is at the start */
    free(ARRAY(mat));

    /* free the pointers into the memory */
    free(mat);

    mat = NULL;
}

int** readMatrix(int N) {
    int** mat = buildMatrix(N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            scanf("%d", &(mat[i][j]));
        }
    }

    return mat;
}

void fillMatrix(int** mat, int value, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            mat[i][j] = value;
        }
    }
}

void copyMatrix(int** source, int** dest, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            dest[i][j] = source[i][j];
        }
    }
}

void printMatrix(int** mat, int N) {
    for (int i = 0; i < N; i++) {
        printf("%d", mat[i][0]);
        for (int j = 1; j < N; j++) {
            printf(" %d", mat[i][j]);
        }
        printf("\n");
    }
}