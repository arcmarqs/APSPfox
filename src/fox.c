#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

#define MIN(A, B) (A < B) ? A : B
#define INF (INT_MAX / 2 - 1)

void checkFox(int p, int rows) {
    int q = sqrt(p);
    if (q * q == p && rows % q == 0) return;
    fprintf(stderr, "Fox algorithm can't be applied with a matrix of size %d and %d processes.\nAborting...\n", rows, p);
    MPI_Abort(MPI_COMM_WORLD, 0);
    exit(1);
}

void replaceBy(int** mat, int N, int old_value, int new_value) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i != j && mat[i][j] == old_value)
                mat[i][j] = new_value;
        }
    }
}

typedef struct __fox_grid {
    int nprocs;         // number of process
    int rows;           // process rows
    int cols;           // process cols
    int rank;           // process number
    int q;              // Sqrt of number of processes
    MPI_Comm comm;      // global comunicator
    MPI_Comm row_comm;  // row    comunicator
    MPI_Comm col_comm;  // columm comunicator
} FOX_GRID;

void Setup_grid(FOX_GRID* f) {
    int dim[2];
    int periods[2];
    int coords[2];
    int varying_coords[2];

    MPI_Comm_size(MPI_COMM_WORLD, &(f->nprocs));

    f->q = (int)sqrt((double)f->nprocs);
    dim[0] = dim[1] = f->q;
    periods[0] = periods[1] = 1;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, periods, 1, &(f->comm));
    MPI_Comm_rank(f->comm, &(f->rank));
    MPI_Cart_coords(f->comm, f->rank, 2, coords);
    f->rows = coords[0];
    f->cols = coords[1];

    varying_coords[0] = 0;
    varying_coords[1] = 1;
    MPI_Cart_sub(f->comm, varying_coords, &(f->row_comm));

    varying_coords[0] = 1;
    varying_coords[1] = 0;
    MPI_Cart_sub(f->comm, varying_coords, &(f->col_comm));
}

void computeScatter(int* sendcount, int* displace, int numprocs, int n_bar) {
    int gridsize = (int)sqrt(numprocs);

    for (int i = 0; i < numprocs; i++) sendcount[i] = 1;
    int disp = 0;
    for (int i = 0; i < gridsize; i++) {
        for (int j = 0; j < gridsize; j++) {
            displace[i * gridsize + j] = disp;
            disp++;
        }
        disp += ((n_bar - 1) * gridsize);
    }
}

void matrixMultiply(int** mat_A, int** mat_B, int** mat_Out, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                mat_Out[i][j] = MIN(mat_Out[i][j], mat_A[i][k] + mat_B[k][j]);
            }
        }
    }
}

void Fox(FOX_GRID* f, int** mat_A, int** mat_B, int** mat_Out, int** tmp_A, int n_bar) {
    int source = (f->rows + 1) % f->q;
    int dest = (f->rows + f->q - 1) % f->q;

    int bcast_root;
    for (int step = 0; step < f->q; step++) {
        bcast_root = (f->rows + step) % f->q;
        if (bcast_root == f->cols) {
            MPI_Bcast(ARRAY(mat_A), n_bar * n_bar, MPI_INT, bcast_root, f->row_comm);
            matrixMultiply(mat_A, mat_B, mat_Out, n_bar);
        } else {
            MPI_Bcast(ARRAY(tmp_A), n_bar * n_bar, MPI_INT, bcast_root, f->row_comm);
            matrixMultiply(tmp_A, mat_B, mat_Out, n_bar);
        }

        MPI_Status status;
        MPI_Sendrecv_replace(ARRAY(mat_B), n_bar * n_bar, MPI_INT, dest, 1, source, 1, f->col_comm, &status);
    }
}

void APSP(FOX_GRID* grid, int** mat_Out, int n_bar) {
    int** mat_A = buildMatrix(n_bar);
    int** mat_B = buildMatrix(n_bar);
    int** tmp_A = buildMatrix(n_bar);

    for (int i = 1; i < n_bar * grid->q - 1; i <<= 1) {  // i *= 2 --> i <<= 1

        copyMatrix(mat_Out, mat_A, n_bar);
        copyMatrix(mat_Out, mat_B, n_bar);

        Fox(grid, mat_A, mat_B, mat_Out, tmp_A, n_bar);
    }

    destroyMatrix(mat_A);
    destroyMatrix(mat_B);
    destroyMatrix(tmp_A);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int numprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int N;
    int** input_matrix;
    if (rank == 0) {
        scanf("%d", &N);  // read n
        input_matrix = readMatrix(N);
    }

#ifdef TIME
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();
#endif

    if (rank == 0) {
        checkFox(numprocs, N);
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    FOX_GRID grid;
    Setup_grid(&grid);

    int n_bar = N / grid.q;

    int* globalptr = NULL;
    if (rank == 0) {
        replaceBy(input_matrix, N, 0, INF);
        globalptr = ARRAY(input_matrix);
    }

    int sizes[2] = {N, N};
    int subsize[2] = {n_bar, n_bar};
    int starts[2] = {0, 0};

    MPI_Datatype type, subarraytype;
    MPI_Type_create_subarray(2, sizes, subsize, starts, MPI_ORDER_C, MPI_INT, &type);
    MPI_Type_create_resized(type, 0, n_bar * sizeof(int), &subarraytype);
    MPI_Type_commit(&subarraytype);

    int sendcount[numprocs];
    int displace[numprocs];
    if (rank == 0) {
        computeScatter(sendcount, displace, numprocs, n_bar);
    }

    int** mat_Out = buildMatrix(n_bar);
    MPI_Scatterv(globalptr, sendcount, displace, subarraytype, ARRAY(mat_Out), n_bar * n_bar, MPI_INT, 0, grid.comm);

    APSP(&grid, mat_Out, n_bar);

    MPI_Gatherv(ARRAY(mat_Out), n_bar * n_bar, MPI_INT, globalptr, sendcount, displace, subarraytype, 0, grid.comm);

    destroyMatrix(mat_Out);

    MPI_Type_free(&subarraytype);

    if (rank == 0) {
        replaceBy(input_matrix, N, INF, 0);
    }

#ifdef TIME
    MPI_Barrier(MPI_COMM_WORLD);
    double finish = MPI_Wtime();
#endif

    if (rank == 0) {
        printMatrix(input_matrix, N);

#ifdef TIME
        printf("Execution time: %f seconds\n", finish - start);
#endif
    }

    MPI_Finalize();
}
