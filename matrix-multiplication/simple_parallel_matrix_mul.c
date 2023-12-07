#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <mpi.h>

int n, p;

int main(int argc, char **argv) {
    double leftMatrix[4] = {1, 2, 3, 4};
    double rightMatrix[4] = {5, 6, 7, 8};
    int squareSize;
    int myrank = 0;
    double *a, *b, *c, *allB, start, sum;
    n = 2;
    p = 4;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (fmod(sqrt(p), n) != 0) {
        printf("Invalid matrix size: %d for number of processors: %d. Sqrt p must be divisible by n. \n", n, p);
        exit(-1);
    }

    squareSize = n / sqrt(p);

    a = malloc(squareSize * n * sizeof(double));
    b = malloc(squareSize * n * sizeof(double));
    c = malloc(squareSize * n * sizeof(double));
    allB = malloc(n * n * sizeof(double));

    for (int i = 0; i < squareSize * squareSize; i++) {
        a[i] = leftMatrix[i * myrank];
        b[i] = rightMatrix[i * myrank];
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (myrank == 0) {
        start = MPI_Wtime();
    }
    for (int i = 0; i < p; i++) {
        MPI_Gather(b, squareSize * n, MPI_DOUBLE, allB, squareSize * n, MPI_DOUBLE,
                   i, MPI_COMM_WORLD);
    }

    for (int i = 0; i < squareSize; i++) {
        for (int j = 0; j < n; j++) {
            sum = 0.;
            for (int k = 0; k < n; k++) {
                sum += a[i * n + k] * allB[k * n + j];
            }
            c[i * n + j] = sum;
        }
    }

    free(allB);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myrank == 0) {
        printf("It took %f seconds to multiply 2 %dx%d matrices.\n",
               MPI_Wtime() - start, n, n);
    }

    MPI_Finalize();

    free(a);
    free(b);
    free(c);
}
