#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <cblas.h>

int main(void) {
    double *resultMatrix;
    struct timeval tv1, tv2;
    struct timezone tz;
    int matrixDimension = 2;
    int iterations = 10;
    double executionTimes[iterations];

    unsigned long matrixMemorySize = matrixDimension * matrixDimension * sizeof(double);

    resultMatrix = malloc(matrixMemorySize);

    if (!resultMatrix) {
        printf("Insufficient memory for matrices of dimension %d.\n", matrixDimension);
        exit(-1);
    }

    double a[2], b[2];
    gettimeofday(&tv1, &tz);

    for (int k = 0; k < 4; k++) {

        if (k == 0){
            a[0] = 1;
            a[1] = 1;
            b[0] = 1;
            b[1] = 1;
        }
        if (k == 1){
            a[0] = 1;
            a[1] = 1;
            b[0] = 0;
            b[1] = 0;
        }
        if (k == 3){
            a[0] = 0;
            a[1] = 0;
            b[0] = 1;
            b[1] = 1;
        }
        if (k == 4){
            a[0] = 0;
            a[1] = 0;
            b[0] = 0;
            b[1] = 0;
        }

        ATL_dgemm(CblasNoTrans,
                  CblasNoTrans,
                  1, //rows in A, C
                  1, //cols in B, C
                  2, //cols in A, rows in B
                  1.0,
                  a,
                  1, //size of first dimension of A
                  b,
                  1, //size of first dimension of B
                  1.0,
                  &resultMatrix[k],
                  1); //size of first dimension of C

        printf("Result matrix value at index %d is: %f\n\n", k, resultMatrix[k]);

    }

    gettimeofday(&tv2, &tz);
    double timeElapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;

    free(resultMatrix);

    printf("Average time taken for simple matrix multiplication on array with %dx%d dimensions: %f.\n\n",
           matrixDimension, matrixDimension, timeElapsed);

    return 0;
}