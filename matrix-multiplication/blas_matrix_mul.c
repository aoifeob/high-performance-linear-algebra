#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <cblas.h>

void initMatrix(int matrixDimension, double matrix[]) {
    for (int i = 0; i < matrixDimension; i++) {
        for (int j = 0; j < matrixDimension; j++) {
            double value = rand();
            matrix[i * matrixDimension + j] = value;
        }
    }
}

int main() {
    double *firstMatrix;
    double *secondMatrix;
    double *resultMatrix;
    struct timeval tv1, tv2;
    struct timezone tz;
    int matrixDimension = 1024; //default value, can be overwritten by user input

    printf("Enter matrix dimension n : \n\n");
    scanf("%d", &matrixDimension);

    unsigned long matrixMemorySize = matrixDimension * matrixDimension * sizeof(double);

    firstMatrix = malloc(matrixMemorySize);
    secondMatrix = malloc(matrixMemorySize);
    resultMatrix = malloc(matrixMemorySize);

    if (!firstMatrix || !secondMatrix || !resultMatrix) {
        printf("Insufficient memory for matrices of dimension %d.\n", matrixDimension);
        exit(-1);
    }

    initMatrix(matrixDimension, firstMatrix);
    initMatrix(matrixDimension, secondMatrix);

    gettimeofday(&tv1, &tz);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, matrixDimension, matrixDimension, matrixDimension, 1.0,
                firstMatrix, matrixDimension, secondMatrix, matrixDimension, 1.0, resultMatrix, matrixDimension);
    gettimeofday(&tv2, &tz);
    double timeElapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;
    printf("Time taken for simple matrix multiplication on array with %dx%d dimensions: %f.\n\n",
           matrixDimension, matrixDimension, timeElapsed);

    free(firstMatrix);
    free(secondMatrix);
    free(resultMatrix);

    return 0;
}