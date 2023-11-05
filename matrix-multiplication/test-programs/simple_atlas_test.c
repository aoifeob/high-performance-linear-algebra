#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>

int main(void) {
    double *resultMatrix;
    int matrixDimension = 2;

    double testFirstMatrix[] = {1, 1, 0, 0};
    double testSecondMatrix[] = {1, 0, 1, 0};

    unsigned long matrixMemorySize = matrixDimension * matrixDimension * sizeof(double);

    resultMatrix = malloc(matrixMemorySize);

    if (!resultMatrix) {
        printf("Insufficient memory for matrices of dimension %d.\n", matrixDimension);
        exit(-1);
    }

    ATL_dgemm(CblasNoTrans,
              CblasNoTrans,
              matrixDimension,
              matrixDimension,
              matrixDimension,
              1.0,
              testFirstMatrix,
              matrixDimension,
              testSecondMatrix,
              matrixDimension,
              1.0,
              resultMatrix,
              matrixDimension);

    for (int j = 0; j < matrixDimension * matrixDimension; j++) {
        printf("Result matrix value at index %d is: %f\n\n", j, resultMatrix[j]);
    }

    free(resultMatrix);

    return 0;
}