#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

void print(int matrixDimension, const double *matrix) {
    for (int col = 0; col < matrixDimension; col++) {
        for (int row = 0; row < matrixDimension; row++) {
            printf("%f\t", matrix[col * matrixDimension + row]);
        }
        printf("\n");
    }
    printf("\n");
}

void serialMultiply(int matrixDimension, const double *leftMatrix, const double *rightMatrix,
                    double *serialMulResultMatrix) {
    for (int col = 0; col < matrixDimension; col++) {
        for (int row = 0; row < matrixDimension; row++) {
            double element = 0;
            for (int k = 0; k < matrixDimension; k++) {

                element += leftMatrix[row + k * matrixDimension] * rightMatrix[k + matrixDimension * col];
            }
            serialMulResultMatrix[col * matrixDimension + row] = element;
        }
    }
}

double calculateSerialNorm(int matrixDimension, double *serialMulResultMatrix) {
    double oneNorm = 0;
    for (int col = 0; col < matrixDimension; col++) {
        double thisColNorm = 0;
        for (int row = 0; row < matrixDimension; row++) {
            double absoluteElementValue = fabs(serialMulResultMatrix[col * matrixDimension + row]);
            thisColNorm += absoluteElementValue;
        }
        if (thisColNorm > oneNorm) {
            oneNorm = thisColNorm;
        }
    }
    return oneNorm;
}

int main(void) {
    double leftMatrix[4] = {1, 2, 3, 4};
    double rightMatrix[4] = {5, 6, 7, 8};
    struct timeval tv1, tv2;
    struct timezone tz;
    int matrixDimension = 2; //default value, can be overwritten by user input

    unsigned long matrixMemorySize = matrixDimension * matrixDimension * sizeof(double);
    double *serialMulResultMatrix = (double *) malloc(matrixMemorySize);

    if (!serialMulResultMatrix) {
        printf("Insufficient memory for matrices of dimension %d.\n", matrixDimension);
        exit(-1);
    }

    // Serial matrix multiplication
    gettimeofday(&tv1, &tz);
    serialMultiply(matrixDimension, leftMatrix, rightMatrix, serialMulResultMatrix);
    double serialNorm = calculateSerialNorm(matrixDimension, serialMulResultMatrix);
    gettimeofday(&tv2, &tz);
    double serialMulTimeElapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;

    print(matrixDimension, serialMulResultMatrix);

    printf("Matrix norm %f\n\n", serialNorm);

    printf("Times taken for matrix multiplication on array with %dx%d dimensions: \n Serial: %f\n\n",
           matrixDimension, matrixDimension, serialMulTimeElapsed);

    free(serialMulResultMatrix);

    return 0;
}