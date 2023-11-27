#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

void initMatrix(int matrixDimension, double matrix[]) {
    for (int col = 0; col < matrixDimension; col++) {
        for (int row = 0; row < matrixDimension; row++) {
            double value = rand();
            matrix[col * matrixDimension + row] = value;
        }
    }
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
    double *leftMatrix, *rightMatrix;
    double *serialMulResultMatrix;
    struct timeval tv1, tv2;
    struct timezone tz;
    int matrixDimension = 2048; //default value, can be overwritten by user input

    printf("Enter matrix dimension n : \n\n");
    scanf("%d", &matrixDimension);

    unsigned long matrixMemorySize = matrixDimension * matrixDimension * sizeof(double);

    leftMatrix = (double *) malloc(matrixMemorySize);
    rightMatrix = (double *) malloc(matrixMemorySize);
    serialMulResultMatrix = (double *) malloc(matrixMemorySize);

    if (!leftMatrix || !rightMatrix || !serialMulResultMatrix) {
        printf("Insufficient memory for matrices of dimension %d.\n", matrixDimension);
        exit(-1);
    }

    initMatrix(matrixDimension, leftMatrix);
    initMatrix(matrixDimension, rightMatrix);

    // Serial matrix multiplication
    gettimeofday(&tv1, &tz);
    serialMultiply(matrixDimension, leftMatrix, rightMatrix, serialMulResultMatrix);
    double serialNorm = calculateSerialNorm(matrixDimension, serialMulResultMatrix);
    gettimeofday(&tv2, &tz);
    double serialMulTimeElapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;

    printf("Times taken for matrix norm calculation on array with %dx%d dimensions: %f\n\n",
           matrixDimension, matrixDimension, serialMulTimeElapsed);

    free(serialMulResultMatrix);
    free(leftMatrix);
    free(rightMatrix);

    return 0;
}