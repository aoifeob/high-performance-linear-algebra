#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <math.h>
#include <cblas.h>

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

void blasMultiply(int matrixDimension, const double *leftMatrix, const double *rightMatrix,
                  double *blasMulResultMatrix){
    cblas_dgemm(CblasColMajor,
                CblasNoTrans,
                CblasNoTrans,
                matrixDimension,
                matrixDimension,
                matrixDimension,
                1.0,
                leftMatrix,
                matrixDimension,
                rightMatrix,
                matrixDimension,
                1.0,
                blasMulResultMatrix,
                matrixDimension);
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

void assertMatricesAreEquivalent(int matrixDimension, const double *serialMulResultMatrix,
                                 const double *parallelResultMatrix) {
    bool matrixValuesAreEqual = true;
    for (int col = 0; col < matrixDimension; col++) {
        for (int row = 0; row < matrixDimension; row++) {
            double serialMulElement = serialMulResultMatrix[col * matrixDimension + row];
            double parallelMulElement = parallelResultMatrix[col * matrixDimension + row];
            if (serialMulElement != parallelMulElement) {
                // print all non-matching values before exiting
                printf("Matrix elements at column %d, row %d are different. \n Serial mul matrix value: %f \n Parallel mul matrix value: %f \n\n",
                       col, row, serialMulElement, parallelMulElement);
                matrixValuesAreEqual = false;
            }
        }
    }
    if (!matrixValuesAreEqual) {
        exit(-1);
    }
}

int main(void) {
    double *leftMatrix, *rightMatrix;
    int matrixDimension = 2; //default value, can be overwritten by user input

    printf("Enter matrix dimension n : \n\n");
    scanf("%d", &matrixDimension);

    unsigned long matrixMemorySize = matrixDimension * matrixDimension * sizeof(double);
    leftMatrix = (double *)malloc(matrixMemorySize);
    rightMatrix = (double *)malloc(matrixMemorySize);
    double *serialMulResultMatrix = (double *) malloc(matrixMemorySize);
    double *blasMulResultMatrix = (double *) malloc(matrixMemorySize);

    if (!leftMatrix || !rightMatrix || !serialMulResultMatrix) {
        printf("Insufficient memory for matrices of dimension %d.\n", matrixDimension);
        exit(-1);
    }

    initMatrix(matrixDimension, leftMatrix);
    initMatrix(matrixDimension, rightMatrix);

    // Serial matrix multiplication
    serialMultiply(matrixDimension, leftMatrix, rightMatrix, serialMulResultMatrix);
    blasMultiply(matrixDimension, leftMatrix, rightMatrix, blasMulResultMatrix);

    assertMatricesAreEquivalent(matrixDimension, serialMulResultMatrix, blasMulResultMatrix);

    printf("Okay\n");

    free(leftMatrix);
    free(rightMatrix);
    free(serialMulResultMatrix);

    return 0;
}
