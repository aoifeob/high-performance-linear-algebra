#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <math.h>
#include <cblas.h>

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
    double leftMatrix[4] = {1, 2, 3, 4};
    double rightMatrix[4] = {5, 6, 7, 8};
    int matrixDimension = 2; //default value, can be overwritten by user input

    unsigned long matrixMemorySize = matrixDimension * matrixDimension * sizeof(double);
    double *serialMulResultMatrix = (double *) malloc(matrixMemorySize);
    double *blasMulResultMatrix = (double *) malloc(matrixMemorySize);

    if (!serialMulResultMatrix) {
        printf("Insufficient memory for matrices of dimension %d.\n", matrixDimension);
        exit(-1);
    }

    // Serial matrix multiplication
    serialMultiply(matrixDimension, leftMatrix, rightMatrix, serialMulResultMatrix);
    blasMultiply(matrixDimension, leftMatrix, rightMatrix, blasMulResultMatrix);

    assertMatricesAreEquivalent(matrixDimension, serialMulResultMatrix, blasMulResultMatrix);

    free(serialMulResultMatrix);

    return 0;
}
