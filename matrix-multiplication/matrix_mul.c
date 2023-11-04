#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

void initMatrix(int matrixDimension, double matrix[]) {
    for (int i = 0; i < matrixDimension; i++) {
        for (int j = 0; j < matrixDimension; j++) {
            double value = rand();
            matrix[i * matrixDimension + j] = value;
        }
    }
}

void initZeroMatrix(int matrixDimension, double matrix[]) {
    for (int i = 0; i < matrixDimension; i++) {
        for (int j = 0; j < matrixDimension; j++) {
            matrix[i * matrixDimension + j] = 0;
        }
    }
}

void multiply(int matrixDimension, double *firstMatrix, double *secondMatrix, double *resultMatrix) {
    for (int i = 0; i < matrixDimension; i++) {
        for (int j = 0; j < matrixDimension; j++) {
            float value = 0;
            for (int k = 0; k < matrixDimension; k++) {
                value += firstMatrix[i * matrixDimension + k] * secondMatrix[k * matrixDimension + j];
            }
            resultMatrix[i * matrixDimension + j] = value;
        }
    }
}

int main(void) {
    double *firstMatrix;
    double *secondMatrix;
    double *resultMatrix;
    struct timeval tv1, tv2;
    struct timezone tz;
    int matrixDimension = 3;
    int iterations = 10;
    double executionTimes[iterations];

    double testFirstMatrix[] = {9, 3, 8, 5, 7, 2, 1, 4, 6};
    double testSecondMatrix[] = {6, 5, 4, 3, 2, 1, 8, 9, 7};

    printf("Enter matrix dimension n : \n\n");
    scanf("%d", &matrixDimension);

    printf("Beginning operations using an %dx%d matrix.\n\n", matrixDimension, matrixDimension);

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

    for (int j = 0; j < iterations; j++) {
        gettimeofday(&tv1, &tz);
//        multiply(matrixDimension, firstMatrix, secondMatrix, resultMatrix);
        multiply(matrixDimension, testFirstMatrix, testSecondMatrix, resultMatrix);
        gettimeofday(&tv2, &tz);
        double timeElapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;
        executionTimes[j] = timeElapsed;
        initZeroMatrix(matrixDimension, resultMatrix);
    }

    free(firstMatrix);
    free(secondMatrix);
    free(resultMatrix);

    double avgExecutionTime = 0;
    for (int j = 0; j < iterations; j++) {
        avgExecutionTime += executionTimes[j];
    }
    avgExecutionTime = avgExecutionTime / iterations;
    printf("Average time taken for simple matrix multiplication on array with %dx%d dimensions: %f.\n\n",
           matrixDimension, matrixDimension, avgExecutionTime);

    return 0;
}