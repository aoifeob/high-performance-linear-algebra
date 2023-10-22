#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

void initMatrix(int matrixDimension, double matrix[]) {
//    printf("Beginning matrix initialization.\n");
    for (int i = 0; i < matrixDimension; i++) {
        for (int j = 0; j < matrixDimension; j++) {
            double value = rand();
//            printf("Inserting value %f at row %d, column %d. \n", value, i, j);
            matrix[i * matrixDimension + j] = value;
        }
    }
//    printf("Matrix initialization complete.\n\n");
}

void initZeroMatrix(int matrixDimension, double matrix[]) {
//    printf("Beginning result matrix initialization.\n");
    for (int i = 0; i < matrixDimension; i++) {
        for (int j = 0; j < matrixDimension; j++) {
            matrix[i * matrixDimension + j] = 0;
        }
    }
//    printf("Matrix result initialization complete.\n\n");
}

void multiply(int matrixDimension, double *firstMatrix, double *secondMatrix, double *resultMatrix) {
    for (int i = 0; i < matrixDimension; i++) {
        for (int j = 0; j < matrixDimension; j++) {
            for (int k = 0; k < matrixDimension; k++) {
                int resultMatrixIndex = i * matrixDimension + j;

                printf("i=%d, j=%d, k=%d\n", i, j, k);

                printf("Result matrix index %d value : %f" , resultMatrixIndex, resultMatrix[i * matrixDimension + j]);

                resultMatrix[i * matrixDimension + j] +=
                        firstMatrix[i * matrixDimension + k] * secondMatrix[k * matrixDimension + j];

                printf("updated to : %f\n" , resultMatrix[i * matrixDimension + j]);

            }
//            printf("Result matrix has value %f at row %d, column %d. \n", resultMatrix[i * matrixDimension + j], i, j);
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
//        printf("Memory allocation required for each matrix is %lu bytes.\n\n", matrixMemorySize);

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
//        printf("Time taken for iteration %d of simple matrix multiplication on array with %dx%d dimensions: %f.\n\n", j,
//               matrixDimension, matrixDimension, timeElapsed);
        initZeroMatrix(matrixDimension, resultMatrix);
    }

    for (int j=0; j < 9; j++){
//        printf("Result matrix value at index %d is: %f\n\n", j, resultMatrix[j]);
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