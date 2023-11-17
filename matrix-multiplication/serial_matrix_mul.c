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

void serialMultiply(int matrixDimension, const double *leftMatrix, const double *rightMatrix, double *serialMulResultMatrix) {
    for (int col = 0; col < matrixDimension; col++) {
        for (int row = 0; row < matrixDimension; row++) {
            for (int k = 0; k < matrixDimension; k++) {
                serialMulResultMatrix[row * matrixDimension + col] +=
                        leftMatrix[col + k * matrixDimension] * rightMatrix[k + matrixDimension * row];
            }
        }
    }
}

double calculateSerialNorm(int matrixDimension, double *serialMulResultMatrix){
    double oneNorm = 0;
    for (int col=0; col < matrixDimension ; col++){
        double thisColNorm = 0;
        for (int row=0; row < matrixDimension ; row++){
            double absoluteElementValue = fabs(serialMulResultMatrix[col * matrixDimension + row]);
            thisColNorm += absoluteElementValue;
        }
        if (thisColNorm > oneNorm){
            oneNorm = thisColNorm;
        }
    }
    return oneNorm;
}

int main(void) {
    double leftMatrix[] = {1, 2, 3, 4};
    double rightMatrix[] = {5, 6, 7, 8};
    double *serialMulResultMatrix;
    double serialNorm;
    struct timeval tv1, tv2;
    struct timezone tz;
    int matrixDimension = 2; //default value, can be overwritten by user input
    int iterations = 1;
    double serialMulExecutionTimes[iterations];

    unsigned long matrixMemorySize = matrixDimension * matrixDimension * sizeof(double);

    for (int i = 0; i < iterations; i++) {
        serialMulResultMatrix = malloc(matrixMemorySize);

        if (!serialMulResultMatrix) {
            printf("Insufficient memory for matrices of dimension %d.\n", matrixDimension);
            exit(-1);
        }

        // Serial matrix multiplication
        gettimeofday(&tv1, &tz);
        serialMultiply(matrixDimension, leftMatrix, rightMatrix, serialMulResultMatrix);
        serialNorm = calculateSerialNorm(matrixDimension, serialMulResultMatrix);
        gettimeofday(&tv2, &tz);
        double serialMulTimeElapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;
        serialMulExecutionTimes[i] = serialMulTimeElapsed;

        printf("Serial norm for matrix is %f\n", serialNorm);

        free(serialMulResultMatrix);
    }

    double avgSerialMulIterationTime = 0;
    for (int i = 0; i < iterations; i++) {
        avgSerialMulIterationTime += serialMulExecutionTimes[i];
    }
    avgSerialMulIterationTime = avgSerialMulIterationTime / iterations;

    printf("Average time taken for matrix norm calculation on matrix with %dx%d dimensions: \n Serial: %f \n\n",
           matrixDimension, matrixDimension, avgSerialMulIterationTime);

    return 0;
}
