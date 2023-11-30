#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>

bool isLastThread(int threadNum, int totalThreads) {
    return threadNum == totalThreads - 1;
}

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

void *multiplySlice(int sliceWidth, int matrixDimension, int resultMatrixSliceStartingIndex, const double *leftMatrix,
                    const double *rightMatrixSlice, double *parallelResultMatrix) {
    int currentResultIndex = 0;

    for (int rightMatrixCol = 0; rightMatrixCol < sliceWidth; rightMatrixCol++) {
        for (int leftMatrixRow = 0; leftMatrixRow < matrixDimension; leftMatrixRow++) {

            //calculate value for a single element of the result matrix by multiplying the single row and single column
            double element = 0;
            for (int k = 0; k < matrixDimension; k++) {
                element += leftMatrix[leftMatrixRow + k * matrixDimension] *
                           rightMatrixSlice[k + matrixDimension * rightMatrixCol];
            }
            parallelResultMatrix[resultMatrixSliceStartingIndex + currentResultIndex] = element;

            currentResultIndex++;
        }
    }

}

void parallelMultiply(int numProcesses, int matrixDimension, double *leftMatrix, double *rightMatrix,
                      double *parallelResultMatrix) {
    int threadNumber;
    int sliceWidth = matrixDimension / numProcesses;
    int elementsInSlice = matrixDimension * sliceWidth;

    //create threads
#pragma omp parallel shared (leftMatrix, rightMatrix, parallelResultMatrix) private (threadNumber)
    {
        //construct thread data
        threadNumber = omp_get_thread_num();
        double *rightMatrixSlice = rightMatrix + threadNumber * elementsInSlice;
        int thisThreadSliceWidth = isLastThread(threadNumber, numProcesses)
                                   ? matrixDimension - (sliceWidth * (numProcesses - 1))
                                   : sliceWidth;
        int resultMatrixSliceStartingIndex = threadNumber * matrixDimension * sliceWidth;

        //calculate slice
        multiplySlice(thisThreadSliceWidth, matrixDimension, resultMatrixSliceStartingIndex,
                      leftMatrix, rightMatrixSlice, parallelResultMatrix);
    } //end parallel

}

void *calculateSliceNorm(int matrixDimension, int sliceWidth, const double *resultMatrixSlice, double *oneNorm) {
    double sliceNorm = 0;

    //iterate through columns of the slice
    for (int col = 0; col < sliceWidth; col++) {
        double thisColNorm = 0;

        //iterate through rows of the column to sum absolute values
        for (int row = 0; row < matrixDimension; row++) {
            double absoluteElementValue = fabs(
                    resultMatrixSlice[col * matrixDimension + row]);
            thisColNorm += absoluteElementValue;
        }

        //if norm of the current column is greater than the current max column norm, update it to the current value
        if (thisColNorm > sliceNorm) {
            sliceNorm = thisColNorm;
        }

    }

    //if norm of the current column is greater than the current max column norm, update it to the current value
    if (sliceNorm > *(oneNorm)) {
        *(oneNorm) = sliceNorm;
        printf("Updating one norm to %f\n\n", sliceNorm);
    }

}

void calculateParallelNorm(int numProcesses, int matrixDimension, double *parallelMulResultMatrix, double *oneNorm) {
    int threadNumber;
    int sliceWidth = matrixDimension / numProcesses;
    int elementsInSlice = matrixDimension * sliceWidth;

#pragma omp parallel shared (parallelMulResultMatrix, oneNorm, sliceWidth, elementsInSlice) private (threadNumber)
    {
        threadNumber = omp_get_thread_num();

        //construct thread data
        double *resultMatrixSlice = parallelMulResultMatrix + threadNumber * elementsInSlice;
        int thisThreadSliceWidth = isLastThread(threadNumber, numProcesses)
                                   ? matrixDimension - (sliceWidth * (numProcesses - 1))
                                   : sliceWidth;

        calculateSliceNorm(matrixDimension, thisThreadSliceWidth, resultMatrixSlice, oneNorm);
    }
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

void assertNormsAreEquivalent(double serialNorm, double parallelNorm) {
    if (serialNorm != parallelNorm) {
        printf("Matrix norms are different. Serial norm: %f \n Parallel norm: %f \n\n", serialNorm, parallelNorm);
        exit(-1);
    }
}

int main(void) {
    double leftMatrix[4] = {1, 2, 3, 4};
    double rightMatrix[4] = {5, 6, 7, 8};
    double *parallelMulResultMatrix;
    double parallelNorm;
    struct timeval tv1, tv2;
    struct timezone tz;
    int matrixDimension = 2; //default value, can be overwritten by user input
    int numThreads = 2;

    omp_set_num_threads(numThreads);

    unsigned long matrixMemorySize = matrixDimension * matrixDimension * sizeof(double);

    parallelMulResultMatrix = (double *) malloc(matrixMemorySize);

    if (!parallelMulResultMatrix) {
        printf("Insufficient memory for matrices of dimension %d.\n", matrixDimension);
        exit(-1);
    }

    // parallel matrix multiplication
    gettimeofday(&tv1, &tz);
    parallelMultiply(numThreads, matrixDimension, leftMatrix, rightMatrix, parallelMulResultMatrix);
    calculateParallelNorm(numThreads, matrixDimension, parallelMulResultMatrix, &parallelNorm);
    gettimeofday(&tv2, &tz);
    double parallelMulTimeElapsed =
            (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;

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

    assertMatricesAreEquivalent(matrixDimension, serialMulResultMatrix, parallelMulResultMatrix);
    print(matrixDimension, serialMulResultMatrix);
    assertNormsAreEquivalent(serialNorm, parallelNorm);

    printf("Times taken for matrix multiplication on array with %dx%d dimensions: \n Serial: %f \n Parallel: %f\n\n",
           matrixDimension, matrixDimension, serialMulTimeElapsed, parallelMulTimeElapsed);

    free(serialMulResultMatrix);
    free(parallelMulResultMatrix);

    return 0;
}