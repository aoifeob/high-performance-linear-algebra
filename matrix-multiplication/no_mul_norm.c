#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>

bool isLastThread(int threadNum, int totalThreads) {
    return threadNum == totalThreads - 1;
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

#pragma omp parallel shared (parallelMulResultMatrix, oneNorm) private (threadNumber, sliceWidth, elementsInSlice)
    {
        threadNumber = omp_get_thread_num();
        printf("Thread %d calculating norm\n\n", threadNumber);

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
    double matrix[4] = {23, 34, 31, 46};
    double parallelNorm;
    struct timeval tv1, tv2;
    struct timezone tz;
    int matrixDimension = 2; //default value, can be overwritten by user input
    int numThreads = 2;

    omp_set_num_threads(numThreads);

    // parallel matrix multiplication
    gettimeofday(&tv1, &tz);
    calculateParallelNorm(numThreads, matrixDimension, matrix, &parallelNorm);
    gettimeofday(&tv2, &tz);
    double parallelMulTimeElapsed =
            (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;

    // Serial matrix multiplication
    gettimeofday(&tv1, &tz);
    double serialNorm = calculateSerialNorm(matrixDimension, matrix);
    gettimeofday(&tv2, &tz);
    double serialMulTimeElapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;

    assertNormsAreEquivalent(serialNorm, parallelNorm);

    printf("Times taken for matrix multiplication on array with %dx%d dimensions: \n Serial: %f \n Parallel: %f\n\n",
           matrixDimension, matrixDimension, serialMulTimeElapsed, parallelMulTimeElapsed);

    return 0;
}