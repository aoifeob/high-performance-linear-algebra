#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>

bool isLastThread(int threadNum, int totalThreads) {
    return threadNum == totalThreads - 1;
}

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

void parallelMultiply(int numProcesses, int matrixDimension, const double *leftMatrix, double *rightMatrix,
                      double *parallelResultMatrix) {
    int threadNumber;
    int sliceWidth = matrixDimension / numProcesses;
    int elementsInSlice = matrixDimension * sliceWidth;

    //create parallel region
#pragma omp parallel shared (sliceWidth, elementsInSlice) private (threadNumber)
    {
        int numThreads = omp_get_num_threads();
        threadNumber = omp_get_thread_num();

        int resultMatrixStartingIndex = threadNumber * elementsInSlice;
        double *rightMatrixSlice = rightMatrix + threadNumber * elementsInSlice;
        int thisThreadSliceWidth = isLastThread(threadNumber, numThreads)
                                   ? matrixDimension - (sliceWidth * (numThreads - 1))
                                   : sliceWidth;

        int currentResultIndex = 0;

#pragma omp for
        //calculate result matrix elements for each slice
        for (int rightMatrixCol = 0; rightMatrixCol < matrixDimension; rightMatrixCol += sliceWidth) {

            for (int colInSlice = 0; colInSlice < thisThreadSliceWidth; colInSlice++) {

                for (int leftMatrixRow = 0; leftMatrixRow < matrixDimension; leftMatrixRow++) {

                    //calculate value for a single element of the result matrix by multiplying the single row and single column
                    double element = 0;
                    for (int k = 0; k < matrixDimension; k++) {
                        element += leftMatrix[leftMatrixRow + k * matrixDimension] *
                                   rightMatrixSlice[k + matrixDimension * colInSlice];
                    }

                    parallelResultMatrix[resultMatrixStartingIndex + currentResultIndex] = element;

                    currentResultIndex++;
                }
            }
        }
    } //end parallel region

}

void calculateParallelNorm(int numProcesses, int matrixDimension, double *parallelMulResultMatrix, double *oneNorm) {
    int threadNumber;
    double sliceNorm = 0;
    int sliceWidth = matrixDimension / numProcesses;
    int elementsInSlice = matrixDimension * sliceWidth;

#pragma omp parallel shared (sliceWidth, elementsInSlice) private (threadNumber, sliceNorm)
    {
        threadNumber = omp_get_thread_num();

        //construct thread data
        int thisThreadSliceWidth = isLastThread(threadNumber, numProcesses)
                                   ? matrixDimension - (sliceWidth * (numProcesses - 1))
                                   : sliceWidth;
        double *resultMatrixSlice = parallelMulResultMatrix + threadNumber * matrixDimension * thisThreadSliceWidth;

#pragma omp for
        //iterate through columns of the slice
        for (int resultMatrixCols = 0; resultMatrixCols < matrixDimension; resultMatrixCols += sliceWidth) {
            for (int colInSlice = 0; colInSlice < thisThreadSliceWidth; colInSlice++) {
                double thisColNorm = 0;

                //iterate through rows of the column to sum absolute values
                for (int row = 0; row < matrixDimension; row++) {
                    thisColNorm += fabs(resultMatrixSlice[colInSlice * matrixDimension + row]);
                }

                //if norm of the current column is greater than the current max column norm, update it to the current value
                if (thisColNorm > sliceNorm) {
                    sliceNorm = thisColNorm;
                }

            }

            //if norm of the current column is greater than the current max column norm, update it to the current value
#pragma omp critical
            if (sliceNorm > *(oneNorm)) {
                *(oneNorm) = sliceNorm;
            }
        }
    } //end parallel region
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
    int maxThreads = omp_get_num_procs();

    double *leftMatrix, *rightMatrix;
    double *parallelMulResultMatrix;
    double parallelNorm;
    struct timeval tv1, tv2;
    struct timezone tz;
    int matrixDimension = 2048; //default value, can be overwritten by user input
    int numThreads = 8;
    int shouldRunSerialProgram = 0;

    printf("This program supports serial and parallel one-norm computation. To disable serial computation, enter 1. Otherwise, enter 0.\n\n");
    scanf("%d", &shouldRunSerialProgram);

    printf("Enter matrix dimension n : \n\n");
    scanf("%d", &matrixDimension);

//    printf("Enter number of working threads p: \n\n");
//    if (scanf("%d", &numThreads) < 1 || numThreads > maxThreads) {
//        printf("Invalid number of threads %d specified", numThreads);
//        exit(-1);
//    }
//
//    if (numThreads > matrixDimension) {
//        printf("Number of threads p: %d should be smaller than matrix dimension n: %d\n\n",
//               matrixDimension, numThreads);
//        exit(-1);
//    }
//
//    if (0 != matrixDimension % numThreads) {
//        printf("Matrix with dimension n: %d and number of threads p: %d will be partitioned into uneven slices\n\n",
//               matrixDimension, numThreads);
//    }

    omp_set_num_threads(numThreads);

    unsigned long matrixMemorySize = matrixDimension * matrixDimension * sizeof(double);

    leftMatrix = (double *) malloc(matrixMemorySize);
    rightMatrix = (double *) malloc(matrixMemorySize);
    parallelMulResultMatrix = (double *) malloc(matrixMemorySize);

    if (!leftMatrix || !rightMatrix || !parallelMulResultMatrix) {
        printf("Insufficient memory for matrices of dimension %d.\n", matrixDimension);
        exit(-1);
    }

    initMatrix(matrixDimension, leftMatrix);
    initMatrix(matrixDimension, rightMatrix);

    // parallel matrix multiplication
    gettimeofday(&tv1, &tz);
    parallelMultiply(numThreads, matrixDimension, leftMatrix, rightMatrix, parallelMulResultMatrix);
    calculateParallelNorm(numThreads, matrixDimension, parallelMulResultMatrix, &parallelNorm);
    gettimeofday(&tv2, &tz);
    double parallelMulTimeElapsed =
            (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;

    if (!shouldRunSerialProgram) {
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
        assertNormsAreEquivalent(serialNorm, parallelNorm);

        printf("Times taken for matrix multiplication on array with %dx%d dimensions: \n Serial: %f \n Parallel: %f\n\n",
               matrixDimension, matrixDimension, serialMulTimeElapsed, parallelMulTimeElapsed);

        free(serialMulResultMatrix);
    } else {
        printf("Time taken for parallel matrix multiplication on array with %dx%d dimensions: %f\n\n",
               matrixDimension, matrixDimension, parallelMulTimeElapsed);
    }

    free(leftMatrix);
    free(rightMatrix);
    free(parallelMulResultMatrix);

    return 0;
}