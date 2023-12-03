#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
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

        printf("Thread %d has starting index %d and slice width %d\n ", threadNumber, resultMatrixStartingIndex,
               thisThreadSliceWidth);

#pragma omp for
        //calculate results for each slice
        for (int rightMatrixCol = 0; rightMatrixCol < matrixDimension; rightMatrixCol += sliceWidth) {

            for (int colInSlice = 0; colInSlice < thisThreadSliceWidth; colInSlice++) {

                printf("Thread %d of %d beginning execution of loop with slice width %d\n ", threadNumber, numThreads,
                       thisThreadSliceWidth);

                for (int leftMatrixRow = 0; leftMatrixRow < matrixDimension; leftMatrixRow++) {

                    //calculate value for a single element of the result matrix by multiplying the single row and single column
                    double element = 0;
                    for (int k = 0; k < matrixDimension; k++) {
                        element += leftMatrix[leftMatrixRow + k * matrixDimension] *
                                   rightMatrixSlice[k + matrixDimension * colInSlice];
                        printf("Thread %d "
                               "calculating result index %d "
                               "by summing product %f of\n "
                               "Left matrix index %d"
                               ": %f\n "
                               "Right matrix index %d"
                               ": %f\n\n",
                               threadNumber,
                               resultMatrixStartingIndex + currentResultIndex,
                               leftMatrix[leftMatrixRow + k * matrixDimension] * rightMatrixSlice[k + matrixDimension * rightMatrixCol],
                               leftMatrixRow + k * matrixDimension,
                               leftMatrix[leftMatrixRow + k * matrixDimension],
                               k + matrixDimension * colInSlice,
                               rightMatrixSlice[k + matrixDimension * colInSlice]);
                    }

                    parallelResultMatrix[resultMatrixStartingIndex + currentResultIndex] = element;

                    currentResultIndex++;
                }
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

int main(void) {
    int maxThreads = omp_get_num_procs();

    double leftMatrix[4] = {1, 2, 3, 4};
    double rightMatrix[4] = {5, 6, 7, 8};
    double *parallelMulResultMatrix;
    struct timeval tv1, tv2;
    struct timezone tz;
    int matrixDimension = 2; //default value, can be overwritten by user input
    int numThreads = 2;
    int shouldRunSerialProgram = 0;

//    printf("This program supports serial and parallel one-norm computation. To disable serial computation, enter 1. Otherwise, enter 0.\n\n");
//    scanf("%d", &shouldRunSerialProgram);
//
//    printf("Enter matrix dimension n : \n\n");
//    scanf("%d", &matrixDimension);
//
//    printf("Enter number of working processes p: \n\n");
//    if (scanf("%d", &numThreads) < 1 || numThreads > maxThreads) {
//        printf("Invalid number of processes %d specified", numThreads);
//        exit(-1);
//    }
//
//    if (numThreads > matrixDimension) {
//        printf("Number of processes p: %d should be smaller than matrix dimension n: %d\n\n",
//               matrixDimension, numThreads);
//        exit(-1);
//    }
//
//    if (0 != matrixDimension % numThreads) {
//        printf("Matrix with dimension n: %d and number of processes p: %d will be partitioned into uneven slices\n\n",
//               matrixDimension, numThreads);
//    }

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
        gettimeofday(&tv2, &tz);
        double serialMulTimeElapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;

        assertMatricesAreEquivalent(matrixDimension, serialMulResultMatrix, parallelMulResultMatrix);

        print(matrixDimension, serialMulResultMatrix);

        printf("Times taken for matrix multiplication on array with %dx%d dimensions: \n Serial: %f \n Parallel: %f\n\n",
               matrixDimension, matrixDimension, serialMulTimeElapsed, parallelMulTimeElapsed);

        free(serialMulResultMatrix);
    } else {
        printf("Time taken for parallel matrix multiplication on array with %dx%d dimensions: %f\n\n",
               matrixDimension, matrixDimension, parallelMulTimeElapsed);
    }

    free(parallelMulResultMatrix);

    return 0;
}