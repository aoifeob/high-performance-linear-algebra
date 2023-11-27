#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <math.h>

typedef struct {
    double *leftMatrix;
    double *rightMatrixSlice;
    double *resultMatrix;
    int matrixDimension;
    int sliceWidth;
    int resultMatrixSliceStartingIndex;
} mul_slice_data;

typedef struct {
    double *resultMatrixSlice;
    double *oneNorm;
    int matrixDimension;
    int sliceWidth;
} norm_slice_data;

bool isLastThread(int threadNum, int totalThreads){
    return threadNum == totalThreads -1;
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

void *multiplySlice(void *arg) {
    mul_slice_data *mul_slice_data = arg;

    int currentResultIndex = 0;

    for (int rightMatrixCol = 0; rightMatrixCol < mul_slice_data->sliceWidth; rightMatrixCol++) {
        for (int leftMatrixRow = 0; leftMatrixRow < mul_slice_data->matrixDimension; leftMatrixRow++) {

            //calculate value for a single element of the result matrix by multiplying the single row and single column
            double element = 0;
            for (int k = 0; k < mul_slice_data->matrixDimension; k++) {
                element += mul_slice_data->leftMatrix[leftMatrixRow + k * mul_slice_data->matrixDimension] *
                           mul_slice_data->rightMatrixSlice[k + mul_slice_data->matrixDimension * rightMatrixCol];
            }
            mul_slice_data->resultMatrix[mul_slice_data->resultMatrixSliceStartingIndex + currentResultIndex] = element;

            currentResultIndex++;
        }
    }

}

void parallelMultiply(int numProcesses, int matrixDimension, double *leftMatrix, double *rightMatrix,
                      double *parallelResultMatrix) {
    void *thread_status;
    mul_slice_data *thread_mul_slice_data;
    int sliceWidth = matrixDimension / numProcesses;
    int elementsInSlice = matrixDimension * sliceWidth;

    thread_mul_slice_data = malloc(numProcesses * sizeof(mul_slice_data));

    //create threads
    for (int thread = 0; thread < numProcesses; thread++) {
        //construct slice data
        thread_mul_slice_data[thread].leftMatrix = leftMatrix;
        thread_mul_slice_data[thread].rightMatrixSlice = rightMatrix + thread * elementsInSlice;
        thread_mul_slice_data[thread].resultMatrix = parallelResultMatrix;
        thread_mul_slice_data[thread].matrixDimension = matrixDimension;
        thread_mul_slice_data[thread].sliceWidth = isLastThread(thread, numProcesses) ? matrixDimension - (sliceWidth * (numProcesses-1)) : sliceWidth;
        thread_mul_slice_data[thread].resultMatrixSliceStartingIndex = thread * matrixDimension * sliceWidth;

        //calculate slice

    }

    free(thread_mul_slice_data);
}

void *calculateSliceNorm(void *arg) {
    norm_slice_data *norm_slice_data = arg;
    double sliceNorm = 0;

    //iterate through columns of the slice
    for (int col = 0; col < norm_slice_data->sliceWidth; col++) {
        double thisColNorm = 0;

        //iterate through rows of the column to sum absolute values
        for (int row = 0; row < norm_slice_data->matrixDimension; row++) {
            double absoluteElementValue = fabs(
                    norm_slice_data->resultMatrixSlice[col * norm_slice_data->matrixDimension + row]);
            thisColNorm += absoluteElementValue;
        }

        //if norm of the current column is greater than the current max column norm, update it to the current value
        if (thisColNorm > sliceNorm) {
            sliceNorm = thisColNorm;
        }

    }

    //if norm of the current column is greater than the current max column norm, update it to the current value
    if (sliceNorm > *(norm_slice_data->oneNorm)) {
        *(norm_slice_data->oneNorm) = sliceNorm;
    }

}

void calculateParallelNorm(int numProcesses, int matrixDimension, double *parallelMulResultMatrix, double *oneNorm) {
    norm_slice_data *thread_norm_slice_data;
    int sliceWidth = matrixDimension / numProcesses;
    int elementsInSlice = matrixDimension * sliceWidth;

    thread_norm_slice_data = malloc(numProcesses * sizeof(mul_slice_data));

    for (int thread = 0; thread < numProcesses; thread++) {
        //construct slice data
        thread_norm_slice_data[thread].resultMatrixSlice = parallelMulResultMatrix + thread * elementsInSlice;
        thread_norm_slice_data[thread].oneNorm = oneNorm;
        thread_norm_slice_data[thread].matrixDimension = matrixDimension;
        thread_norm_slice_data[thread].sliceWidth = isLastThread(thread, numProcesses) ? matrixDimension - (sliceWidth * (numProcesses-1)) : sliceWidth;

        //calculate norm for cols in slice
    }

    free(thread_norm_slice_data);
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
    int matrixDimension = 2048; //default value, can be overwritten by user input
    int numThreads = 8;
    int shouldRunSerialProgram = 0;

    printf("This program supports serial and parallel one-norm computation. To disable serial computation, enter 1. Otherwise, enter 0.\n\n");
    scanf("%d", &shouldRunSerialProgram);

    printf("Enter matrix dimension n : \n\n");
    scanf("%d", &matrixDimension);

    printf("Enter number of working processes p: \n\n");
    if (scanf("%d", &numThreads) < 1) {
        printf("Invalid number of processes %d specified", numThreads);
        exit(-1);
    }

    if (numThreads > matrixDimension) {
        printf("Number of processes p: %d should be smaller than matrix dimension n: %d\n\n",
               matrixDimension, numThreads);
        exit(-1);
    }

    if (0 != matrixDimension % numThreads) {
        printf("Matrix with dimension n: %d and number of processes p: %d will be partitioned into uneven slices\n\n",
               matrixDimension, numThreads);
    }

    unsigned long matrixMemorySize = matrixDimension * matrixDimension * sizeof(double);

    parallelMulResultMatrix = (double *)malloc(matrixMemorySize);

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

    if (!shouldRunSerialProgram){
        double *serialMulResultMatrix = (double *)malloc(matrixMemorySize);

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

    free(parallelMulResultMatrix);

    return 0;
}