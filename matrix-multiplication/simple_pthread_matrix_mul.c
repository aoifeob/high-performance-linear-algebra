#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <math.h>
#include <pthread.h>
//#include <cblas.h>

#define MAX_THREADS 124

typedef struct {
    double *leftMatrix;
    double *rightMatrixSlice;
    double *resultMatrix;
    int matrixDimension;
    int sliceWidth;
} mul_slice_data;

void serialMultiply(int matrixDimension, const double *leftMatrix, const double *rightMatrix,
                    double *serialMulResultMatrix) {
    for (int col = 0; col < matrixDimension; col++) {
        for (int row = 0; row < matrixDimension; row++) {
            double element = 0;
            for (int k = 0; k < matrixDimension; k++) {

                 element += leftMatrix[row + k * matrixDimension] * rightMatrix[k + matrixDimension * col];
            }
            printf("Setting value %f on result matrix at col %d row %d. Result stored in index %d\n\n",
                   element,
                   col,
                   row,
                   col * matrixDimension + row);
            serialMulResultMatrix[col * matrixDimension + row] = element;
            printf("Result matrix col %d row %d has value %f\n\n", col, row, serialMulResultMatrix[row * matrixDimension + col]);
        }
    }
}

void *pThreadMultiplySlice(void *arg) {
    //TODO: implement slice multiplication logic
    mul_slice_data *mul_slice_data = arg;
    int matrixDimension = mul_slice_data->matrixDimension;

    int resultMatrixIndex = 0;

    for (int rightMatrixCol = 0; rightMatrixCol < mul_slice_data->sliceWidth; rightMatrixCol++) {
        for (int leftMatrixRow = 0; leftMatrixRow < matrixDimension; leftMatrixRow++) {
            for (int leftMatrixCol = 0; leftMatrixCol < matrixDimension; leftMatrixCol++) {
                //TODO: create slice of left matrix containing a single row
                double *leftMatrixSingleRowSlice;
                leftMatrixSingleRowSlice = malloc(mul_slice_data->matrixDimension * sizeof(double));
                leftMatrixSingleRowSlice[leftMatrixCol] = mul_slice_data->leftMatrix[leftMatrixRow + leftMatrixCol * matrixDimension];

//                cblas_dgemm(CblasColMajor,
//                            CblasNoTrans,
//                            CblasNoTrans,
//                            matrixDimension, //rows in A, C
//                            mul_slice_data->sliceWidth, //cols in B, C
//                            matrixDimension, //cols in A, rows in B
//                            1.0,
//                            leftMatrixSingleRowSlice,
//                            matrixDimension, //stride of A
//                            mul_slice_data->rightMatrixSlice,
//                            matrixDimension, //stride of B
//                            1.0,
//                            &(mul_slice_data->resultMatrix)[resultMatrixIndex],
//                            matrixDimension); //stride of C

                free(leftMatrixSingleRowSlice);

                resultMatrixIndex++;
            }
        }
    }

    pthread_exit(NULL);
}

void pThreadMultiply(int numThreads, int matrixDimension, double *leftMatrix, const double *rightMatrix,
                     double *pthreadResultMatrix) {
    pthread_t *working_thread;
    void *thread_status;
    mul_slice_data *thread_mul_slice_data;
    int sliceWidth = matrixDimension / numThreads;
    double rightMatrixSlice[numThreads][matrixDimension * sliceWidth];

    working_thread = malloc(numThreads * sizeof(pthread_t));
    thread_mul_slice_data = malloc(numThreads * sizeof(mul_slice_data));

    //create a vertical slice of rightMatrix per thread
    //TODO: support slice width > 1
    for (int thread = 0; thread < numThreads; thread++) {
        for (int col = 0; col < matrixDimension; col++) {
            for (int row = 0; row < matrixDimension; row++) {
                rightMatrixSlice[thread][col * thread + row] = rightMatrix[col * thread + row];
            }
        }
    }

    //create threads
    for (int thread = 0; thread < numThreads; thread++) {
        //construct slice data
        thread_mul_slice_data[thread].leftMatrix = leftMatrix;
        thread_mul_slice_data[thread].rightMatrixSlice = rightMatrixSlice[thread];
        thread_mul_slice_data[thread].resultMatrix = pthreadResultMatrix;
        thread_mul_slice_data[thread].matrixDimension = matrixDimension;
        thread_mul_slice_data[thread].sliceWidth = sliceWidth;

        //create thread to calculate slice
        pthread_create(&working_thread[thread], NULL, pThreadMultiplySlice,
                       (void *) &thread_mul_slice_data[thread]);
    }

    //join threads
    for (int thread = 0; thread < numThreads; thread++) {
        pthread_join(working_thread[thread], &thread_status);
    }

    free(working_thread);
    free(thread_mul_slice_data);
}

void assertMatricesAreEquivalent(int matrixDimension, const double *serialMulResultMatrix,
                                 const double *pthreadResultMatrix) {
    bool matrixValuesAreEqual = true;
    for (int col = 0; col < matrixDimension; col++) {
        for (int row = 0; row < matrixDimension; row++) {
            double serialMulElement = serialMulResultMatrix[col * matrixDimension + row];
            double pthreadMulElement = pthreadResultMatrix[col * matrixDimension + row];
            if (serialMulElement != pthreadMulElement) {
                // print all non-matching values before exiting
                printf("Matrix elements at column %d, row %d are different. \n Serial mul matrix value: %f \n Pthread mul matrix value: %f \n\n",
                       col, row, serialMulElement, pthreadMulElement);
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
    double *serialMulResultMatrix;
    double *pthreadMulResultMatrix;
    struct timeval tv1, tv2;
    struct timezone tz;
    int matrixDimension = 2;
    int numThreads = 2;

    unsigned long matrixMemorySize = matrixDimension * matrixDimension * sizeof(double);

    serialMulResultMatrix = malloc(matrixMemorySize);
    pthreadMulResultMatrix = malloc(matrixMemorySize);

    if (!serialMulResultMatrix || !pthreadMulResultMatrix) {
        printf("Insufficient memory for matrices of dimension %d.\n", matrixDimension);
        exit(-1);
    }

    // Serial matrix multiplication
    gettimeofday(&tv1, &tz);
    serialMultiply(matrixDimension, leftMatrix, rightMatrix, serialMulResultMatrix);
    gettimeofday(&tv2, &tz);
    double serialMulTimeElapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;

    // pthread matrix multiplication
    gettimeofday(&tv1, &tz);
//    pThreadMultiply(numThreads, matrixDimension, leftMatrix, rightMatrix, pthreadMulResultMatrix);
    gettimeofday(&tv2, &tz);
    double pthreadMulTimeElapsed =
            (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;

    assertMatricesAreEquivalent(matrixDimension, serialMulResultMatrix, pthreadMulResultMatrix);

    free(serialMulResultMatrix);
    free(pthreadMulResultMatrix);

    printf("Average times taken for matrix multiplication on array with %dx%d dimensions: \n Serial: %f \n Pthread: %f\n\n",
           matrixDimension, matrixDimension, serialMulTimeElapsed, pthreadMulTimeElapsed);

    return 0;
}