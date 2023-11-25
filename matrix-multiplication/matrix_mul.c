#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <math.h>
#include <pthread.h>

#define MAX_THREADS 124

typedef struct {
    double *leftMatrix;
    double *rightMatrixSlice;
    double *resultMatrix;
    int matrixDimension;
    int sliceWidth;
    int resultMatrixSliceStartingIndex;
    int threadNumber;
} mul_slice_data;

typedef struct {
    double *resultMatrixSlice;
    double *oneNorm;
    int matrixDimension;
    int sliceWidth;
    pthread_mutex_t *mutex;
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

void *pThreadMultiplySlice(void *arg) {
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

    pthread_exit(NULL);
}

void pThreadMultiply(int numThreads, int matrixDimension, double *leftMatrix, double *rightMatrix,
                     double *pthreadResultMatrix) {
    pthread_t *working_thread;
    void *thread_status;
    mul_slice_data *thread_mul_slice_data;
    int sliceWidth = matrixDimension / numThreads;
    int elementsInSlice = matrixDimension * sliceWidth;

    working_thread = malloc(numThreads * sizeof(pthread_t));
    thread_mul_slice_data = malloc(numThreads * sizeof(mul_slice_data));

    //create threads
    for (int thread = 0; thread < numThreads; thread++) {
        //construct slice data
        thread_mul_slice_data[thread].leftMatrix = leftMatrix;
        thread_mul_slice_data[thread].rightMatrixSlice = rightMatrix + thread * elementsInSlice;
        thread_mul_slice_data[thread].resultMatrix = pthreadResultMatrix;
        thread_mul_slice_data[thread].matrixDimension = matrixDimension;
        thread_mul_slice_data[thread].sliceWidth = isLastThread(thread, numThreads) ? matrixDimension - (sliceWidth * (numThreads-1)) : sliceWidth;
        thread_mul_slice_data[thread].resultMatrixSliceStartingIndex = thread * matrixDimension * sliceWidth;
        thread_mul_slice_data[thread].threadNumber = thread;

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

void *pThreadCalculateSliceNorm(void *arg) {
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

    //lock mutex before reading
    pthread_mutex_lock(norm_slice_data->mutex);

    //if norm of the current column is greater than the current max column norm, update it to the current value
    if (sliceNorm > *(norm_slice_data->oneNorm)) {
        *(norm_slice_data->oneNorm) = sliceNorm;
    }

    //unlock mutex after read/write
    pthread_mutex_unlock(norm_slice_data->mutex);

    pthread_exit(NULL);
}

void calculatePthreadNorm(int numThreads, int matrixDimension, double *pthreadMulResultMatrix, double *oneNorm) {
    pthread_t *working_thread;
    pthread_mutex_t *mutex_one_norm;
    void *thread_status;
    norm_slice_data *thread_norm_slice_data;
    int sliceWidth = matrixDimension / numThreads;
    int elementsInSlice = matrixDimension * sliceWidth;

    working_thread = malloc(numThreads * sizeof(pthread_t));
    thread_norm_slice_data = malloc(numThreads * sizeof(mul_slice_data));
    mutex_one_norm = malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(mutex_one_norm, NULL);

    for (int thread = 0; thread < numThreads; thread++) {
        //construct slice data
        thread_norm_slice_data[thread].resultMatrixSlice = pthreadMulResultMatrix + thread * elementsInSlice;
        thread_norm_slice_data[thread].oneNorm = oneNorm;
        thread_norm_slice_data[thread].matrixDimension = matrixDimension;
        thread_norm_slice_data[thread].sliceWidth = isLastThread(thread, numThreads) ? matrixDimension - (sliceWidth * (numThreads-1)) : sliceWidth;
        thread_norm_slice_data[thread].mutex = mutex_one_norm;

        //create thread to calculate norm for cols in slice
        pthread_create(&working_thread[thread], NULL, pThreadCalculateSliceNorm,
                       (void *) &thread_norm_slice_data[thread]);
    }

    //join threads
    for (int thread = 0; thread < numThreads; thread++) {
        pthread_join(working_thread[thread], &thread_status);
    }

    free(working_thread);
    free(thread_norm_slice_data);
    pthread_mutex_destroy(mutex_one_norm);
    free(mutex_one_norm);
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

void assertNormsAreEquivalent(double serialNorm, double pThreadNorm) {
    if (serialNorm != pThreadNorm) {
        printf("Matrix norms are different. Serial norm: %f \n Pthread norm: %f \n\n", serialNorm, pThreadNorm);
        exit(-1);
    }
}

int main(void) {
    double *leftMatrix, *rightMatrix;
    double *serialMulResultMatrix;
    double *pthreadMulResultMatrix;
    double serialNorm, pThreadNorm;
    struct timeval tv1, tv2;
    struct timezone tz;
    int matrixDimension = 2048; //default value, can be overwritten by user input
    int numThreads = 8;

    printf("Enter matrix dimension n : \n\n");
    scanf("%d", &matrixDimension);

    printf("Enter number of working threads p: \n\n");
    if (scanf("%d", &numThreads) < 1 || numThreads > MAX_THREADS) {
        printf("Invalid number of threads %d specified", numThreads);
        exit(-1);
    }

    if (numThreads > matrixDimension) {
        printf("Number of threads p: %d should be smaller than matrix dimension n: %d\n\n",
               matrixDimension, numThreads);
        exit(-1);
    }

    if (0 != matrixDimension % numThreads) {
        printf("Matrix with dimension n: %d and number of threads p: %d will be partitioned into uneven slices\n\n",
               matrixDimension, numThreads);
    }

    unsigned long matrixMemorySize = matrixDimension * matrixDimension * sizeof(double);

    leftMatrix = (double *)malloc(matrixMemorySize);
    rightMatrix = (double *)malloc(matrixMemorySize);
    serialMulResultMatrix = (double *)malloc(matrixMemorySize);
    pthreadMulResultMatrix = (double *)malloc(matrixMemorySize);

    if (!leftMatrix || !rightMatrix || !serialMulResultMatrix || !pthreadMulResultMatrix) {
        printf("Insufficient memory for matrices of dimension %d.\n", matrixDimension);
        exit(-1);
    }

    initMatrix(matrixDimension, leftMatrix);
    initMatrix(matrixDimension, rightMatrix);

    // Serial matrix multiplication
    gettimeofday(&tv1, &tz);
    serialMultiply(matrixDimension, leftMatrix, rightMatrix, serialMulResultMatrix);
    serialNorm = calculateSerialNorm(matrixDimension, serialMulResultMatrix);
    gettimeofday(&tv2, &tz);
    double serialMulTimeElapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;


    // pthread matrix multiplication
    gettimeofday(&tv1, &tz);
    pThreadMultiply(numThreads, matrixDimension, leftMatrix, rightMatrix, pthreadMulResultMatrix);
    calculatePthreadNorm(numThreads, matrixDimension, pthreadMulResultMatrix, &pThreadNorm);
    gettimeofday(&tv2, &tz);
    double pthreadMulTimeElapsed =
            (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;

    assertMatricesAreEquivalent(matrixDimension, serialMulResultMatrix, pthreadMulResultMatrix);
    assertNormsAreEquivalent(serialNorm, pThreadNorm);

    free(leftMatrix);
    free(rightMatrix);
    free(serialMulResultMatrix);
    free(pthreadMulResultMatrix);

    printf("Times taken for matrix multiplication on array with %dx%d dimensions: \n Serial: %f \n Pthread: %f\n\n",
           matrixDimension, matrixDimension, serialMulTimeElapsed, pthreadMulTimeElapsed);

    return 0;
}