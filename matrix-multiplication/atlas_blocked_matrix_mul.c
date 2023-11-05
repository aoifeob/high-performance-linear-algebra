#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <cblas.h>

void initMatrix(int matrixDimension, double matrix[]) {
    for (int i = 0; i < matrixDimension; i++) {
        for (int j = 0; j < matrixDimension; j++) {
            double value = rand();
            matrix[i * matrixDimension + j] = value;
        }
    }
}

int main(void) {
//    double *leftMatrix;
//    double *rightMatrix;
    double leftMatrix[4] = {1, 1, 0, 0};
    double rightMatrix[4] = {1, 0, 1, 0};
    double *resultMatrix;

    struct timeval tv1, tv2;
    struct timezone tz;
    double blockMultiplicationTimeElapsed = 0;

    int resultMatrixIndex = 0;
    int matrixDimension = 2;
    int iterations = 10;

    unsigned long matrixMemorySize = matrixDimension * matrixDimension * sizeof(double);

//    leftMatrix = malloc(matrixMemorySize);
//    rightMatrix = malloc(matrixMemorySize);
    resultMatrix = malloc(matrixMemorySize);

//    initMatrix(matrixDimension, leftMatrix);
//    initMatrix(matrixDimension, rightMatrix);

//    if (!leftMatrix || !rightMatrix || !resultMatrix) {
//        printf("Insufficient memory for matrices of dimension %d.\n", matrixDimension);
//        exit(-1);
//    }

    if (!resultMatrix) {
        printf("Insufficient memory for matrices of dimension %d.\n", matrixDimension);
        exit(-1);
    }

    //loop iterating through each block of the right matrix
    for (int current_block_num_of_right_matrix = 0;
         current_block_num_of_right_matrix < matrixDimension;
         current_block_num_of_right_matrix++) {

        //create var to hold a block of the right matrix
        double *block_of_right_matrix;
        block_of_right_matrix = malloc(matrixDimension * sizeof(double));

        //build the block of the right matrix
        for (int current_col_in_right_matrix = 0;
             current_col_in_right_matrix < matrixDimension; current_col_in_right_matrix++) {
            block_of_right_matrix[current_col_in_right_matrix] = rightMatrix[current_block_num_of_right_matrix +
                                                                             current_col_in_right_matrix *
                                                                             matrixDimension];
        }

        //loop iterating through each block of the left matrix
        for (int current_block_num_of_left_matrix = 0;
             current_block_num_of_left_matrix < matrixDimension;
             current_block_num_of_left_matrix++) {

            //create var to hold a block of the left matrix
            double *block_of_left_matrix;
            block_of_left_matrix = malloc(matrixDimension * sizeof(double));

            //build the block of the left matrix
            for (int current_row_in_left_matrix = 0;
                 current_row_in_left_matrix < matrixDimension;
                 current_row_in_left_matrix++) {

                block_of_left_matrix[current_row_in_left_matrix] = leftMatrix[current_row_in_left_matrix +
                                                                              current_block_num_of_left_matrix *
                                                                              matrixDimension];

            }

            gettimeofday(&tv1, &tz);

            ATL_dgemm(CblasNoTrans,
                      CblasNoTrans,
                      1, //rows in A, C
                      1, //cols in B, C
                      matrixDimension, //cols in A, rows in B
                      1.0,
                      block_of_left_matrix,
                      1, //stride of A
                      block_of_right_matrix,
                      1, //stride of B
                      1.0,
                      &resultMatrix[resultMatrixIndex],
                      1); //stride of C

            printf("Result matrix has value %f at index %d\n", resultMatrix[resultMatrixIndex], resultMatrixIndex);

            gettimeofday(&tv2, &tz);
            double timeElapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;
            blockMultiplicationTimeElapsed += timeElapsed;

            //increment var tracking the index of the result matrix to be calculated
            resultMatrixIndex++;

            free(block_of_left_matrix);
        }

        free(block_of_right_matrix);
    }
    printf("\nTotal time taken for simple matrix multiplication on array with %dx%d dimensions: %f.\n\n",
           matrixDimension, matrixDimension, blockMultiplicationTimeElapsed);

//    free(leftMatrix);
//    free(rightMatrix);
    free(resultMatrix);

    return 0;
}