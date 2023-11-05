#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
//#include <cblas.h>

void initMatrix(int matrixDimension, double matrix[]) {
    for (int i = 0; i < matrixDimension; i++) {
        for (int j = 0; j < matrixDimension; j++) {
            double value = rand();
            matrix[i * matrixDimension + j] = value;
        }
    }
}

int main(void) {
    double *leftMatrix;
    double *rightMatrix;
    double *resultMatrix;

    struct timeval tv1, tv2;
    struct timezone tz;
    double blockMultiplicationTimeElapsed = 0;

    int matrixDimension = 1024; //default value, can be overwritten by user input
    int blockSize = 2; //default value, can be overwritten by user input

//    printf("Enter matrix dimension n : \n\n");
//    scanf("%d", &matrixDimension);
//
//    printf("Enter matrix block width b : \n\n");
//    scanf("%d", &blockSize);

    unsigned long matrixMemorySize = matrixDimension * matrixDimension * sizeof(double);

    leftMatrix = malloc(matrixMemorySize);
    rightMatrix = malloc(matrixMemorySize);
    resultMatrix = malloc(matrixMemorySize);

    initMatrix(matrixDimension, leftMatrix);
    initMatrix(matrixDimension, rightMatrix);

    if (!leftMatrix || !rightMatrix || !resultMatrix) {
        printf("Insufficient memory for matrices of dimension %d.\n", matrixDimension);
        exit(-1);
    }

    int resultMatrixIndex = 0;
    //loop iterating through each block of the right matrix
    for (int current_block_num_of_right_matrix = 0;
         current_block_num_of_right_matrix < matrixDimension / blockSize;
         current_block_num_of_right_matrix++) {

        //create var to hold a block of the right matrix to be stored in cache
        double *cached_block_of_right_matrix;
        cached_block_of_right_matrix = malloc(matrixDimension * blockSize * sizeof(double));

        //build the block of the right matrix to be stored in cache
        for (int current_col_in_right_block = 0;
             current_col_in_right_block < blockSize;
             current_col_in_right_block++) {
            for (int current_col_in_right_matrix = 0;
                 current_col_in_right_matrix < matrixDimension * blockSize;
                 current_col_in_right_matrix++) {
                cached_block_of_right_matrix[matrixDimension * current_col_in_right_block +
                                             current_col_in_right_matrix] = rightMatrix[
                        current_block_num_of_right_matrix + current_col_in_right_block +
                        current_col_in_right_matrix *
                        matrixDimension];
            }
        }

        double *block_of_right_matrix_to_multiply;
        block_of_right_matrix_to_multiply = malloc(matrixDimension * sizeof(double));

        if (blockSize > 1) {
            //build the block of the right matrix to use for multiplication
            for (int current_col_in_right_block = 0;
                 current_col_in_right_block < matrixDimension;
                 current_col_in_right_block++) {
                block_of_right_matrix_to_multiply[current_col_in_right_block] = cached_block_of_right_matrix[
                        current_block_num_of_right_matrix +
                        current_col_in_right_block * matrixDimension];
            }
        }

        //loop iterating through each block of the left matrix
        for (int current_block_num_of_left_matrix = 0;
             current_block_num_of_left_matrix < matrixDimension / blockSize;
             current_block_num_of_left_matrix++) {

            //create var to hold a block of the left matrix to be stored in cache
            double *cached_block_of_left_matrix;
            cached_block_of_left_matrix = malloc(matrixDimension * blockSize * sizeof(double));

            //build the block of the left matrix to be stored in cache
            for (int current_row_in_left_block = 0;
                 current_row_in_left_block < blockSize;
                 current_row_in_left_block++) {
                for (int current_row_in_left_matrix = 0;
                     current_row_in_left_matrix < matrixDimension * blockSize;
                     current_row_in_left_matrix++) {

                    cached_block_of_left_matrix[current_row_in_left_block + matrixDimension +
                                                current_row_in_left_matrix] = leftMatrix[current_row_in_left_matrix +
                                                                                         current_block_num_of_left_matrix *
                                                                                         matrixDimension];

                }
            }

            double *block_of_left_matrix_to_multiply;
            block_of_left_matrix_to_multiply = malloc(matrixDimension * sizeof(double));

            //if the cached block contains more than one row/col, isolate the row and col to be used for this multiplication
            if (blockSize > 1) {
                //build the block of the left matrix to use for multiplication
                for (int current_row_in_left_block = 0;
                     current_row_in_left_block < matrixDimension;
                     current_row_in_left_block++) {
                    block_of_left_matrix_to_multiply[current_row_in_left_block] = cached_block_of_left_matrix[
                            current_block_num_of_left_matrix * matrixDimension +
                            current_row_in_left_block];
                }

                gettimeofday(&tv1, &tz);

//                ATL_dgemm(CblasNoTrans,
//                          CblasNoTrans,
//                          1, //rows in A, C
//                          1, //cols in B, C
//                          matrixDimension, //cols in A, rows in B
//                          1.0,
//                          block_of_left_matrix_to_multiply,
//                          1, //stride of A
//                          block_of_right_matrix_to_multiply,
//                          1, //stride of B
//                          1.0,
//                          &resultMatrix[resultMatrixIndex],
//                          1); //stride of C

                gettimeofday(&tv2, &tz);
            } else {
                gettimeofday(&tv1, &tz);

//                ATL_dgemm(CblasNoTrans,
//                          CblasNoTrans,
//                          1, //rows in A, C
//                          1, //cols in B, C
//                          matrixDimension, //cols in A, rows in B
//                          1.0,
//                          cached_block_of_left_matrix,
//                          1, //stride of A
//                          cached_block_of_right_matrix,
//                          1, //stride of B
//                          1.0,
//                          &resultMatrix[resultMatrixIndex],
//                          1); //stride of C

                gettimeofday(&tv2, &tz);
            }

            double timeElapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;
            blockMultiplicationTimeElapsed += timeElapsed;

            //increment var tracking the index of the result matrix to be calculated
            resultMatrixIndex++;

            free(block_of_left_matrix_to_multiply);
            free(cached_block_of_left_matrix);
        }

        free(block_of_right_matrix_to_multiply);
        free(cached_block_of_right_matrix);
    }

    free(leftMatrix);
    free(rightMatrix);
    free(resultMatrix);

    printf("Time taken for simple matrix multiplication on array with %dx%d dimensions: %f.\n\n",
           matrixDimension, matrixDimension, blockMultiplicationTimeElapsed);

    return 0;
}