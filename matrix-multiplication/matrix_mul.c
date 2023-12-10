#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <sys/time.h>
#include <mpi.h>

int matrixDimension, p;

void initMatrix(int dim, double matrix[]) {
    for (int col = 0; col < dim; col++) {
        for (int row = 0; row < dim; row++) {
            double value = rand();
            matrix[col * dim + row] = value;
        }
    }
}

void serialMultiply(int dim, const double *leftMatrix, const double *rightMatrix,
                    double *serialMulResultMatrix) {
    for (int col = 0; col < dim; col++) {
        for (int row = 0; row < dim; row++) {
            double element = 0;
            for (int k = 0; k < dim; k++) {
                element += leftMatrix[row + k * dim] * rightMatrix[k + dim * col];
            }
            serialMulResultMatrix[col * dim + row] = element;
        }
    }
}

void assertMatricesAreEquivalent(int dim, const double *serialMulResultMatrix,
                                 const double *parallelResultMatrix) {
    bool matrixValuesAreEqual = true;
    for (int col = 0; col < dim; col++) {
        for (int row = 0; row < dim; row++) {
            double serialMulElement = serialMulResultMatrix[col * dim + row];
            double parallelMulElement = parallelResultMatrix[col * dim + row];
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

int main(int argc, char **argv) {
    int thisProcRank, squareSize;
    double *leftMatrix, *rightMatrix, *serialResultMatrix;
    double *leftSubMatrix, *rightSubMatrix, *resultSubMatrix, *leftRows, *rightCols, *parallelResultMatrix, parallelMulStartTime, parallelMulTimeElapsed;
    struct timeval tv1, tv2;
    struct timezone tz;
    int shouldRunSerialProgram = 0;

    if (argc < 2) {
        printf("Invalid number of arguments supplied. \n");
        exit(-1);
    }

    matrixDimension = atoi(argv[1]);

    printf("This program supports serial and parallel one-norm computation. To disable serial computation, enter 1. Otherwise, enter 0.\n\n");
    scanf("%d", &shouldRunSerialProgram);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    int rootP = (int) sqrt(p);

    if (rootP - sqrt(p) != 0.0) {
        printf("Invalid number of processes: %d. P must be a square number. \n", p);
        exit(-1);
    }

    if (matrixDimension % rootP != 0) {
        printf("Invalid matrix size: %d for number of processes: %d. Sqrt p (%d) must be divisible by n. \n",
               matrixDimension, p, rootP);
        exit(-1);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &thisProcRank);

    squareSize = matrixDimension / rootP;

    leftMatrix = malloc(matrixDimension * matrixDimension * sizeof(double));
    rightMatrix = malloc(matrixDimension * matrixDimension * sizeof(double));
    leftSubMatrix = malloc(squareSize * squareSize * sizeof(double));
    rightSubMatrix = malloc(squareSize * squareSize * sizeof(double));
    resultSubMatrix = malloc(squareSize * squareSize * sizeof(double));
    leftRows = malloc(matrixDimension * squareSize * sizeof(double));
    rightCols = malloc(matrixDimension * squareSize * sizeof(double));

    if (!leftMatrix || !rightMatrix || !leftSubMatrix || !rightSubMatrix || !resultSubMatrix || !leftRows || !rightCols) {
        printf("Insufficient memory for matrices of size: %d", matrixDimension);
        exit(-1);
    }

    initMatrix(matrixDimension, leftMatrix);
    initMatrix(matrixDimension, rightMatrix);

    //determine where the squares the process is responsible for are located in the whole matrix
    int thisProcCol = (int) floor(((double)thisProcRank / rootP));
    int thisProcRow = thisProcRank % rootP;
    int fullMatrixStartingIndex = thisProcRow * squareSize + (matrixDimension * squareSize * thisProcCol);

    // populate values of left/right subMatrices for process
    for (int col = 0; col < squareSize; col++) {
        for (int row = 0; row < squareSize; row++) {
            leftSubMatrix[col * squareSize + row] = leftMatrix[fullMatrixStartingIndex + (matrixDimension * col) + row];
            rightSubMatrix[col * squareSize + row] = rightMatrix[fullMatrixStartingIndex + (matrixDimension * col) + row];
        }
    }

    // synchronise to ensure all processes know their squares
    MPI_Barrier(MPI_COMM_WORLD);

    // create communicators for communicating rows and cols
    MPI_Comm colComm, rowComm;
    MPI_Comm_split(MPI_COMM_WORLD, thisProcCol, thisProcRank, &colComm);
    MPI_Comm_split(MPI_COMM_WORLD, thisProcRow, thisProcRank, &rowComm);

    // begin tracking time
    if (thisProcRank == 0) {
        parallelMulStartTime = MPI_Wtime();
    }

    // all processes send required square data to other processes
    MPI_Allgather(rightSubMatrix,
                  squareSize * squareSize,
                  MPI_DOUBLE,
                  rightCols,
                  squareSize * squareSize,
                  MPI_DOUBLE,
                  colComm);

    MPI_Allgather(leftSubMatrix,
                  squareSize * squareSize,
                  MPI_DOUBLE,
                  leftRows,
                  squareSize * squareSize,
                  MPI_DOUBLE,
                  rowComm);

    // calculate result sub matrix values
    for (int col = 0; col < squareSize; col++) {
        for (int row = 0; row < squareSize; row++) {
            double sum = 0;
            for (int k = 0; k < matrixDimension; k++) {
                sum += leftRows[row + k * squareSize] * rightCols[k + matrixDimension * col];
            }
            resultSubMatrix[col * squareSize + row] = sum;
        }
    }

    free(leftRows);
    free(rightCols);

    //synchronise to ensure all processes have finished calculating their result submatrix
    MPI_Barrier(MPI_COMM_WORLD);

    if (thisProcRank == 0) {
        parallelMulTimeElapsed = MPI_Wtime() - parallelMulStartTime;
        parallelResultMatrix = malloc(matrixDimension * matrixDimension * sizeof(double));
    }

    //gather full result matrix
    MPI_Gather(resultSubMatrix,
               squareSize * squareSize,
               MPI_DOUBLE,
               parallelResultMatrix,
               squareSize * squareSize,
               MPI_DOUBLE,
               0,
               MPI_COMM_WORLD);

    //synchronise to ensure all processes have finished sending their result sub matrices
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();

    if (thisProcRank == 0) {
        if (shouldRunSerialProgram){
            if (!serialResultMatrix){
                printf("Insufficient memory for matrices of size: %d", matrixDimension);
                exit(-1);
            }

            // Serial matrix multiplication
            gettimeofday(&tv1, &tz);
            serialMultiply(matrixDimension, leftMatrix, rightMatrix, serialResultMatrix);
            gettimeofday(&tv2, &tz);
            double serialMulTimeElapsed = (double) (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) * 1.e-6;

            assertMatricesAreEquivalent(matrixDimension, serialResultMatrix, parallelResultMatrix);

            printf("Times taken for matrix multiplication on array with %dx%d dimensions: \n Serial: %f \n Parallel: %f\n\n",
                   matrixDimension, matrixDimension, serialMulTimeElapsed, parallelMulTimeElapsed);

            free(serialResultMatrix);
        }
        else {
            printf("Time taken for parallel matrix multiplication on array with %dx%d dimensions: %f\n\n",
                   matrixDimension, matrixDimension, parallelMulTimeElapsed);
        }
        free(parallelResultMatrix);
    }

    free(leftMatrix);
    free(rightMatrix);
    free(leftSubMatrix);
    free(rightSubMatrix);
    free(resultSubMatrix);

    return 0;
}
