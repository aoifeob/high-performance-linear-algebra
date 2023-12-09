#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <mpi.h>

int matrixDimension, p;

void printMatrix(int dim, double matrix[]) {
    for (int col = 0; col < dim; col++) {
        for (int row = 0; row < dim; row++) {
            printf("Matrix col: %d, row: %d has value %f\n", col, row, matrix[col * dim + row]);
        }
        printf("\n");
    }
    printf("\n");
}

void printSlice(int rows, int cols, int dim, double matrix[]) {
    for (int col = 0; col < cols; col++) {
        for (int row = 0; row < rows; row++) {
            printf("Slice col: %d, row: %d has value %f\n", col, row, matrix[col * dim + row]);
        }
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char **argv) {
    double leftMatrix[4] = {1, 2, 3, 4};
    double rightMatrix[4] = {5, 6, 7, 8};
    int thisProcRank, squareSize, neighboursToReceiveFrom, numElementsToReceive;
    double *leftSubMatrix, *rightSubMatrix, *resultSubMatrix, *leftRows, *rightCols, *parallelResultMatrix, parallelMulStartTime, parallelMulTimeElapsed;

    if (argc < 2) {
        printf("Invalid number of arguments supplied. \n");
        exit(-1);
    }

    matrixDimension = atoi(argv[1]);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    int rootP = (int) sqrt(p);

    if (rootP - sqrt(p) != 0.0) {
        printf("Invalid number of processes: %d. P must be a square number. \n", p);
        exit(-1);
    }

    if (rootP % matrixDimension != 0) {
        printf("Invalid matrix size: %d for number of processes: %d. Sqrt p must be divisible by n. \n",
               matrixDimension, p);
        exit(-1);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &thisProcRank);

    squareSize = matrixDimension / rootP;
    neighboursToReceiveFrom = rootP - 1;
    numElementsToReceive = (matrixDimension * matrixDimension) / p;

    leftSubMatrix = malloc(squareSize * squareSize * sizeof(double));
    rightSubMatrix = malloc(squareSize * squareSize * sizeof(double));
    resultSubMatrix = malloc(squareSize * squareSize * sizeof(double));
    leftRows = malloc(matrixDimension * squareSize * sizeof(double));
    rightCols = malloc(matrixDimension * squareSize * sizeof(double));

    if (!leftSubMatrix || !rightSubMatrix || !resultSubMatrix || !leftRows || !rightCols) {
        printf("Insufficient memory for matrices of size: %d", matrixDimension);
        exit(-1);
    }

    //determine where the squares the process is responsible for are located in the whole matrix
    int thisProcCol = (int) floor((thisProcRank / squareSize));
    int thisProcRow = thisProcRank % squareSize;
    int fullMatrixStartingIndex = thisProcRow * squareSize + (matrixDimension * squareSize * thisProcCol);

    // populate values of left/right subMatrices for process
    int resultSubMatrixIndex = 0;
    for (int col = fullMatrixStartingIndex; col < fullMatrixStartingIndex + squareSize; col++) {
        for (int row = fullMatrixStartingIndex; row < fullMatrixStartingIndex + squareSize; row++) {
            leftSubMatrix[resultSubMatrixIndex] = leftMatrix[fullMatrixStartingIndex + (matrixDimension * col) + row];
            rightSubMatrix[resultSubMatrixIndex] = rightMatrix[fullMatrixStartingIndex + (matrixDimension * col) + row];
            resultSubMatrixIndex++;
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

    printf("Process %d finished gathering cols\n", thisProcRank);
    printSlice(squareSize, matrixDimension, matrixDimension, rightCols);

    MPI_Allgather(leftSubMatrix,
                  squareSize * squareSize,
                  MPI_DOUBLE,
                  leftRows,
                  squareSize * squareSize,
                  MPI_DOUBLE,
                  rowComm);

    printf("Process %d finished gathering rows\n", thisProcRank);
    printSlice(matrixDimension, squareSize, matrixDimension, rightCols);

    // calculate result sub matrix values
    for (int col = 0; col < squareSize; col++) {
        for (int row = 0; row < squareSize; row++) {
            double element = 0;
            for (int k = 0; k < squareSize; k++) {
                element += leftRows[row + k * squareSize] * rightCols[k + squareSize * col];
            }
            resultSubMatrix[col * squareSize + row] = element;
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
        printf("Time taken for parallel matrix multiplication on array with %dx%d dimensions: %f\n\n",
               matrixDimension, matrixDimension, parallelMulTimeElapsed);
        printMatrix(matrixDimension, parallelResultMatrix);
        free(parallelResultMatrix);
    }
    free(leftSubMatrix);
    free(rightSubMatrix);
    free(resultSubMatrix);
}
