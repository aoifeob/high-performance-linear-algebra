#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <mpi.h>

int matrixDimension, p;

int main(int argc, char **argv) {
    double leftMatrix[4] = {1, 2, 3, 4};
    double rightMatrix[4] = {5, 6, 7, 8};
    int thisProcRank, squareSize, neighboursToReceiveFrom, numElementsToReceive;
    double *leftSubMatrix, *rightSubMatrix, *resultSubMatrix, *leftRow, *rightCol, *parallelResultMatrix, parallelMulStartTime, parallelMulTimeElapsed;

    if (argc < 2) {
        printf("Invalid number of arguments supplied. \n");
        exit(-1);
    }

    matrixDimension = atoi(argv[1]);

//    p = 16;
//    matrixDimension = 4;
//    thisProcRank = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    int rootP = (int) sqrt(p);

    if (rootP - sqrt(p) != 0.0) {
        printf("Invalid number of processors: %d. P must be a square number. \n", p);
        exit(-1);
    }

    if (rootP % matrixDimension != 0) {
        printf("Invalid matrix size: %d for number of processors: %d. Sqrt p must be divisible by n. \n",
               matrixDimension, p);
        exit(-1);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &thisProcRank);

    squareSize = matrixDimension / rootP;
    neighboursToReceiveFrom = rootP - 1;
    numElementsToReceive = (matrixDimension * matrixDimension) / p;

    leftSubMatrix = malloc(squareSize * matrixDimension * sizeof(double));
    rightSubMatrix = malloc(squareSize * matrixDimension * sizeof(double));
    resultSubMatrix = malloc(squareSize * matrixDimension * sizeof(double));
    leftRow = malloc(matrixDimension * squareSize * sizeof(double));
    rightCol = malloc(matrixDimension * squareSize * sizeof(double));

    if (!leftSubMatrix || !rightSubMatrix || !resultSubMatrix || !leftRow || !rightCol) {
        printf("Insufficient memory for matrices of size: %d", matrixDimension);
        exit(-1);
    }

    //determine where the squares the processor is responsible for are located in the whole matrix
    int thisProcCol = (int) floor((thisProcRank / squareSize));
    int thisProcRow = thisProcRank % squareSize;
    int fullMatrixStartingIndex = thisProcRow * squareSize + (matrixDimension * squareSize * thisProcCol);

    // populate values of left/right subMatrices for processor
    int resultSubMatrixIndex = 0;
    for (int col = fullMatrixStartingIndex; col < fullMatrixStartingIndex + squareSize; col++) {
        for (int row = fullMatrixStartingIndex; row < fullMatrixStartingIndex + squareSize; row++) {
            leftSubMatrix[resultSubMatrixIndex] = leftMatrix[fullMatrixStartingIndex + (matrixDimension * col) + row];
            rightSubMatrix[resultSubMatrixIndex] = rightMatrix[fullMatrixStartingIndex + (matrixDimension * col) + row];
            resultSubMatrixIndex++;
        }
    }

    // synchronise to ensure all processors know their squares
    MPI_Barrier(MPI_COMM_WORLD);

    // create communicators for communicating rows and cols
    MPI_Comm colComm, rowComm;
    MPI_Comm_split(MPI_COMM_WORLD, thisProcCol, thisProcRank, &colComm);
    MPI_Comm_split(MPI_COMM_WORLD, thisProcRow, thisProcRank, &rowComm);

    // begin tracking time
    if (thisProcRank == 0) {
        parallelMulStartTime = MPI_Wtime();
    }

    // all processors send required square data to other processors
    MPI_Allgather(leftSubMatrix,
                  squareSize * squareSize,
                  MPI_DOUBLE,
                  leftRow,
                  squareSize * squareSize,
                  MPI_DOUBLE,
                  rowComm);
    MPI_Allgather(rightSubMatrix,
                  squareSize * squareSize,
                  MPI_DOUBLE,
                  rightCol,
                  squareSize * squareSize,
                  MPI_DOUBLE,
                  colComm);

    // calculate result sub matrix values
    for (int col = 0; col < squareSize; col++) {
        for (int row = 0; row < squareSize; row++) {
            double element = 0;
            for (int k = 0; k < squareSize; k++) {
                element += leftRow[row + k * squareSize] * rightCol[k + squareSize * col];
            }
            resultSubMatrix[col * squareSize + row] = element;
        }
    }

    free(leftRow);
    free(rightCol);

    //synchronise to ensure all processors have finished calculating their result submatrix
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

    printf("Time taken for parallel matrix multiplication on array with %dx%d dimensions: %f\n\n",
           matrixDimension, matrixDimension, parallelMulTimeElapsed);

    MPI_Finalize();

    if (thisProcRank == 0) {
        free(parallelResultMatrix);
    }
    free(leftSubMatrix);
    free(rightSubMatrix);
    free(resultSubMatrix);
}
