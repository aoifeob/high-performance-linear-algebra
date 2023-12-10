//#include <stdio.h>
//#include <stdlib.h>
//#include <mpi.h>
//
//double dot_product(double *x, double *y, int m, MPI_Comm comm) {
//    double local_dot_prod;
//    double global_dot_prod; /* result (at process zero) */
//    int i;
//    local_dot_prod = 0.;
//    for(i=0; i<m; i++)
//        local_dot_prod += x[i]*y[i];
//    MPI_Reduce(&local_dot_prod,
//               &global_dot_prod,
//               1,
//               MPI_DOUBLE,
//               MPI_SUM,
//               0,
//               comm);
//    return global_dot_prod;
//}
//
//int main(void) {
//    int numProcesses = 8;
//    int sumTotal = 0;
//    int my_global_rank, my_rank_in_a;
//    int global_size;
//
//    MPI_Init(&argc, &argv);
//    MPI_Comm_size(MPI_COMM_WORLD, &global_size);
//    MPI_Comm_rank(MPI_COMM_WORLD, &my_global_rank);
//
//    dot_product(numProcesses, &sumTotal);
//
//    printf("Total sum is %d \n", sumTotal);
//
//    if(MPI_COMM_WORLD != MPI_COMM_NULL)
//        MPI_Comm_free(MPI_COMM_WORLD);
//
//    MPI_Finalize();
//
//    return 0;
//}