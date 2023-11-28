#include <stdio.h>
#include <omp.h>

int main(void) {
    int numProcesses = omp_get_num_procs();
    int sumTotal = 0;

    printf("Number of processes is %d \n", numProcesses);

    omp_set_num_threads(numProcesses);

    #pragma omp parallel for reduction(+: sumTotal)
    for (int i = 0; i < numProcesses; i++) {
        sumTotal += i;
    }

    printf("Total sum is %d \n", sumTotal);

    return 0;
}