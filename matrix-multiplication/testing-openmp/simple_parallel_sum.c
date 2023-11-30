#include <stdio.h>
#include <omp.h>

int main(void) {
    int numProcesses = omp_get_num_procs();
    int sumTotal = 0;

    omp_set_num_threads(numProcesses);
    printf("Number of threads is %d \n", omp_get_num_threads());

#pragma omp parallel shared(myInt)
    {

        printf("Thread %d beginning execution\n", omp_get_thread_num());

    #pragma omp parallel for reduction(+: sumTotal)
        for (int i = 0; i < numProcesses; i++) {
            sumTotal += i;
        }

        printf("Total sum is %d \n", sumTotal);
    }

    return 0;
}

