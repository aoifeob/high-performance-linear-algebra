#include <stdio.h>
#include <omp.h>

int main(void) {
    int numProcesses = omp_get_num_procs();
    int sumTotal = 0;

    printf("Number of processes is %d \n", numProcesses);

    omp_set_num_threads(numProcesses);

    int myInt = 0;

#pragma omp parallel shared(myInt)
    {
        printf("Thread %d beginning execution\n", omp_get_thread_num());

#pragma omp parallel for reduction(+: sumTotal)
        for (int i = 0; i < numProcesses; i++) {
            printf("Thread %d incrementing total by %d \n", omp_get_thread_num(), i);
            sumTotal += i;
        }

        printf("Total sum is %d \n", sumTotal);

#pragma omp critical
        myInt++;
        printf("Thread %d incremented myInt to %d\n", omp_get_thread_num(), myInt);
    }

    return 0;
}