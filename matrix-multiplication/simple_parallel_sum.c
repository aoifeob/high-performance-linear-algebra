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

        int threadNum = omp_get_thread_num();

        printf("Thread %d beginning execution\n", threadNum);

#pragma omp parallel for reduction(+: sumTotal)
        int innerThreadNum = omp_get_thread_num();

        for (int i = 0; i < numProcesses; i++) {
                printf("Thread %d incrementing total by %d \n", innerThreadNum, i);
                sumTotal += i;
            }

        printf("Total sum is %d \n", sumTotal);

#pragma omp critical
        myInt++;
        printf("Thread %d incremented myInt to %d\n", threadNum, myInt);
    }

    return 0;
}