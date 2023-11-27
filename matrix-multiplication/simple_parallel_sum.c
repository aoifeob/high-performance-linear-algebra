#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int valueToSum;
    int *totalSum;
} sum_data;

void *sum(void *arg) {
    sum_data *sum_data = arg;

    *(sum_data->totalSum) += sum_data->valueToSum;
}

void parallelSum(int numProcesses, int *totalSum) {
    sum_data *process_sum_data;

    process_sum_data = malloc(numProcesses * sizeof(sum_data));

    //create processes
    for (int processes = 0; processes < numProcesses; processes++) {
        process_sum_data[processes].valueToSum = processes;
        process_sum_data[processes].totalSum = totalSum;

        //calculate slice
    }

    free(process_sum_data);
}

int main(void) {
    int numProcesses = 8;
    int sumTotal = 0;

    parallelSum(numProcesses, &sumTotal);

    printf("Total sum is %d \n", sumTotal);

    return 0;
}