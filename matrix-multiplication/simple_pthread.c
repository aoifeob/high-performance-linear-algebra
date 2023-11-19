#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#define MAX_THREADS 124

typedef struct {
    int valueToSum;
    int *totalSum;
    pthread_mutex_t *mutex;
} sum_data;

void *sum(void *arg) {
    sum_data *sum_data = arg;

    //lock mutex before reading
    pthread_mutex_lock(sum_data->mutex);

    *(sum_data->totalSum) += sum_data->valueToSum;

    //unlock mutex after read/write
    pthread_mutex_unlock(sum_data->mutex);

    pthread_exit(NULL);
}

void pThreadSum(int numThreads, int *totalSum) {
    pthread_t *working_thread;
    void *thread_status;
    sum_data *thread_sum_data;
    pthread_mutex_t *mutex_sum;

    working_thread = malloc(numThreads * sizeof(pthread_t));
    thread_sum_data = malloc(numThreads * sizeof(sum_data));

    mutex_sum = malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(mutex_sum, NULL);

    //create threads
    for (int thread = 0; thread < numThreads; thread++) {
        thread_sum_data[thread].valueToSum = thread;
        thread_sum_data[thread].totalSum = totalSum;
        thread_sum_data[thread].mutex = mutex_sum;

        //create thread to calculate slice
        pthread_create(&working_thread[thread], NULL, sum,
                       (void *) &thread_sum_data[thread]);
    }

    //join threads
    for (int thread = 0; thread < numThreads; thread++) {
        pthread_join(working_thread[thread], &thread_status);
    }

    free(working_thread);
    free(thread_sum_data);
}

int main(void) {
    int numThreads = 8;
    int sumTotal = 0;

    pThreadSum(numThreads, &sumTotal);

    printf("Total sum is %d \n", sumTotal);

    return 0;
}