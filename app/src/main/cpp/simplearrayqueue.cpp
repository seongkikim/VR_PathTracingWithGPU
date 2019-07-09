//
// Created by user on 2019-07-09.
//

#include <stdlib.h>
#include <string.h>
#include "include/native-lib.h"
#include "include/CLBVH.h"

extern short width;
extern short height;

typedef struct {
    int genthread;
    Vec *colors;
    unsigned int *samples;
    unsigned int *samplesDiff;
    ToDiffInfo *tdis;
} QMEMBER;

#define MAX 10

QMEMBER queue[MAX];
int front, rear;

int queue_length() {
    int arear = rear;

    while(front > arear)
        arear += MAX;

    return arear - front;
}

void init_queue() {
    front = rear = 0;

    for(int i = 0; i < MAX; i++) {
        queue[i].genthread = -1;
        queue[i].colors = (Vec *) malloc(sizeof(Vec) * width * height);
        queue[i].samples = (unsigned int *) malloc(sizeof(unsigned int) * width * height);
        queue[i].samplesDiff = (unsigned int *) malloc(sizeof(unsigned int) * width * height);
        queue[i].tdis = (ToDiffInfo *) malloc(sizeof(ToDiffInfo) * width * height);
    }
}

void clear_queue() {
    front = rear;
}

void release_queue() {
    clear_queue();

    for(int i = 0; i < MAX; i++) {
        free(queue[i].colors);
        free(queue[i].samples);
        free(queue[i].samplesDiff);
        free(queue[i].tdis);
    }
}

void put_queue(int thread, Vec *colors, unsigned int *samples, unsigned int *samplesDiff, ToDiffInfo *tdis) {
    if ((rear + 1) % MAX == front){
        LOGI("Queue overflow.");
        return;
    }

    queue[rear].genthread = thread;
    memcpy(queue[rear].colors, colors, sizeof(Vec) * width * height);
    memcpy(queue[rear].samples, samples, sizeof(unsigned int) * width * height);
    memcpy(queue[rear].samplesDiff, samplesDiff, sizeof(unsigned int) * width * height);
    memcpy(queue[rear].tdis, tdis, sizeof(ToDiffInfo) * width * height);

    rear = ++rear % MAX;
}

void get_queue(int *thread, Vec **colors, unsigned int **samples, unsigned int **samplesDiff, ToDiffInfo **tdis) {
    if (front == rear) {
        LOGI("Queue underflow.");
        return;
    }

    *thread = queue[front].genthread;
    *colors = queue[front].colors;
    *samples = queue[front].samples;
    *samplesDiff = queue[front].samplesDiff;
    *tdis = queue[front].tdis;
    /*memcpy(colors, queue[front].colors, sizeof(Vec) * width * height);
    memcpy(samples, queue[front].samples, sizeof(unsigned int) * width * height);
    memcpy(samplesDiff, queue[front].samplesDiff, sizeof(unsigned int) * width * height);
    memcpy(tdis, queue[front].tdis, sizeof(ToDiffInfo) * width * height);*/

    front = ++front % MAX;
}

int peek_thread_queue() {
    if (front == rear) {
        LOGI("Queue underflow.");
        return -1;
    }

    return queue[front].genthread;
}