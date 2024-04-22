#include <stdio.h>
#include <stdlib.h>

typedef struct mystruct {
    int *ptr1, *ptr2;
} mystruct;

int main() {
    mystruct *p_test = calloc(1, sizeof(mystruct));
    if (p_test->ptr2 == NULL)
        printf("nullptr\n");
}