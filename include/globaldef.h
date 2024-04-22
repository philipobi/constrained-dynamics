#pragma once

#define min(a, b) ((a) < (b) ? (a) : (b))
#define square(a) ((a) * (a))
#define absolute(a) ((a) > 0 ? (a) : -(a))
#define swap_ptr(a, b, temp)                                                  \
    temp = b;                                                                 \
    b = a;                                                                    \
    a = temp;
#define EPSILON 1e-6
#define DEBUG 1

typedef double sfloat;

enum {
    SUCCESS,
    ERROR
};