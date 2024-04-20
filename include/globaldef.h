#pragma once

#define min(a, b) ((a) < (b) ? (a) : (b))
#define square(a) ((a) * (a))
#define absolute(a) ((a) > 0 ? (a) : -(a))
#define EPSILON 1e-6

typedef double sfloat;

enum {
    SUCCESS,
    ERROR
};