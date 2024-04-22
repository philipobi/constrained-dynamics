#pragma once
#include <globaldef.h>
#include <sparse_linalg.h>
#include <constraint.h>

typedef struct simulation {
    sfloat *q, *dq, *C, *dJdq;
    sparse_matrix *J;
    constraints *constraints;
} simulation;