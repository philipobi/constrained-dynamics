#pragma once
#include <constraint.h>
#include <globaldef.h>
#include <sparse_linalg.h>

typedef struct simulation {
    int n, m;
    sfloat *q, *dq, *C, *dJdq;
    sparse_matrix *J;
    constraints *constraints;
} simulation;

simulation *init_simulation(int n, int m, sfloat L);

void destruct_simulation(simulation *sim);