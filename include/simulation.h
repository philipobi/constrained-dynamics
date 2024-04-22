#pragma once
#include <constraint.h>
#include <globaldef.h>
#include <sparse_linalg.h>

typedef struct simulation {
    int n, m;
    sfloat *q, *dq, *d2q, *Q, *C, *dC, *dJdq, *lambda, *b;
    sparse_matrix *J, *JJT;
    constraints *constraints;
} simulation;

simulation *init_simulation(int n, int m, sfloat L);

void destruct_simulation(simulation *sim);