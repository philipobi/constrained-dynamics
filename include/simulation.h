#pragma once
#include <constraint.h>
#include <globaldef.h>
#include <sparse_linalg.h>

typedef struct simulation {
    int n, m;
    sfloat *q, *dq, *d2q, *Q, *C, *dC, *dJdq, *x, *b, *zeros;
    sparse_matrix *J, *JWJT;
    constraints *constraints;
} simulation;

simulation *init_simulation(int n, int m, sfloat L);

void destruct_simulation(simulation *sim);

void propagate_simulation(simulation *sim, const sfloat dt);

void output_positions(simulation *sim);

void test_simulation(simulation *sim, sfloat dt);