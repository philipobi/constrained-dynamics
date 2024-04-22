#include <constraint.h>
#include <globaldef.h>
#include <linalg.h>
#include <simulation.h>
#include <solver.h>
#include <sparse_linalg.h>
#include <stddef.h>
#include <stdlib.h>

#define KS 0
#define KD 0

simulation *init_simulation(int n, int m, sfloat L)
// Initializes Simulation of (n x m) grid with distances L
{
    simulation *sim = NULL;
    constraints *constr = NULL;
    sfloat *p_vec;
    int n_constr;

    // allocate memory
    if (                                                                                 //
        !(sim = calloc(1, sizeof(simulation))) ||                                        //
        !(constr = sim->constraints = init_constraints(n * (m - 1) + m * (n - 1), m)) || //
        !(sim->q = malloc(3 * n * m * sizeof(sfloat))) ||                                //
        !(sim->dq = calloc(3 * n * m, sizeof(sfloat))) ||                                //
        !(sim->d2q = calloc(3 * n * m, sizeof(sfloat))) ||                               //
        !(sim->Q = calloc(3 * n * m, sizeof(sfloat))) ||                                 //
        !(sim->C = malloc(n_constr = constr->n_constr * sizeof(sfloat))) ||              //
        !(sim->dC = malloc(n_constr * sizeof(sfloat))) ||                                //
        !(sim->dJdq = malloc(n_constr * sizeof(sfloat))) ||                              //
        !(sim->lambda = calloc(n_constr, sizeof(sfloat))) ||                             //
        !(sim->b = malloc(n_constr * sizeof(sfloat))) ||                                 //
        !(sim->J = init_sparse(n_constr, 3 * n * m, 6)) ||                               //
        !(sim->JJT = init_sparse(n_constr, n_constr, 10))                                //
    ) {
        destruct_simulation(sim);
        return NULL;
    };

    int i, j;
    distance_constraint *p_distance_c = constr->distance_constraints;
    // setup horizontal constraints
    for (i = 0; i < n; i++) {
        for (j = 0; j < m - 1; j++, p_distance_c++) {
            p_distance_c->L = L;
            p_distance_c->q2 = (p_distance_c->q1 = 3 * (i * m + j)) + 3;
        }
    }
    // setup vertical constraints
    for (j = 0; j < m; j++) {
        for (i = 0; i < n - 1; i++, p_distance_c++) {
            p_distance_c->L = L;
            p_distance_c->q2 = (p_distance_c->q1 = 3 * (i * m + j)) + 3 * m;
        }
    }

    // setup fixpoint constraints
    fixpoint_constraint *p_fixpoint_c = constr->fixpoint_constraints;
    for (j = 0, p_vec = constr->fixpoints; j < m; j++) {
        p_fixpoint_c->q1 = 3 * j;
        p_fixpoint_c->point = p_vec;
        *p_vec++ = 0;
        *p_vec++ = j * L;
        *p_vec++ = 0;
    }

    // setup location vector q
    p_vec = sim->q;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            *p_vec++ = i * L; // x
            *p_vec++ = j * L; // y
            *p_vec++ = 0;     // z
        }
    }

    sim->n = n;
    sim->m = m;
    return sim;
}

void destruct_simulation(simulation *sim) {
    if (!sim)
        return;
    free(sim->q);
    free(sim->dq);
    free(sim->d2q);
    free(sim->Q);
    free(sim->C);
    free(sim->dC);
    free(sim->dJdq);
    free(sim->lambda);
    free(sim->b);
    destruct_constraints(sim->constraints);
    destruct_sparse(sim->J);
    destruct_sparse(sim->JJT);
    free(sim);
}

void propagate_simulation(simulation *sim, sfloat dt) {

    int n = sim->n, m = sim->m;
    reset_sparse(sim->J);
    reset_sparse(sim->JJT);
    eval_constraints(sim);
    sparse_mul_transpose(sim->J, sim->JJT);
    sparse_mul_vec(sim->J, sim->dq, sim->dC);
    sfloat *p_JQ = sim->b;
    sparse_mul_vec(sim->J, sim->Q, p_JQ);

    sfloat *p_dJdq = sim->dJdq, *p_C = sim->C, *p_dC = sim->dC;
    for (int i = 0; i < sim->constraints->n_constr; i++, p_JQ++) {
        (*p_JQ) = -(*p_dJdq++) - (*p_JQ) - KS * (*p_C++) - KD * (*p_dC++);
    }

    if (minres_solve(sim->JJT, sim->b, sim->lambda) == ERROR)
        return;

    sparse_transpose_mul_vec(sim->J, sim->lambda, sim->d2q);
    add_vec_inplace(sim->d2q, 1, sim->Q, 3 * n * m);
    add_vec_inplace(sim->dq, dt, sim->d2q, 3 * n * m);
    add_vec_inplace(sim->q, dt, sim->dq, 3 * n * m);
}