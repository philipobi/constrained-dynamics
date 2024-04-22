#include <constraint.h>
#include <globaldef.h>
#include <simulation.h>
#include <sparse_linalg.h>
#include <stddef.h>

simulation *init_simulation(int n, int m, sfloat L)
// Initializes Simulation of (n x m) grid with distances L
{
    simulation *sim = NULL;
    constraints *constr = NULL;
    sfloat *p_vec;

    // allocate memory
    if (!(sim = calloc(1, sizeof(simulation))) || !(constr = sim->constraints = init_constraints(n * (m - 1) + m * (n - 1), m)) ||
        !(sim->q = malloc(3 * n * m * sizeof(sfloat))) || !(sim->dq = calloc(3 * n * m, sizeof(sfloat))) ||
        !(sim->C = malloc(constr->n_constr * sizeof(sfloat))) || !(sim->dJdq = malloc(constr->n_constr * sizeof(sfloat))) ||
        !(sim->J = init_sparse(constr->n_constr, 3 * n * m, 6))) {
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

    return sim;
}

void destruct_simulation(simulation *sim) {
    if (!sim)
        return;
    free(sim->q);
    free(sim->dq);
    free(sim->C);
    free(sim->dJdq);
    destruct_constraints(sim->constraints);
    destruct_sparse(sim->J);
}