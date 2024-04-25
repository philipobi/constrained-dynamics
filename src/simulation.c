#include <constraint.h>
#include <globaldef.h>
#include <linalg.h>
#include <simulation.h>
#include <solver.h>
#include <sparse_linalg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#define K_S 0.1
#define K_D 0.1
#define F_G -9.81

simulation *init_simulation(int n, int m, sfloat L)
// Initializes Simulation of (n x m) grid with distances L
{
    simulation *sim = NULL;
    constraints *constr = NULL;
    int n_constr;
    int n_coord;

    // allocate memory
    if (                                                                                 //
        !(sim = calloc(1, sizeof(simulation))) ||                                        // simulation object
        !(constr = sim->constraints = init_constraints(n * (m - 1) + m * (n - 1), m)) || // constraint objects
        !(sim->q = malloc(sizeof(sfloat) * (n_coord = 3 * n * m))) ||                    // positions
        !(sim->dq = calloc(n_coord, sizeof(sfloat))) ||                                  // initialize velocities to 0
        !(sim->d2q = calloc(n_coord, sizeof(sfloat))) ||                                 // initialize accelerations to 0
        !(sim->Q = malloc(n_coord * sizeof(sfloat))) ||                                  // real forces
        !(sim->C = malloc(sizeof(sfloat) * (n_constr = constr->n_constr))) ||            // constraint values
        !(sim->dC = malloc(n_constr * sizeof(sfloat))) ||                                // derivative of constraint values
        !(sim->dJdq = malloc(n_constr * sizeof(sfloat))) ||                              // dJdq
        !(sim->x = malloc(n_constr * sizeof(sfloat))) ||                                 // lagrange multiplier
        !(sim->zeros = calloc(n_constr, sizeof(sfloat))) ||                              // zeros for lagrange multiplier
        !(sim->b = malloc(n_constr * sizeof(sfloat))) ||                                 // eq right-hand side
        !(sim->J = init_sparse(n_constr, n_coord, 6)) ||                                 // jacobian
        !(sim->JWJT = init_sparse(n_constr, n_constr, 10))                               // jacobian * W * jacobian_T
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

    sfloat *p_vec;

    // setup fixpoint constraints
    fixpoint_constraint *p_fixpoint_c = constr->fixpoint_constraints;
    p_vec = constr->fixpoints;
    for (j = 0; j < constr->n_fixpoint_c; j++, p_fixpoint_c++) {
        p_fixpoint_c->q1 = 3 * j;
        p_fixpoint_c->fixpoint = p_vec;
        *p_vec++ = 0;     // fixpoint_x
        *p_vec++ = j * L; // fixpoint_y
        *p_vec++ = 0;     // fixpoint_z
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

    // setup force vector Q
    p_vec = sim->Q;
    for (i = 0; i < n_coord; i++) {
        *p_vec++ = 0;   // F_x
        *p_vec++ = 0;   // F_y
        *p_vec++ = F_G; // F_z
    }

    sim->n = n;
    sim->m = m;
    return sim;
}

void destruct_simulation(simulation *sim) {
    if (!sim)
        return;
    destruct_constraints(sim->constraints);
    free(sim->q);
    free(sim->dq);
    free(sim->d2q);
    free(sim->Q);
    free(sim->C);
    free(sim->dC);
    free(sim->dJdq);
    free(sim->x);
    free(sim->zeros);
    free(sim->b);
    destruct_sparse(sim->J);
    destruct_sparse(sim->JWJT);
    free(sim);
}

void propagate_simulation(simulation *sim, const sfloat dt) {

    int i, n_coord = 3 * sim->n * sim->m, n_constr = sim->constraints->n_constr;

    // reset matrices to empty
    destruct_sparse(sim->J);
    sim->J = init_sparse(n_constr, n_coord, 6);
    destruct_sparse(sim->JWJT);
    sim->JWJT = init_sparse(n_constr, n_constr, 10);

    if (!sim->J || !sim->JWJT) {
        printf("error");
        return;
    }
    // reset vectors
    copy_vec(sim->zeros, sim->x, n_constr);
    copy_vec(sim->zeros, sim->dC, n_constr);
    copy_vec(sim->zeros, sim->C, n_constr);
    copy_vec(sim->zeros, sim->dJdq, n_constr);
    copy_vec(sim->zeros, sim->b, n_constr);
    for (i = 0; i < n_coord; i++)
        sim->d2q[i] = 0;

    // evaluate constraints for current iteration and build J, C, dJdq
    eval_constraints(sim);

    // J * W * J_T
    sparse_mul_transpose(sim->J, sim->JWJT);

    // dC = J * dq
    sparse_mul_vec(sim->J, sim->dq, sim->dC);

    // b <- J * W * Q
    sparse_mul_vec(sim->J, sim->Q, sim->b);

    // b = - dJdq - J * W * Q - K_S * C - K_D * dC
    // b <- - dJdq -b - K_S * C - K_D * dC
    for (i = 0; i < sim->constraints->n_constr; i++) {
        sim->b[i] = -sim->dJdq[i] - sim->b[i] - K_S * sim->C[i] - K_D * sim->dC[i];
    }

    // solve (J * W * J_T) * x = b
    if (cgs_solve(sim->JWJT, sim->b, sim->x) == -1)
        return;

    // d2q <- Q_c = J_T * x
    sparse_transpose_mul_vec(sim->J, sim->x, sim->d2q);

    // d2q = W * (Q + Q_c)
    // d2q <- Q + d2q
    for (i = 0; i < n_coord; i++) {
        sim->d2q[i] += sim->Q[i];
    }

    // Integration:
    // dq <- dq + d2q * dt
    add_vec_inplace(sim->dq, dt, sim->d2q, n_coord);
    // q <- q + dq * dt
    add_vec_inplace(sim->q, dt, sim->dq, n_coord);

#if DEBUG == 1
    int n_constr = sim->constraints->n_constr;
    pprint_array(sim->b, 1, n_constr, "b");
    pprint_array(sim->x, 1, n_constr, "lambda");
    printf("J\n");
    print_sparse(sim->J);
    printf("JJT\n");
    print_sparse(sim->JWJT);
#endif
}

void output_positions(simulation *sim) { print_array(sim->q, sim->n * sim->m, 3); }

void test_simulation(simulation *sim, sfloat dt) {
    FILE *f_q = NULL, *f_dq = NULL, *f_Q = NULL;

    if (                                  //
        !(f_q = fopen("q.txt", "r")) ||   //
        !(f_dq = fopen("dq.txt", "r")) || //
        !(f_Q = fopen("F.txt", "r"))      //
    ) {
        printf("error occurred\n");
        goto exit;
    }
    int n, m;
    sfloat *q, *dq, *Q;
    parse_array(f_q, &q, &n, &m);
    parse_array(f_dq, &dq, &n, &m);
    parse_array(f_Q, &Q, &n, &m);

    copy_vec(q, sim->q, m);
    copy_vec(dq, sim->dq, m);
    copy_vec(Q, sim->Q, m);

    printf("\n\nNew Iteration\n\n");
    printf("\n\nq\n\n");
    print_array(sim->q, 27, 1);
    printf("\n\ndq\n\n");
    print_array(sim->dq, 27, 1);

    propagate_simulation(sim, dt);

    printf("\n\nq\n\n");
    print_array(sim->q, 27, 1);
    printf("\n\ndq\n\n");
    print_array(sim->dq, 27, 1);

exit:
    free(q);
    free(dq);
    free(Q);
    fclose(f_q);
    fclose(f_dq);
    fclose(f_Q);
}