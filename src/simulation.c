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
#define F_G 9.81

simulation *init_simulation(int n, int m, sfloat L)
// Initializes Simulation of (n x m) grid with distances L
{
    simulation *sim = NULL;
    constraints *constr = NULL;
    sfloat *mem = NULL;
    int n_constr;
    int n_coord = 3 * n * m;
    size_t size = sizeof(sfloat);

    // allocate memory
    if (                                                                                 //
        !(sim = calloc(1, sizeof(simulation))) ||                                        // simulation object
        !(constr = sim->constraints = init_constraints(n * (m - 1) + m * (n - 1), m)) || // constraint objects

        !(sim->q = malloc(n_coord * size)) ||   // positions
        !(sim->dq = malloc(n_coord * size)) ||  // velocities
        !(sim->d2q = malloc(n_coord * size)) || // accelerations
        !(sim->Q = malloc(n_coord * size)) ||   // real forces

        !(sim->C = malloc((n_constr = constr->n_constr) * size)) || // constraint values
        !(sim->dC = malloc(n_constr * size)) ||                     // derivative of constraint values
        !(sim->dJdq = malloc(n_constr * size)) ||                   // dJdq
        !(sim->b = malloc(n_constr * size)) ||                      // eq right-hand side
        !(sim->x = malloc(n_constr * size)) ||                      // lagrange multiplier

        !(sim->J = init_sparse(n_constr, n_coord, 6)) ||   // jacobian
        !(sim->JWJT = init_sparse(n_constr, n_constr, 10)) // jacobian * W * jacobian_T
    ) {
        destruct_simulation(sim);
        return NULL;
    };

    sim->n = n;
    sim->m = m;

    sim->Q_c = sim->d2q;
    sim->JWQ = sim->b;

    int i, j;
    sfloat *p_vec;

    // setup location vector q
    p_vec = sim->q;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            *p_vec++ = i * L; // x
            *p_vec++ = j * L; // y
            *p_vec++ = 0;     // z
        }
    }

    // set velocities to zero
    reset_vec(sim->dq, n_coord);

    // set accelerations to zero
    reset_vec(sim->d2q, n_coord);

    // setup force vector Q
    p_vec = sim->Q;
    for (i = 0; i < n_coord; i++) {
        *p_vec++ = 0;    // F_x
        *p_vec++ = 0;    // F_y
        *p_vec++ = -F_G; // F_z
    }

    // setup distance constraints
    distance_constraint *p_distance_c = constr->distance_constraints;
    // horizontal
    for (i = 0; i < n; i++) {
        for (j = 0; j < m - 1; j++, p_distance_c++) {
            p_distance_c->L = L;
            p_distance_c->q2 = (p_distance_c->q1 = 3 * (i * m + j)) + 3;
        }
    }
    // vertical
    for (j = 0; j < m; j++) {
        for (i = 0; i < n - 1; i++, p_distance_c++) {
            p_distance_c->L = L;
            p_distance_c->q2 = (p_distance_c->q1 = 3 * (i * m + j)) + 3 * m;
        }
    }

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
    free(sim->b);
    free(sim->x);

    destruct_sparse(sim->J);
    destruct_sparse(sim->JWJT);

    free(sim);
}

void propagate_simulation(simulation *sim, const sfloat dt) {

    int i;
    int n_coord = 3 * sim->n * sim->m;
    int n_constr = sim->constraints->n_constr;

    // reset matrices to empty
    reset_sparse(sim->J);
    reset_sparse(sim->JWJT);

    // evaluate constraints for current iteration and build J, C, dJdq
    eval_constraints(sim);

    // J * W * J_T
    sparse_mul_transpose(sim->J, sim->JWJT);

    // dC = J * dq
    sparse_mul_vec(sim->J, sim->dq, sim->dC);

    // J * W * Q
    sparse_mul_vec(sim->J, sim->Q, sim->JWQ);

    // b = - dJdq - J * W * Q - K_S * C - K_D * dC
    for (i = 0; i < sim->constraints->n_constr; i++) {
        sim->b[i] = -sim->dJdq[i] - sim->JWQ[i] - K_S * sim->C[i] - K_D * sim->dC[i];
    }

    // solve (J * W * J_T) * x = b
    if (cgs_solve(sim->JWJT, sim->b, NULL, sim->x) == -1)
        return;

    // Q_c = J_T * x
    sparse_transpose_mul_vec(sim->J, sim->x, sim->Q_c);

    for (i = 0; i < n_coord; i++) {
        // d2q = W * (Q + Q_c)
        sim->d2q[i] = sim->Q[i] + sim->Q_c[i];

        // Integration:
        // dq <- dq + d2q * dt
        sim->dq[i] += dt * sim->d2q[i];
        // q <- q + dq * dt
        sim->q[i] += dt * sim->dq[i];
    }

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