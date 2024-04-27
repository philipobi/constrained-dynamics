#include <constraint.h>
#include <globaldef.h>
#include <linalg.h>
#include <simulation.h>
#include <sparse_linalg.h>
#include <stdlib.h>

constraints *init_constraints(int n_distance_c, int n_fixpoint_c) {
    constraints *constr = NULL;
    if (                                                                                        //
        !(constr = calloc(1, sizeof(constraints))) ||                                           //
        !(constr->distance_constraints = malloc(n_distance_c * sizeof(distance_constraint))) || //
        !(constr->fixpoint_constraints = malloc(n_fixpoint_c * sizeof(fixpoint_constraint))) || //
        !(constr->fixpoints = malloc(3 * n_fixpoint_c * sizeof(sfloat)))                        //
    ) {
        destruct_constraints(constr);
        return NULL;
    }

    constr->n_constr = (constr->n_distance_c = n_distance_c) + (constr->n_fixpoint_c = n_fixpoint_c);

    return constr;
}

void destruct_constraints(constraints *constr) {
    if (constr == NULL)
        return;
    free(constr->distance_constraints);
    free(constr->fixpoint_constraints);
    free(constr->fixpoints);
    free(constr);
}

static void eval_distance_constraint(const distance_constraint *constraint, const sfloat *p_q, const sfloat *p_dq, sparse_row *J_row, sfloat *p_C,
                                     sfloat *p_dJdq, sfloat *diff) {

    int i, q1 = constraint->q1, q2 = constraint->q2;
    const sfloat *p_q1 = p_q + q1, *p_q2 = p_q + q2;

    for (i = 0; i < 3; i++) {
        J_row->col[i] = q1 + i;
        J_row->col[i + 3] = q2 + i;
        J_row->data[i + 3] = -(J_row->data[i] = 2 * (diff[i] = (p_q1[i] - p_q2[i])));
    }
    J_row->nnz = 6;

    *p_C = norm2(diff, 3) - square(constraint->L);

    p_q1 = p_dq + q1;
    p_q2 = p_dq + q2;
    for (i = 0; i < 3; i++)
        diff[i] = p_q1[i] - p_q2[i];

    *p_dJdq = 2 * norm2(diff, 3);
}

static void eval_fixpoint_constraint(const fixpoint_constraint *constraint, const sfloat *p_q, const sfloat *p_dq, sparse_row *J_row, sfloat *p_C,
                                     sfloat *p_dJdq, sfloat *diff) {
    int i, q1 = constraint->q1;
    const sfloat *p_q1 = p_q + q1;

    for (i = 0; i < 3; i++) {
        J_row->col[i] = q1 + i;
        J_row->data[i] = 2 * (diff[i] = p_q1[i] - constraint->fixpoint[i]);
    }
    J_row->nnz = 3;

    *p_C = norm2(diff, 3);

    *p_dJdq = 2 * norm2(p_dq + q1, 3);
}

void eval_constraints(simulation *sim) {
    sfloat diff[3];

    int i;
    sparse_row *p_row = sim->J->rows;
    sfloat *p_C = sim->C, *p_dJdq = sim->dJdq;

    const distance_constraint *p_distance_constraint = sim->constraints->distance_constraints;
    for (i = 0; i < sim->constraints->n_distance_c; i++)
        eval_distance_constraint(p_distance_constraint++, sim->q, sim->dq, p_row++, p_C++, p_dJdq++, diff);

    const fixpoint_constraint *p_fixpoint_constraint = sim->constraints->fixpoint_constraints;
    for (i = 0; i < sim->constraints->n_fixpoint_c; i++)
        eval_fixpoint_constraint(p_fixpoint_constraint++, sim->q, sim->dq, p_row++, p_C++, p_dJdq++, diff);
}