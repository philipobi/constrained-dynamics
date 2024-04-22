#include <constraint.h>
#include <globaldef.h>
#include <linalg.h>
#include <simulation.h>
#include <sparse_linalg.h>

void
eval_constraints (simulation *sim)
{
    sfloat diff[3];

    int i, n;
    sparse_row *p_row = sim->J->rows;
    sfloat *p_C = sim->C, *p_dJdq = sim->dJdq;

    const distance_constraint *p_distance_constraint
        = sim->constraints->distance_constraints;
    for (i = 0, n = sim->constraints->n_distance_c; i < n; i++)
        eval_distance_constraint (p_distance_constraint++, sim->q, sim->dq,
                                  p_row++, p_C++, p_dJdq++, diff);

    const fixpoint_constraint *p_fixpoint_constraint
        = sim->constraints->fixpoint_constraints;
    for (i = 0, n = sim->constraints->n_fixpoint_c; i < n; i++)
        eval_fixpoint_constraint (p_fixpoint_constraint++, sim->q, sim->dq,
                                  p_row++, p_C++, p_dJdq++, diff);
}

static void
eval_distance_constraint (const distance_constraint *constraint,
                          const sfloat *p_q, const sfloat *p_dq,
                          sparse_row *J_row, sfloat *p_C, sfloat *p_dJdq,
                          sfloat *diff)
{
    int i, q1 = constraint->q1, q2 = constraint->q2;
    sfloat *row_data1 = J_row->data, *row_data2 = row_data1 + 3,
           *row_col1 = J_row->col, *row_col2 = row_col1 + 3, *p_diff;
    const sfloat *p_q1 = p_q + q1, *p_q2 = p_q + q2;

    for (i = 0, p_diff = diff; i < 3; i++)
        {
            *row_col1++ = q1 + i;
            *row_col2++ = q2 + i;
            *row_data2++
                = -(*row_data1++ = 2 * (*p_diff++ = *p_q1++ - *p_q2++));
        }

    *p_C = norm2 (diff, 3) - square (constraint->L);

    p_q1 = p_dq + q1;
    p_q2 = p_dq + q2;
    for (i = 0, p_diff = diff; i < 3; i++)
        *p_diff++ = *p_q1++ - *p_q2++;
    *p_dJdq = 2 * norm2 (diff, 3);
}

static void
eval_fixpoint_constraint (const fixpoint_constraint *constraint,
                          const sfloat *p_q, const sfloat *p_dq,
                          sparse_row *J_row, sfloat *p_C, sfloat *p_dJdq,
                          sfloat *diff)
{
    int i, q1 = constraint->q1;
    sfloat *row_data = J_row->data, *row_col = J_row->col, *p_diff;
    const sfloat *p_q1 = p_q + q1, *p_point = constraint->point;

    for (i = 0, p_diff = diff; i < 3; i++)
        {
            *row_col++ = q1 + i;
            *row_data++ = 2 * (*p_diff++ = *p_q1++ - *p_point++);
        }

    *p_C = norm2 (diff, 3);

    *p_dJdq = 2 * norm2 (p_dq + q1, 3);
}