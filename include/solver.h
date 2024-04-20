#pragma once
#include <globaldef.h>
#include <sparse_linalg.h>

int minres_solve (const sparse_matrix *p_A, const sfloat *p_b,
                  const sfloat *p_x0, sfloat *p_x);