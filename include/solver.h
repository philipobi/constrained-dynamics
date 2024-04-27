#pragma once
#include <globaldef.h>
#include <sparse_linalg.h>

int cr_solve(const sparse_matrix *A, const sfloat *b, const sfloat *x_0, sfloat *x);

int cgs_solve(const sparse_matrix *A, const sfloat *b, const sfloat *x_0, sfloat *x);