#pragma once
#include <globaldef.h>
#include <sparse_linalg.h>

int minres_solve(const sparse_matrix *A, const sfloat *b, sfloat *x);

int cgs_solve(const sparse_matrix *A, const sfloat *b, sfloat *x);