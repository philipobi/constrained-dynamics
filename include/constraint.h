#pragma once
#include <globaldef.h>
#include <simulation.h>
#include <sparse_linalg.h>

typedef struct constraints
{
    distance_constraint *distance_constraints;
    fixpoint_constraint *fixpoint_constraints;
    int n_distance_c, n_fixpoint_c;
} constraints;

void eval_constraints (simulation *sim);

typedef struct distance_constraint
{
    sfloat L;
    int q1, q2;
} distance_constraint;

typedef struct fixpoint_constraint
{
    int q1;
    sfloat *point;
} fixpoint_constraint;