#include <globaldef.h>
#include <simulation.h>
#include <constraint.h>
#include <sparse_linalg.h>
#include <stddef.h>

simulation *
init (int n, int m, sfloat L)
// Initializes Simulation of (n x m) grid with distances L
{
    simulation *sim = NULL;
    constraints *constr = NULL;
    distance_constraint *distance_c = NULL;
    fixpoint_constraint *fixpoint_c = NULL;
    
}