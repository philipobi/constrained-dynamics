#include <constraint.h>
#include <globaldef.h>
#include <linalg.h>
#include <simulation.h>
#include <solver.h>
#include <sparse_linalg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

void run_simulation() {
    simulation *sim = init_simulation(3, 3, 5);
    if (!sim)
        return;

    clock_t start = clock();
    for (int i = 0; i < 10000; i++) {
        propagate_simulation(sim, .01);
        output_positions(sim);
        // usleep(1e6);
    }

    float seconds = (float)(clock() - start) / CLOCKS_PER_SEC;
    // printf("%f seconds\n", seconds);

    destruct_simulation(sim);
}

void test() {
    int n, m, k, l;
    sfloat *A_ = NULL, *b = NULL, *x = NULL;
    FILE *f_A = NULL, *f_b = NULL;
    sparse_matrix *A = NULL;
    if (!(f_A = fopen("matrix.txt", "r")) || !(f_b = fopen("vector.txt", "r")))
        goto exit;
    parse_array(f_A, &A_, &n, &m);
    parse_array(f_b, &b, &k, &l);

    if (                                       //
        !A_ ||                                 //
        !(n * m * k * l) ||                    //
        !(n == m && m == l && k == 1) ||       //
        !(x = calloc(m, sizeof(sfloat))) ||    //
        !(array_to_sparse(A_, n, m, &A), A) || //
        cgs_solve(A, b, NULL, x) == -1)        //
    {
        printf("Error occurred.\n");
        goto exit;
    }

    print_array(x, m, 1);

exit:
    free(A_);
    free(b);
    free(x);
    fclose(f_A);
    fclose(f_b);
    destruct_sparse(A);
}

void test_sim() {
    simulation *sim = NULL;
    if (!(sim = init_simulation(3, 3, 5)))
        goto exit;

    while (1) {
        test_simulation(sim, 0.01);
    }

exit:
    destruct_simulation(sim);
}

int main() {
    run_simulation();
    return 0;
}