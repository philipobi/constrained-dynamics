#include <globaldef.h>
#include <linalg.h>
#include <solver.h>
#include <sparse_linalg.h>
#include <stdio.h>
#include <stdlib.h>

void test_parse() {

    sfloat *p_arr = NULL, *p_b = NULL, *p_x = NULL;
    sparse_matrix *p_A = NULL;

    FILE *f_A = fopen("matrix.txt", "r");
    FILE *f_b = fopen("vector.txt", "r");
    if (!f_A || !f_b)
        goto cleanup;

    int n, m, i, j;
    parse_array(f_A, &p_arr, &n, &m);
    parse_array(f_b, &p_b, &i, &j);

    if (m != j) {
        printf("Dimensions dont match\n");
        goto cleanup;
    }

    if (!p_arr || !p_b)
        goto cleanup;

    print_array(p_arr, n, m, "Matrix A");
    print_array(p_b, j, i, "Vector b");

    array_to_sparse(p_arr, n, m, &p_A);
    p_x = calloc(m, sizeof(sfloat));
    if (!p_A || !p_x)
        goto cleanup;

    if (minres_solve(p_A, p_b, p_x) == ERROR) {
        printf("solver error\n");
        goto cleanup;
    }

    sparse_mul_vec(p_A, p_x, p_b);
    print_array(p_x, p_A->m, 1, "Solution x");
    print_array(p_b, p_A->n, 1, "A*x");

cleanup:
    free(p_arr);
    free(p_b);
    free(p_x);
    destruct_sparse(p_A);
    fclose(f_A);
    fclose(f_b);
}

int main() {
    test_parse();
    return 0;
}