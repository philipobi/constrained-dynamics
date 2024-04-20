#include <stdio.h>
#include <globaldef.h>
#include <sparse_linalg.h>
#include <linalg.h>
#include <solver.h>
#include <stdlib.h>

void
test_parse (int argc, char **argv)
{
    if (argc < 3) return;
    
    sfloat *p_arr_A = NULL, *p_b = NULL, *p_x0 = NULL, *p_x = NULL;
    sparse_matrix *p_A = NULL;
    
    FILE *f_A = fopen(argv[1], "r");
    FILE *f_b = fopen(argv[2], "r");
    if (!f_A || !f_b) goto cleanup;
    
    int n, m, i, j;
    parse_array (f_A, &p_arr_A, &n, &m);
    parse_array (f_b, &p_b, &i, &j);

    if(!p_arr_A || !p_b) goto cleanup;
    
    array_to_sparse(p_arr_A, n, m, &p_A);
    p_x0 = calloc(m, sizeof(sfloat));
    p_x = malloc(m * sizeof(sfloat));
    if(!p_A || !p_x0 || !p_x) goto cleanup;
    
    if(minres_solve(p_A, p_b, p_x0, p_x) == ERROR) {
        goto cleanup;
        }

    print_array(p_x, p_A->m, 1);

    cleanup:
    free(p_arr_A);
    free(p_b);
    free(p_x0);
    free(p_x);
    destruct_matrix(p_A);
    fclose(f_A);
    fclose(f_b);
}

int main(int argc, char *argv[])
{
    test_parse(argc, argv);
    return 0;
}