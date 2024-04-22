#include <globaldef.h>
#include <linalg.h>
#include <solver.h>
#include <sparse_linalg.h>
#include <stdio.h>
#include <stdlib.h>
#define MAXITER 100

int
minres_solve (const sparse_matrix *A, const sfloat *b, sfloat *x)
{
    int n = A->n;

    if (n != A->m) return ERROR;

    sfloat tol2 = square(EPSILON);

    sfloat *r = NULL, *p = NULL, *s = NULL;

    r = calloc (n, sizeof (sfloat));
    p = calloc (3 * n, sizeof (sfloat));
    s = calloc (3 * n, sizeof (sfloat));

    if (!r || !p || !s)
        {
            free (r);
            free (p);
            free (s);
            return ERROR;
        };

    sfloat alpha = 0, beta1 = 0, beta2 = 0, distance = 0, *temp, *Ax0,
        *p_k = p, *p_k_1 = p + n, *p_k_2 = p + 2 * n, *s_k = s, *s_k_1 = s + n,
        *s_k_2 = s + 2 * n;

    sfloat s2k[2];
    sfloat *s2k_1 = s2k, *s2k_2 = s2k + 1;

    // setup:
    // p_0 = r_0 = b - A * x_0
    sparse_mul_vec (A, x, (Ax0 = s_k));
    int i = 0;
    for (sfloat *p_0 = p_k, *r_0 = r; i < n; i++)
        *p_0++ = *r_0++ = *b++ - *Ax0++;

    sparse_mul_vec (A, p_k, s_k);
    
    // k: 0 -> 1
    swap_ptr (p_k, p_k_1, temp);
    swap_ptr (s_k, s_k_1, temp);
    swap_ptr (s2k_1, s2k_2, temp);

    // k = 1
    alpha = dot (r, s_k_1, n) / (*s2k_1 = norm2 (s_k_1, n));
    add_vec_inplace (x, alpha, p_k_1, n);
    add_vec_inplace (r, -alpha, s_k_1, n);
    copy_vec (s_k_1, p_k, n);
    sparse_mul_vec (A, s_k_1, s_k);
    beta1 = dot (s_k, s_k_1, n) / *s2k_1;
    add_vec_inplace (p_k, -beta1, p_k_1, n);
    add_vec_inplace (s_k, -beta1, s_k_1, n);
    
    // k: 1 -> 2
    swap_ptr (p_k_1, p_k_2, temp);
    swap_ptr (p_k, p_k_1, temp);
    swap_ptr (s_k_1, s_k_2, temp);
    swap_ptr (s_k, s_k_1, temp);
    swap_ptr (s2k_1, s2k_2, temp);

    // k > 1
    for (int i = 2; i < MAXITER; i++)
        {
            alpha = dot (r, s_k_1, n) / (*s2k_1 = norm2 (s_k_1, n));
            add_vec_inplace (x, alpha, p_k_1, n);
            add_vec_inplace (r, -alpha, s_k_1, n);
            if ((distance = norm2 (r, n)) < tol2)
                break;
            copy_vec (s_k_1, p_k, n);
            sparse_mul_vec (A, s_k_1, s_k);
            beta1 = dot (s_k, s_k_1, n) / *s2k_1;
            beta2 = dot (s_k, s_k_2, n) / *s2k_2;
            add_2vec_inplace (p_k, -beta1, p_k_1, -beta2, p_k_2, n);
            add_2vec_inplace (s_k, -beta1, s_k_1, -beta2, s_k_2, n);

            #if DEBUG == 1
            printf("Iteration: %d\tDistance: %f\n", i, distance);
            #endif
            
            #if DEBUG == 2
            printf("Iteration %d\n", i);
            printf("alpha: %f\n", alpha);
            printf("beta1: %f\n", beta1);
            printf("beta2: %f\n", beta2);
            print_array(x, 1, n, "x_k");
            print_array(r, 1, n, "r_k");
            print_array(p_k, 1, n, "p_k");
            print_array(p_k_1, 1, n, "p_k-1");
            print_array(p_k_2, 1, n, "p_k-1");
            print_array(s_k, 1, n, "s_k");
            print_array(s_k_1, 1, n, "s_k-1");
            print_array(s_k_2, 1, n, "s_k-1");
            #endif

            // k: k -> k + 1
            swap_ptr (p_k_1, p_k_2, temp);
            swap_ptr (p_k, p_k_1, temp);
            swap_ptr (s_k_1, s_k_2, temp);
            swap_ptr (s_k, s_k_1, temp);
            swap_ptr (s2k_1, s2k_2, temp);
        }

    free (r);
    free (p);
    free (s);
    return SUCCESS;
}