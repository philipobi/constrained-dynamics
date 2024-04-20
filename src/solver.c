#include <globaldef.h>
#include <linalg.h>
#include <sparse_linalg.h>
#include <solver.h>
#include <stdlib.h>
#include <stdio.h>

int
minres_solve (const sparse_matrix *p_A, const sfloat *p_b, const sfloat *p_x0,
              sfloat *p_x)
{
    
    int n = p_A->n;

    sfloat *p_r = NULL, *p_p = NULL, *p_Ar = NULL, *p_Ap = NULL;

    p_r = malloc (n * sizeof (sfloat));
    p_p = malloc (n * sizeof (sfloat));
    p_Ar = malloc (n * sizeof (sfloat));
    p_Ap = malloc (n * sizeof (sfloat));

    if (!p_r || !p_p || !p_Ar || !p_Ap)
        {
            free (p_r);
            free (p_p);
            free (p_Ar);
            free (p_Ap);
            return ERROR;
        }

    sfloat *p_vec, *p_vec_, r_A_r, r1_A_r1, a, b;

    const sfloat *p_vec1, *p_vec2;
    int i;
    
    // p_0 = r_0 = b - A * x_0
    sparse_mul_vec (p_A, p_x0, p_Ar);
    for (i = 0, p_vec_ = p_p, p_vec = p_r, p_vec1 = p_b, p_vec2 = p_Ar; i < n;
         i++)
        *p_vec_++ = *p_vec++ = (*p_vec1++ - *p_vec2++);

    // a_0 = (r_0 * A * r_0)/(A * p_0)^2
    sparse_mul_vec (p_A, p_r, p_Ar);
    r_A_r = dot (p_r, p_Ar, n);
    a = r_A_r / norm2 (p_Ar, n);

    // r_1 = r_0 - a_0 * A * p_0
    for (i = 0, p_vec = p_r, p_vec1 = p_Ar; i < n; i++)
        *p_vec++ -= a * *p_vec1++;

    // b_0 = (r_1 * A * r_1)/(r_0 * A * r_0)
    sparse_mul_vec (p_A, p_r, p_Ar);
    r1_A_r1 = dot (p_r, p_Ar, n);
    b = r1_A_r1 / r_A_r;

    // x_1 = x_0 + a_0 * p_0
    for (i = 0, p_vec = p_x, p_vec1 = p_x0, p_vec2 = p_p; i < n; i++)
        *p_vec++ = *p_vec1++ + a * *p_vec2++;

    // p_1 = r_1 + b_0 * p_0
    for (i = 0, p_vec = p_p, p_vec1 = p_r; i < n; i++, p_vec++)
        *p_vec = *p_vec1++ * b * *p_vec;

    r_A_r = r1_A_r1;

    int iteration = 1;
    sfloat distance;
    while ((distance = norm2 (p_r, n)) > EPSILON)
        {
            // a_k = (r_k * A * r_k)/(A * p_k)^2
            sparse_mul_vec (p_A, p_p, p_Ap);
            a = r_A_r / norm2 (p_Ap, n);

            // r_k+1 = r_k - a_k * A * p_k
            for (i = 0, p_vec = p_r, p_vec1 = p_Ap; i < n; i++)
                *p_vec++ -= a * *p_vec1++;

            // b_k = (r_k+1 * A * r_k+1)/(r_k * A * r_k)
            sparse_mul_vec (p_A, p_r, p_Ar);
            r1_A_r1 = dot (p_r, p_Ar, n);
            b = r1_A_r1 / r_A_r;

            // x_k+1 = x_k + a_k * p_k
            for (i = 0, p_vec = p_x, p_vec1 = p_p; i < n; i++)
                *p_vec++ += a * *p_vec1++;

            // p_k+1 = r_k+1 + b_k * p_k
            for (i = 0, p_vec = p_p, p_vec1 = p_r; i < n; i++, p_vec++)
                *p_vec = *p_vec1++ * b * *p_vec;

            r_A_r = r1_A_r1;

            iteration++;
            printf("Iteration: %d, Distance: %f\n", iteration, distance);
        }

    free (p_r);
    free (p_p);
    free (p_Ar);
    free (p_Ap);
    return SUCCESS;
}