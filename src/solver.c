#include <globaldef.h>
#include <linalg.h>
#include <solver.h>
#include <sparse_linalg.h>
#include <stdio.h>
#include <stdlib.h>
#define MAXITER 100

int minres_solve(const sparse_matrix *A, const sfloat *b, sfloat *x) {
    // solve A*x = b for x
    // A: (n x n)-matrix
    // x: n-dim. vector
    // b: n-dim: vector

    int i, n = A->n;

    sfloat tol2 = square(EPSILON);

    sfloat *r = NULL, *p = NULL, *s = NULL;

    if (                                         //
        !(r = malloc(n * sizeof(sfloat))) ||     // residual vector
        !(p = malloc(3 * n * sizeof(sfloat))) || // search directions p_k, p_k-1, p_k-2
        !(s = malloc(3 * n * sizeof(sfloat)))    // image vectors s_k, s_k-1, s_k-2
    ) {
        free(r);
        free(p);
        free(s);
        return -1;
    };

    sfloat alpha = 0, beta1 = 0, beta2 = 0, distance = 0;
    sfloat *temp, *Ax0;
    sfloat *p_k = p, *p_k_1 = p + n, *p_k_2 = p + 2 * n;
    sfloat *s_k = s, *s_k_1 = s + n, *s_k_2 = s + 2 * n;

    sfloat s2k[2];
    sfloat *s2k_1 = s2k, *s2k_2 = s2k + 1;

    // #######
    //  k = 0
    // #######
    int k = 0;

    // setup:
    // p_0 = r_0 = b - A * x_0
    sparse_mul_vec(A, x, (Ax0 = s_k));
    for (i = 0; i < n; i++)
        p_k[i] = r[i] = (b[i] - Ax0[i]);

    if ((distance = norm2(r, n)) < tol2)
        goto exit;

    // s_0 = A*p_0
    sparse_mul_vec(A, p_k, s_k);

    // k: 0 -> 1
    swap_ptr(p_k, p_k_1, temp);
    swap_ptr(s_k, s_k_1, temp);
    swap_ptr(s2k_1, s2k_2, temp);
    k++;

    // #######
    //  k = 1
    // #######

    // a_0 = (r_0 * s_0) / (s_0)^2
    alpha = dot(r, s_k_1, n) / (*s2k_1 = norm2(s_k_1, n));
    // x_1 = x_0 + a_0 * p_0
    add_vec_inplace(x, alpha, p_k_1, n);
    // r_1 = r_0 - a_0 * s_0
    add_vec_inplace(r, -alpha, s_k_1, n);

    if ((distance = norm2(r, n)) < tol2)
        goto exit;

    // p_1 <- s_0
    copy_vec(s_k_1, p_k, n);
    // s_1 = A*s_0
    sparse_mul_vec(A, s_k_1, s_k);
    // b_1_1 = (s_1 * s_0) / (s_0)^2
    beta1 = dot(s_k, s_k_1, n) / *s2k_1;
    // p_1 <- p_1 - b_1_1 * p_0
    add_vec_inplace(p_k, -beta1, p_k_1, n);
    // s_1 <- s_1 - b_1_1 * s_0
    add_vec_inplace(s_k, -beta1, s_k_1, n);

    // k: 1 -> 2
    swap_ptr(p_k_1, p_k_2, temp);
    swap_ptr(p_k, p_k_1, temp);
    swap_ptr(s_k_1, s_k_2, temp);
    swap_ptr(s_k, s_k_1, temp);
    swap_ptr(s2k_1, s2k_2, temp);
    k++;

    // #######
    //  k > 1
    // #######

    for ( //
        ; k < MAXITER;
        // k: k -> k+1
        swap_ptr(p_k_1, p_k_2, temp), //
        swap_ptr(p_k, p_k_1, temp),   //
        swap_ptr(s_k_1, s_k_2, temp), //
        swap_ptr(s_k, s_k_1, temp),   //
        swap_ptr(s2k_1, s2k_2, temp), //
        k++                           //
    ) {
        // a_k-1 = (r_k-1 * s_k-1) / (s_k-1)^2
        alpha = dot(r, s_k_1, n) / (*s2k_1 = norm2(s_k_1, n));
        // x_k = x_k-1 + a_k-1 * p_k-1
        add_vec_inplace(x, alpha, p_k_1, n);
        // r_k = r_k-1 - a_k-1 * s_k-1
        add_vec_inplace(r, -alpha, s_k_1, n);

        if ((distance = norm2(r, n)) < tol2)
            break;

        // p_k <- s_k-1
        copy_vec(s_k_1, p_k, n);
        // s_k = A*s_k-1
        sparse_mul_vec(A, s_k_1, s_k);
        // b_k_1 = (s_k * s_k-1) / (s_k-1)^2
        beta1 = dot(s_k, s_k_1, n) / *s2k_1;
        // b_k_2 = (s_k * s_k-2) / (s_k-2)^2
        beta2 = dot(s_k, s_k_2, n) / *s2k_2;
        // p_k <- p_k - b_k_1 * p_k-1 - b_k_2 * p_k-2
        add_2vec_inplace(p_k, -beta1, p_k_1, -beta2, p_k_2, n);
        // s_k <- s_k - b_k_1 * s_k-1 - b_k_2 * s_k-2
        add_2vec_inplace(s_k, -beta1, s_k_1, -beta2, s_k_2, n);

#if DEBUG >= 1
        printf("Iteration: %d\tDistance: %f\n", k, distance);
#endif

#if DEBUG >= 2
        printf("Iteration %d\n", i);
        printf("alpha: %f\n", alpha);
        printf("beta1: %f\n", beta1);
        printf("beta2: %f\n", beta2);
        pprint_array(x, n, 1, "x_k");
        pprint_array(r, n, 1, "r_k");
        pprint_array(p_k, n, 1, "p_k");
        pprint_array(p_k_1, n, 1, "p_k-1");
        pprint_array(p_k_2, n, 1, "p_k-1");
        pprint_array(s_k, n, 1, "s_k");
        pprint_array(s_k_1, n, 1, "s_k-1");
        pprint_array(s_k_2, n, 1, "s_k-1");
#endif
    }

exit:
    free(r);
    free(p);
    free(s);
    return k;
}

int cgs_solve(const sparse_matrix *A, const sfloat *b, sfloat *x) {
    sfloat tol2 = square(EPSILON);
    int i, k, n = A->n;

    sfloat *mem;
    sfloat *r, *rt, *p, *q, *u, *hat;

    if (!(mem = malloc(6 * n * sizeof(sfloat)))) {
        return -1;
    }

    sfloat *vhat, *uhat, *qhat, *Ax0;

    r = mem;
    rt = r + n;
    p = rt + n;
    q = p + n;
    u = q + n;
    hat = u + n;

    vhat = uhat = hat;
    qhat = u;
    Ax0 = r;

    sfloat beta, alpha;
    sfloat rho[2];
    sfloat *rho1 = rho, *rho2 = rho + 1, *temp;

    sparse_mul_vec(A, x, Ax0);
    for (i = 0; i < n; i++)
        rt[i] = r[i] = (b[i] - Ax0[i]);

    for (                   //
        k = 1; k < MAXITER; //
        // k -> k+1
        k++,                       //
        swap_ptr(rho1, rho2, temp) //
    ) {
        *rho1 = dot(rt, r, n);
        if (*rho1 == 0) {
            k = 0;
            break;
        }

        if (k == 1)
            for (i = 0; i < n; i++)
                p[i] = u[i] = r[i];
        else {
            beta = *rho1 / *rho2;
            for (i = 0; i < n; i++)
                p[i] = (u[i] = r[i] + beta * q[i]) + beta * (q[i] + beta * p[i]);
        }

        sparse_mul_vec(A, p, vhat);

        alpha = *rho1 / dot(rt, vhat, n);

        for (i = 0; i < n; i++) {
            q[i] = u[i] - alpha * vhat[i];
            uhat[i] = u[i] + q[i];
            x[i] += alpha * uhat[i];
        }

        sparse_mul_vec(A, uhat, qhat);

        add_vec_inplace(r, -alpha, qhat, n);

        if (norm2(r, n) < tol2)
            break;
    }

    free(mem);
    return k;
}