#include <globaldef.h>
#include <linalg.h>
#include <sparse_linalg.h>
#include <stdio.h>
#include <stdlib.h>

static void destruct_rows(sparse_row *p_row, int n) {
    for (int i = 0; i < n; i++, p_row++) {
        free(p_row->col);
        p_row->col = NULL;
        free(p_row->data);
        p_row->data = NULL;
    }
}

static void free_matrix(sparse_matrix *p_matrix) {
    free(p_matrix->rows);
    p_matrix->rows = NULL;
    free(p_matrix);
}

void destruct_sparse(sparse_matrix *p_matrix) {
    if (!p_matrix)
        return;
    destruct_rows(p_matrix->rows, p_matrix->n);
    free_matrix(p_matrix);
}

sparse_matrix *init_sparse(int n, int m, int max_nnz) {
    sparse_matrix *p_matrix = malloc(sizeof(sparse_matrix));
    sparse_row *rows = malloc(n * sizeof(sparse_row));
    if (!p_matrix || !rows) {
        free(p_matrix);
        free(rows);
        return NULL;
    }

    p_matrix->n = n;
    p_matrix->m = m;
    p_matrix->rows = rows;

    int i;
    sparse_row *p_row;
    for (i = 0, p_row = p_matrix->rows; i < n; i++, p_row++) {
        p_row->nnz = 0;
        p_row->col = malloc(max_nnz * sizeof(int));
        p_row->data = malloc(max_nnz * sizeof(sfloat));
        if (!p_row->col || !p_row->data) {
            destruct_rows(p_matrix->rows, i + 1);
            free_matrix(p_matrix);
            return NULL;
        }
    }

    return p_matrix;
}

void reset_sparse(sparse_matrix *p_matrix) {
    int i;
    sparse_row *p_row;
    for (i = 0, p_row = p_matrix->rows; i < p_matrix->n; i++, p_row++)
        p_row->nnz = 0;
}

void sparse_to_array(const sparse_matrix *p_matrix, sfloat *arr) {
    int i, k;
    sparse_row *p_row;
    for (i = 0, p_row = p_matrix->rows; i < p_matrix->n; i++, p_row++) {
        for (k = 0; k < p_row->nnz; k++)
            arr[i * p_matrix->m + p_row->col[k]] = p_row->data[k];
    }
}

void print_sparse(const sparse_matrix *p_matrix) {
    sfloat *arr = calloc(p_matrix->n * p_matrix->m, sizeof(sfloat));
    if (!arr)
        return;
    sparse_to_array(p_matrix, arr);
    print_array(arr, p_matrix->n, p_matrix->m, "");
    free(arr);
}

static sfloat sparse_mul_rows(const sparse_row *p_row_i, const sparse_row *p_row_j) {
    int n, k, m = 0;
    sfloat *data_i, *data_j, sum = 0;

    for (n = 0; n < p_row_i->nnz; n += BLKSIZE) {
        for (; m < p_row_j->nnz; m += BLKSIZE) {
            if (p_row_i->col[n] < p_row_j->col[m])
                break;
            if (p_row_i->col[n] > p_row_j->col[m])
                continue;
            data_i = p_row_i->data + n;
            data_j = p_row_j->data + m;
            for (k = 0; k < BLKSIZE; k++)
                sum += (*data_i++) * (*data_j++);
        }
    }

    return sum;
}

void sparse_mul_transpose(const sparse_matrix *p_matrix, sparse_matrix *p_matrix_result) {
    int i, j, k, n = p_matrix->n;
    sparse_row *p_row_i, *p_row_j, *p_rrow;
    sfloat sum;

    for (i = 0, p_row_i = p_matrix->rows; i < n; i++, p_row_i++) {
        if (p_row_i->nnz == 0)
            continue;

        // j < i
        for (j = 0, p_row_j = p_matrix->rows; j < i; j++, p_row_j++) {
            if (p_row_j->nnz == 0)
                continue;
            sum = sparse_mul_rows(p_row_i, p_row_j);
            if (absolute(sum) < EPSILON)
                continue;
            p_rrow = p_matrix_result->rows + i;
            p_rrow->col[p_rrow->nnz] = j;
            p_rrow->data[p_rrow->nnz++] = sum;

            p_rrow = p_matrix_result->rows + j;
            p_rrow->col[p_rrow->nnz] = i;
            p_rrow->data[p_rrow->nnz++] = sum;
        }

        // j = i
        sum = 0;
        for (k = 0; k < p_row_i->nnz; k++)
            sum += square(p_row_i->data[k]);
        if (absolute(sum) < EPSILON)
            continue;
        p_rrow = p_matrix_result->rows + i;
        p_rrow->col[p_rrow->nnz] = i;
        p_rrow->data[p_rrow->nnz++] = sum;
    }
}

static sfloat sparse_mul_row_vec(const sparse_row *p_row, const sfloat *p_vec) {
    int i, *p_col;
    sfloat *p_data, sum = 0;
    for (i = 0, p_col = p_row->col, p_data = p_row->data; i < p_row->nnz; i++)
        sum += *p_data++ * p_vec[*p_col++];
    if (absolute(sum) < EPSILON)
        sum = 0;
    return sum;
}

void sparse_mul_vec(const sparse_matrix *p_matrix, const sfloat *p_vec, sfloat *p_vec_result) {
    int i;
    sparse_row *p_row;
    for (i = 0, p_row = p_matrix->rows; i < p_matrix->n; i++)
        *p_vec_result++ = sparse_mul_row_vec(p_row++, p_vec);
}

void array_to_sparse(const sfloat *arr, int n, int m, sparse_matrix **p_mat) {
    sparse_matrix *p_matrix = init_sparse(n, m, m);
    if (!p_matrix) {
        *p_mat = NULL;
        return;
    }

    int i, j;
    sparse_row *p_row = p_matrix->rows;
    for (i = 0; i < n; i++, p_row++) {
        for (j = 0; j < m; j++, arr++) {
            if (absolute(*arr) < EPSILON)
                continue;
            p_row->col[p_row->nnz] = j;
            p_row->data[p_row->nnz++] = *arr;
        }
    }

    *p_mat = p_matrix;
}
