#include <globaldef.h>
#include <linalg.h>
#include <sparse_linalg.h>
#include <stdio.h>
#include <stdlib.h>

static void destruct_rows(sparse_row *rows) {
    if (!rows)
        return;
    free(rows[0].col);
    free(rows[0].data);
    free(rows);
}

void destruct_sparse(sparse_matrix *p_matrix) {
    if (!p_matrix)
        return;
    destruct_rows(p_matrix->rows);
    free(p_matrix);
}

sparse_row *init_rows(int n, int max_nnz) {
    int total_nnz = n * max_nnz;
    sparse_row *rows = NULL;
    int *cols = NULL;
    sfloat *data = NULL;

    if (                                                            //
        !(rows = calloc(n, sizeof(sparse_row))) ||                  //
        !(cols = rows[0].col = malloc(total_nnz * sizeof(int))) ||  //
        !(data = rows[0].data = malloc(total_nnz * sizeof(sfloat))) //
    ) {
        destruct_rows(rows);
        return NULL;
    }

    for (int i = 0; i < n; i++, cols += max_nnz, data += max_nnz) {
        rows[i].col = cols;
        rows[i].data = data;
    }

    return rows;
}

sparse_matrix *init_sparse(int n, int m, int max_nnz) {
    sparse_matrix *p_matrix = NULL;
    sparse_row *rows = malloc(n * sizeof(sparse_row));
    if (                                                  //
        !(p_matrix = calloc(1, sizeof(sparse_matrix))) || //
        !(rows = p_matrix->rows = init_rows(n, max_nnz))  //
    ) {
        destruct_sparse(p_matrix);
        return NULL;
    }

    p_matrix->n = n;
    p_matrix->m = m;
    return p_matrix;
}

void reset_sparse(sparse_matrix *p_matrix) {
    for (int i = 0; i < p_matrix->n; i++)
        p_matrix->rows[i].nnz = 0;
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
    print_array(arr, p_matrix->n, p_matrix->m);
    free(arr);
}

static sfloat sparse_mul_rows(const sparse_row *p_row_i, const sparse_row *p_row_j, const int blksize) {
    int n, k, m = 0;
    sfloat *data_i, *data_j, sum = 0;

    for (n = 0; n < p_row_i->nnz; n += blksize) {
        for (; m < p_row_j->nnz; m += blksize) {
            if (p_row_i->col[n] < p_row_j->col[m])
                break;
            if (p_row_i->col[n] > p_row_j->col[m])
                continue;
            data_i = p_row_i->data + n;
            data_j = p_row_j->data + m;
            for (k = 0; k < blksize; k++)
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
        // j < i
        for (j = 0, p_row_j = p_matrix->rows; j < i; j++, p_row_j++) {
            if ((sum = sparse_mul_rows(p_row_i, p_row_j, BLKSIZE)) == 0)
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
        if (sum == 0)
            continue;
        p_rrow = p_matrix_result->rows + i;
        p_rrow->col[p_rrow->nnz] = i;
        p_rrow->data[p_rrow->nnz++] = sum;
    }
}

static sfloat sparse_mul_row_vec(const sparse_row *p_row, const sfloat *p_vec) {
    int i;
    int *p_col = p_row->col;
    sfloat *p_data, sum = 0;
    for (i = 0; i < p_row->nnz; i++)
        sum += p_row->data[i] * p_vec[p_row->col[i]];
    return sum;
}

void sparse_mul_vec(const sparse_matrix *p_matrix, const sfloat *p_vec, sfloat *p_vec_result) {
    sparse_row *p_row = p_matrix->rows;
    for (int i = 0; i < p_matrix->n; i++)
        *p_vec_result++ = sparse_mul_row_vec(p_row++, p_vec);
}

void sparse_transpose_mul_vec(const sparse_matrix *p_matrix, const sfloat *p_vec, sfloat *p_vec_result) {
    int i, j;
    const sparse_row *p_row = p_matrix->rows;
    const int *p_col;
    const sfloat *p_data;
    for (i = 0; i < p_matrix->n; i++, p_row++, p_vec++) {
        p_col = p_row->col;
        p_data = p_row->data;
        for (j = 0; j < p_row->nnz; j++) {
            p_vec_result[*p_col++] += *p_vec * *p_data++;
        }
    }
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
