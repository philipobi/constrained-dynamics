#include <globaldef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#define min(a,b) ((a) < (b) ? (a) : (b))
#define square(a) ((a) * (a))
#define EPSILON 0.01



typedef struct sparse_row {
    int nnz;
    int *col;
    sfloat *data;
} sparse_row;

typedef struct sparse_matrix {
    int n, m;
    sparse_row **rows;
    sparse_row *row_mem;
} sparse_matrix;

sparse_matrix *init(int n, int m)
{
    sparse_matrix *matrix = malloc(sizeof(sparse_matrix)); 
    sparse_row **p_rows = malloc(n * sizeof(sparse_row *));
    sparse_row *rows = malloc(n * sizeof(sparse_row));

    if (!matrix || !p_rows || !rows) goto error;
    
    int i;
    sparse_row **p_rows_c = p_rows;
    sparse_row *p_row = rows;
    for (i = 0; i < n; i++, p_row++, p_rows_c++) {
        p_row->col = malloc(m * sizeof(int));
        p_row->data = malloc(m * sizeof(sfloat));
        if (!p_row->col || !p_row->data) goto cleanup;
        p_row->nnz = 0;
        *p_rows_c = p_row;
    }
    
    matrix->rows = p_rows;
    matrix->row_mem = rows;
    return matrix;

    cleanup:
    p_rows_c = p_rows;
    p_row = rows;
    for (int k = 0; k < i; k++, p_row++, p_rows_c++) {
        *p_rows_c = NULL;
        free(p_row->col);
        free(p_row->data);
    }
    matrix->rows = NULL;
    
    error:
    free(p_rows);
    free(rows);
    free(matrix);
    return NULL;
}

void unref(sparse_matrix *matrix)
{
    sparse_row **p_rows = matrix->rows;
    sparse_row *p_row;
    for(int i = 0; i < matrix->n; i++) {
        if ((p_row = *p_rows) == NULL) continue;
        free(p_row->col);
        free(p_row->data);
        *p_rows++ = NULL;
    }
    free(matrix->rows[0]);
}

sfloat sparse_mul_rows(const sparse_row *p_row_i, const sparse_row *p_row_j)
{
    int n, k, m=0;
    sfloat 
    sum = 0,
    *data1, *data2;
    for (n = 0; n < p_row_i->nnz; n += 3) {
        for (; m < p_row_j->nnz; m += 3) {
            if (p_row_i->col[n] < p_row_j->col[m]) break;
            if (p_row_i->col[n] > p_row_j->col[m]) continue;
            data1 = p_row_i->data + n;
            data2 = p_row_j->data + m;
            for(k = 0; k < 3; k++) sum += (*data1++) * (*data2++);
        }
    }
}

void sparse_mul_transpose(const sparse_matrix *matrix, sparse_matrix *result_matrix)
{
    int i, j, k, p, n_blk,
    n = matrix->n;
    sparse_row 
    *p_row_i, *p_row_j,
    *p_rrow_i, *p_rrow_j,
    **rows = matrix->rows,
    **rrows = result_matrix->rows;
    sfloat sum;

    for(i = 0; i < n; i++) {
        if ((p_row_i = rows[i]) == NULL) continue; // row i empty
        // j < i
        for(j = 0; j < i; j++) {
            if ((p_row_j = rows[j]) == NULL) continue; // row j empty
            sum = sparse_mul_rows(p_row_i, p_row_j);
            if (sum < EPSILON || -sum < EPSILON) continue; // -> 0
            p_rrow_i = rrows[i];
            p_rrow_i->col[p_rrow_i->nnz] = j;
            p_rrow_i->data[p_rrow_i->nnz++] = sum;
            
            p_rrow_j = rrows[j];
            p_rrow_j->col[p_rrow_j->nnz] = i;
            p_rrow_j->data[p_rrow_j->nnz++] = sum;
            
        }
        // j = i
        sum = 0;
        for (k = 0; k < p_row_i->nnz; k++) sum += square(p_row_i->data[k]);
        if (sum < EPSILON || -sum < EPSILON) continue; // -> 0
        p_rrow_i = rrows[i];
        p_rrow_i->col[p_rrow_i->nnz] = i;
        p_rrow_i->data[p_rrow_i->nnz++] = sum;
    }
}

