#include <globaldef.h>
#include <stdlib.h>
#include <stdio.h>
#define min(a,b) ((a) < (b) ? (a) : (b))
#define square(a) ((a) * (a))
#define abs(a) ((a) > 0 ? (a) : -(a))
#define EPSILON 0.01
#define BLKSIZE 3
#define NBLK 0


typedef struct sparse_row {
    int nnz;
    int *col;
    sfloat *data;
} sparse_row;

typedef struct sparse_matrix {
    int n, m;
    sparse_row *rows;
} sparse_matrix;

void destruct_rows(sparse_row *p_row, int n)
{
    for(int i = 0; i < n; i++, p_row++) {
        free(p_row->col);
        p_row->col = NULL;
        free(p_row->data);
        p_row->data = NULL;
    }
}

void free_matrix(sparse_matrix *p_matrix)
{
    free(p_matrix->rows);
    p_matrix->rows = NULL;
    free(p_matrix);
}

void destruct_matrix(sparse_matrix *p_matrix)
{
    destruct_rows(p_matrix->rows, p_matrix->n);
    free_matrix(p_matrix);
}

sparse_matrix *init(int n, int m, int max_nnz)
{
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
        int *col = malloc(max_nnz * sizeof(int));
        sfloat *data = malloc(max_nnz * sizeof(sfloat));
        if (!col || !data) {
            destruct_rows(p_matrix->rows, i + 1);
            free_matrix(p_matrix);
            return NULL;
        }
        
        p_row->nnz = 0;
        p_row->col = col;
        p_row->data = data;
    }

    return p_matrix;
}



void sparse_to_array(const sparse_matrix *p_matrix, sfloat *arr)
{
    int i, k;
    sparse_row *p_row;
    for (i = 0, p_row = p_matrix->rows; i < p_matrix->n; i++, p_row++) {
        for (k = 0; k < p_row->nnz; k++) arr[i*p_matrix->m + p_row->col[k]] = p_row->data[k];
    }
}

void print_array(sfloat *arr, int n, int m)
{
    int i,j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) printf("%3.0f, ", arr[i*m + j]);
        printf("\n");
    }
}


void print_sparse(const sparse_matrix *p_matrix) 
{  
    sfloat *arr = malloc(p_matrix->n * p_matrix->m * sizeof(sfloat));
    if (!arr) return;
    sparse_to_array(p_matrix, arr);
    print_array(arr, p_matrix->n, p_matrix->m);
    free(arr);
}

sfloat sparse_mul_rows(const sparse_row *p_row_i, const sparse_row *p_row_j)
{
    int n, k, m = 0;
    sfloat 
    sum = 0,
    *data_i, *data_j;
    
    for (n = 0; n < p_row_i->nnz; n += BLKSIZE) {
        for (; m < p_row_j->nnz; m += BLKSIZE) {
            if (p_row_i->col[n] < p_row_j->col[m]) break;
            if (p_row_i->col[n] > p_row_j->col[m]) continue;
            data_i = p_row_i->data + n;
            data_j = p_row_j->data + m;
            for(k = 0; k < BLKSIZE; k++) sum += (*data_i++) * (*data_j++);
        }
    }
    
    return sum;
}

void sparse_mul_transpose(const sparse_matrix *p_matrix, sparse_matrix *p_result_matrix)
{
    int i, j, k,
    n = p_matrix->n;
    sparse_row 
    *p_row_i, *p_row_j,
    *p_rrow;
    sfloat sum;

    for(i = 0, p_row_i = p_matrix->rows; i < n; i++, p_row_i++) {
        if (p_row_i->nnz == 0) continue;
        
        // j < i
        for(j = 0, p_row_j = p_matrix->rows; j < i; j++, p_row_j++) {
            if (p_row_j->nnz == 0) continue;
            sum = sparse_mul_rows(p_row_i, p_row_j);
            if (abs(sum) < EPSILON) continue;
            p_rrow = p_result_matrix->rows + i;
            p_rrow->col[p_rrow->nnz] = j;
            p_rrow->data[p_rrow->nnz++] = sum;
            
            p_rrow = p_result_matrix->rows + j;
            p_rrow->col[p_rrow->nnz] = i;
            p_rrow->data[p_rrow->nnz++] = sum;
            
        }
        
        // j = i
        sum = 0;
        for (k = 0; k < p_row_i->nnz; k++) sum += square(p_row_i->data[k]);
        if (abs(sum) < EPSILON) continue;
        p_rrow = p_result_matrix->rows + i;
        p_rrow->col[p_rrow->nnz] = i;
        p_rrow->data[p_rrow->nnz++] = sum;
    }
}

void test()
{
    sparse_matrix A = {.n = 6, .m = 12};
    sparse_row rows[6];
    
    sfloat values[] = {1,1,1,1,1,1};

    int data0[] = {0,1,2,6,7,8};
    rows[0].col = data0;
    rows[0].data = values;
    rows[0].nnz = 6;
    
    int data1[] = {3,4,5,9,10,11};
    rows[1].col = data1;
    rows[1].data = values;
    rows[1].nnz = 6;
    
    int data2[] = {0,1,2,9,10,11};
    rows[2].col = data2;
    rows[2].data = values;
    rows[2].nnz = 6;

    rows[3].nnz = 0;
    
    int data4[] = {0,1,2,3,4,5};
    rows[4].col = data4;
    rows[4].data = values;
    rows[4].nnz = 6;
    
    int data5[] = {6,7,8,9,10,11};
    rows[5].col = data5;
    rows[5].data = values;
    rows[5].nnz = 6;

    A.rows = rows;

    int 
    n = 6,
    m = 6;
    sparse_matrix *p_result = init(n, m, m);
    sparse_mul_transpose(&A, p_result);
    print_sparse(p_result);
    destruct_matrix(p_result);
}


int main ()
{
    test();
    return 0;
}