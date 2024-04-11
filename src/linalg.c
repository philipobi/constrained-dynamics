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


void sparse_to_array(const sparse_matrix *matrix, sfloat *array)
{
    sparse_row *p_row;
    for (int i=0; i<matrix->n; i++) {
        if (!(p_row = matrix->rows[i])) continue;
        for (int k=0; k<p_row->nnz; k++) array[i*matrix->m + p_row->col[k]] = p_row->data[k];
    }
}

void print_array(sfloat *array, int n, int m)
{
    int i,j;
    for (i=0; i<n; i++) {
        for (j=0; j<m; j++) printf("%5.0f", array[i*m + j]);
        printf("\n");
    }
}


void print_sparse(const sparse_matrix *matrix) 
{  
    sfloat *arr = malloc(matrix->n * matrix->m * sizeof(sfloat));
    if (!arr) return;
    sparse_to_array(matrix, arr);
    print_array(arr, matrix->n, matrix->m);
    free(arr);
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
    return sum;
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

void test()
{
    sparse_matrix A = {.n = 6, .m = 12};
    sparse_row rows[5];
    
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

    int data3[] = {0,1,2,3,4,5};
    rows[3].col = data3;
    rows[3].data = values;
    rows[3].nnz = 6;
    
    int data4[] = {6,7,8,9,10,11};
    rows[4].col = data4;
    rows[4].data = values;
    rows[4].nnz = 6;

    sparse_row *p_rows[] = {rows, rows+1, rows+2, NULL, rows+3, rows+4};
    A.rows = p_rows;

    print_sparse(&A);
}


int main ()
{
    test();
    return 0;
}