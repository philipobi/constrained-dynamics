#include <globaldef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef struct coo_mat {
    int n, m, nnz;
    int *row_coords, *col_coords;
    sfloat *data;
} coo_mat;

void coo_to_arr(const coo_mat* matrix, sfloat *arr)
{
    int 
    *p_row_coord = matrix->row_coords,
    *p_col_coord = matrix->col_coords,
    m = matrix->m;
    sfloat *p_data = matrix->data;

    for(int k = 0; k < matrix->nnz; k++, p_row_coord++, p_col_coord++, p_data++) {
        arr[*p_row_coord * m + *p_col_coord] = *p_data;
    }
}

void print_arr(const sfloat *arr, int n, int m)
{
    int i,j;
    for(i = 0; i < n; i++) {
        for(j = 0; j < m; j++) {
            printf("%2.0f    ", arr[i*m + j]);
        }
        printf("\n");
    }
}

void print_coo(const coo_mat *matrix)
{
    int
    n = matrix->n,
    m = matrix->m;
    sfloat *arr = (sfloat *) malloc(n * m * sizeof(sfloat));
    if(arr == NULL) return;

    coo_to_arr(matrix, arr);
    print_arr(arr, n, m);

    free(arr);
}

void test()
{
    int row_coords[] = {0,0,1,1,2,2,3,3};
    int col_coords[] = {0,1,0,1,2,3,2,3};
    sfloat data[] = {1,2,3,4,5,6,7,8};
    coo_mat matrix = {
        .n = 4,
        .m = 4,
        .nnz = 8,
        .row_coords = row_coords,
        .col_coords = col_coords,
        .data = data
    };
    print_coo(&matrix);
}


int main()
{
    test();
    return 0;
}