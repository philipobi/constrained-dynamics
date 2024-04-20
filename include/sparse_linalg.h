#pragma once
#include <globaldef.h>

#define BLKSIZE 3
#define NBLK 0

typedef struct sparse_row
{
    int nnz;
    int *col;
    sfloat *data;
} sparse_row;

typedef struct sparse_matrix
{
    int n, m;
    sparse_row *rows;
} sparse_matrix;

sparse_matrix *init (int n, int m, int max_nnz);

void destruct_matrix (sparse_matrix *p_matrix);

void sparse_mul_transpose (const sparse_matrix *p_matrix,
                           sparse_matrix *p_matrix_result);

void sparse_mul_vec (const sparse_matrix *p_matrix, const sfloat *p_vec,
                     sfloat *p_vec_result);

void array_to_sparse (const sfloat *arr, int n, int m, sparse_matrix **p_mat);