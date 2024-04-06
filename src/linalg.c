#include <definitions.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef struct block {
    int i;
    int j;
    int ilength;
    int jlength;
    sfloat *data;
} block;

typedef struct bs_mat {
    int n;
    int m;
    int n_blocks;
    block *blocks;
} bs_mat;


void bs_to_mat(const bs_mat *matrix, sfloat *data)
{
    block *p_block = matrix->blocks;
    for(int n = 0; n < matrix->n_blocks; n++, p_block++) {
        for(int i=0; i < p_block->ilength; i++) {
            memcpy(
                data + (matrix->m) * (p_block->i + i) + p_block->j,
                p_block->data + i * p_block->jlength,
                p_block->jlength * sizeof(sfloat)
                );
        }
    }
}

void print_mat(const sfloat *mat, int n, int m)
{   
    for(int i=0; i<n; i++) {
        for(int j=0; j<m; j++) {
            printf("%3.0f    ", mat[i*m + j]);
        }
        printf("\n");
    }
}

void print_bs_mat(const bs_mat *matrix)
{
    sfloat *data = (sfloat *) malloc(matrix->n * matrix-> m * sizeof(sfloat));
    if (data == NULL) return;
    bs_to_mat(matrix, data);
    print_mat(data, matrix->n, matrix->m);
    free(data);
}

void test()
{
    sfloat data[] = {1, 2, 3, 4};
    block blocks[] = {
        {0,0,2,2,data},
        {2,2,2,2,data},
        {4,4,2,2,data},
        {6,6,2,2,data},
        {8,8,2,2,data}
    };

    bs_mat matrix = {10, 10, 5, blocks};

    print_bs_mat(&matrix);
}


int main()
{
    test();
    return 0;
}