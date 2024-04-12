#include <globaldef.h>
#include <stdio.h>
#include <stdlib.h>
#define min(a, b) ((a) < (b) ? (a) : (b))
#define square(a) ((a) * (a))
#define abs(a) ((a) > 0 ? (a) : -(a))
#define EPSILON 0.01
#define BLKSIZE 3
#define NBLK 0
#define BUFSIZE_MATRIX 10000
#define BUFSIZE_NUMBER 100

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

void
destruct_rows (sparse_row *p_row, int n)
{
    for (int i = 0; i < n; i++, p_row++)
        {
            free (p_row->col);
            p_row->col = NULL;
            free (p_row->data);
            p_row->data = NULL;
        }
}

void
free_matrix (sparse_matrix *p_matrix)
{
    free (p_matrix->rows);
    p_matrix->rows = NULL;
    free (p_matrix);
}

void
destruct_matrix (sparse_matrix *p_matrix)
{
    destruct_rows (p_matrix->rows, p_matrix->n);
    free_matrix (p_matrix);
}

sparse_matrix *
init (int n, int m, int max_nnz)
{
    sparse_matrix *p_matrix = malloc (sizeof (sparse_matrix));
    sparse_row *rows = malloc (n * sizeof (sparse_row));
    if (!p_matrix || !rows)
        {
            free (p_matrix);
            free (rows);
            return NULL;
        }

    p_matrix->n = n;
    p_matrix->m = m;
    p_matrix->rows = rows;

    int i;
    sparse_row *p_row;
    for (i = 0, p_row = p_matrix->rows; i < n; i++, p_row++)
        {
            int *col = malloc (max_nnz * sizeof (int));
            sfloat *data = malloc (max_nnz * sizeof (sfloat));
            if (!col || !data)
                {
                    destruct_rows (p_matrix->rows, i + 1);
                    free_matrix (p_matrix);
                    return NULL;
                }

            p_row->nnz = 0;
            p_row->col = col;
            p_row->data = data;
        }

    return p_matrix;
}

void
reset (sparse_matrix *p_matrix)
{
    int i;
    sparse_row *p_row;
    for (i = 0, p_row = p_matrix->rows; i < p_matrix->n; i++, p_row++)
        p_row->nnz = 0;
}

void
sparse_to_array (const sparse_matrix *p_matrix, sfloat *arr)
{
    int i, k;
    sparse_row *p_row;
    for (i = 0, p_row = p_matrix->rows; i < p_matrix->n; i++, p_row++)
        {
            for (k = 0; k < p_row->nnz; k++)
                arr[i * p_matrix->m + p_row->col[k]] = p_row->data[k];
        }
}

void
print_array (sfloat *arr, int n, int m)
{
    int i, j;
    for (i = 0; i < n; i++)
        {
            for (j = 0; j < m; j++)
                printf ("%3.0f, ", arr[i * m + j]);
            printf ("\n");
        }
}

void
print_sparse (const sparse_matrix *p_matrix)
{
    sfloat *arr = malloc (p_matrix->n * p_matrix->m * sizeof (sfloat));
    if (!arr)
        return;
    sparse_to_array (p_matrix, arr);
    print_array (arr, p_matrix->n, p_matrix->m);
    free (arr);
}

sfloat
sparse_mul_rows (const sparse_row *p_row_i, const sparse_row *p_row_j)
{
    int n, k, m = 0;
    sfloat *data_i, *data_j, sum = 0;

    for (n = 0; n < p_row_i->nnz; n += BLKSIZE)
        {
            for (; m < p_row_j->nnz; m += BLKSIZE)
                {
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

void
sparse_mul_transpose (const sparse_matrix *p_matrix,
                      sparse_matrix *p_matrix_result)
{
    int i, j, k, n = p_matrix->n;
    sparse_row *p_row_i, *p_row_j, *p_rrow;
    sfloat sum;

    for (i = 0, p_row_i = p_matrix->rows; i < n; i++, p_row_i++)
        {
            if (p_row_i->nnz == 0)
                continue;

            // j < i
            for (j = 0, p_row_j = p_matrix->rows; j < i; j++, p_row_j++)
                {
                    if (p_row_j->nnz == 0)
                        continue;
                    sum = sparse_mul_rows (p_row_i, p_row_j);
                    if (abs (sum) < EPSILON)
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
                sum += square (p_row_i->data[k]);
            if (abs (sum) < EPSILON)
                continue;
            p_rrow = p_matrix_result->rows + i;
            p_rrow->col[p_rrow->nnz] = i;
            p_rrow->data[p_rrow->nnz++] = sum;
        }
}

sfloat
sparse_mul_row_vec (const sparse_row *p_row, const sfloat *p_vec)
{
    int i, *p_col;
    sfloat *p_data, sum = 0;
    for (i = 0, p_col = p_row->col, p_data = p_row->data; i < p_row->nnz; i++)
        sum += *p_data++ * p_vec[*p_col++];
    if (abs (sum) < EPSILON)
        sum = 0;
    return sum;
}

void
sparse_mul_vec (const sparse_matrix *p_matrix, const sfloat *p_vec,
                sfloat *p_vec_result)
{
    int i;
    sparse_row *p_row;
    for (i = 0, p_row = p_matrix->rows; i < p_matrix->n; i++)
        *p_vec_result++ = sparse_mul_row_vec (p_row++, p_vec);
}

void
test ()
{
    sparse_matrix A = { .n = 6, .m = 12 };
    sparse_row rows[6];

    sfloat values[] = { 1, 1, 1, 1, 1, 1 };

    int data0[] = { 0, 1, 2, 6, 7, 8 };
    rows[0].col = data0;
    rows[0].data = values;
    rows[0].nnz = 6;

    int data1[] = { 3, 4, 5, 9, 10, 11 };
    rows[1].col = data1;
    rows[1].data = values;
    rows[1].nnz = 6;

    int data2[] = { 0, 1, 2, 9, 10, 11 };
    rows[2].col = data2;
    rows[2].data = values;
    rows[2].nnz = 6;

    rows[3].nnz = 0;

    int data4[] = { 0, 1, 2, 3, 4, 5 };
    rows[4].col = data4;
    rows[4].data = values;
    rows[4].nnz = 6;

    int data5[] = { 6, 7, 8, 9, 10, 11 };
    rows[5].col = data5;
    rows[5].data = values;
    rows[5].nnz = 6;

    A.rows = rows;

    int n = 6, m = 6;
    sparse_matrix *p_result = init (n, m, m);
    sparse_mul_transpose (&A, p_result);
    print_sparse (p_result);
    destruct_matrix (p_result);
}

int 
parse_row (sfloat *arr, char *str, int bufsize_number, int * p_ncol)
{
    int in_row = 0, in_num = 0, ncol = 0;
    char c, *p_str = str;
    while ((c = getchar()) != EOF && c != '\n' && p_str - str < bufsize_number) {
        if (!in_num && c != ' ') {
            in_num = 1;
            *p_str++ = c;
            continue;
        }
        if (in_num && c == ' ') {
            in_num = 0;
            *p_str++ = '\0';
            *arr++ = strtof(str, NULL);
            ncol++;
            p_str = str;
            continue;
        }
        *p_str++ = c;
    }
    
    if (in_num && p_str - str < bufsize_number) {
        *p_str++ = '\0';
        *arr++ = strtof(str, NULL);
        ncol++;
    }

    *p_ncol = ncol;
    
    return c == EOF ? 1 : 0;
}

sfloat * 
parse_array (int *p_n, int *p_m)
{
    sfloat *arr = malloc(BUFSIZE_MATRIX * sizeof(sfloat));
    char *str = malloc(BUFSIZE_NUMBER * sizeof(char));
    if (!arr || !str) {
        free(arr);
        free(str);
        *p_n = 0;
        *p_m = 0;
        return NULL;
    };

    sfloat *p_arr = arr;
    int n = 0, m = 0, ncol, eof_flag;
    do {
        eof_flag = parse_row(p_arr, str, BUFSIZE_NUMBER, &ncol);
        if(ncol) {
            n++;
            p_arr += ncol;
            m = ncol;
            }

    } while (!eof_flag);

    *p_n = n;
    *p_m = m;
    return arr;
}

void
test_parse (void)
{
    int n, m;
    sfloat *arr = parse_array(&n, &m);
    print_array(arr, n, m);
    free (arr);
}

int
main ()
{
    test_parse ();
    return 0;
}