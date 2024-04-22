#include <globaldef.h>
#include <linalg.h>
#include <stdio.h>
#include <stdlib.h>
#define BUFSIZE_MATRIX 10000
#define BUFSIZE_LINE 1024

static int
parse_row (FILE *p_file, sfloat *arr, char *line, int bufsize, int *p_ncol)
{
    char c, *p_line = line;
    while ((c = fgetc (p_file)) != EOF && c != '\n'
           && p_line - line < bufsize - 1)
        *p_line++ = c;
    *p_line++ = '\0';

    sfloat val, *p_arr = arr;
    char *p1 = line, *p2;
    while ((val = strtof (p1, &p2)) || p2 - p1)
        {
            p1 = p2;
            *p_arr++ = val;
        }
    *p_ncol = p_arr - arr;

    return c == EOF ? 1 : 0;
}

void
parse_array (FILE *p_file, sfloat **p_parsed_arr, int *p_n, int *p_m)
{
    sfloat *arr = malloc (BUFSIZE_MATRIX * sizeof (sfloat));
    char *line = malloc (BUFSIZE_LINE * sizeof (char));
    if (!arr || !line)
        {
            free (arr);
            free (line);
            *p_n = 0;
            *p_m = 0;
            *p_parsed_arr = NULL;
            return;
        };

    sfloat *p_arr = arr;
    int n = 0, m = 0, ncol, eof_flag;
    do
        {
            eof_flag = parse_row (p_file, p_arr, line, BUFSIZE_LINE, &ncol);
            if (ncol)
                {
                    n++;
                    p_arr += ncol;
                    m = ncol;
                }
        }
    while (!eof_flag);

    free (line);
    *p_n = n;
    *p_m = m;
    *p_parsed_arr = arr;
}

void
print_array (const sfloat *p_arr, const int n, const int m, const char *str)
{
    printf ("(%d, %d)\t%s\n", n, m, str);
    int i, j;
    for (i = 0; i < n; i++)
        {
            for (j = 0; j < m; j++)
                printf ("%3.0f, ", p_arr[i * m + j]);
            printf ("\n");
        }
}

sfloat
dot (const sfloat *p_vec1, const sfloat *p_vec2, const int n)
{
    int i;
    sfloat sum;
    for (i = 0, sum = 0; i < n; i++)
        sum += *p_vec1++ * *p_vec2++;
    return sum;
}

sfloat
norm2 (const sfloat *p_vec, const int n)
{
    int i;
    sfloat sum;
    for (i = 0, sum = 0; i < n; i++, p_vec++)
        sum += square (*p_vec);
    return sum;
}

void 
add_vec (const sfloat *p_x, const sfloat a, const sfloat *p_v, sfloat *p_y, const int n)
// performs y = x + a*v
{
    for (int i = 0; i < n; i++) *p_y++ = *p_x++ + a * *p_v++;
}

void
add_vec_inplace (sfloat *p_x, const sfloat a, const sfloat *p_v, const int n)
// performs x = x + a*v
{
    for (int i = 0; i < n; i++)
        *p_x++ += a * *p_v++;
}

void
add_2vec_inplace (sfloat *p_x, const sfloat a, const sfloat *p_u, const sfloat b,
          const sfloat *p_v, const int n)
// performs x = x + a*u + b*v
{
    for (int i = 0; i < n; i++, p_x++)
        *p_x = *p_x + a * *p_u++ + b * *p_v++;
}

void
copy_vec (const sfloat *p_src, sfloat *p_dest, const int n)
{
    for (int i = 0; i < n; i++) *p_dest++ = *p_src++;
}