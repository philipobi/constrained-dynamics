#include <globaldef.h>
#include <linalg.h>
#include <stdio.h>
#include <stdlib.h>
#define BUFSIZE_MATRIX 10000
#define BUFSIZE_LINE 1024

static int
parse_row (FILE* p_file, sfloat *arr, char *line, int bufsize, int *p_ncol)
{
    char c, *p_line = line;
    while ((c = fgetc (p_file)) != EOF && c != '\n' && p_line - line < bufsize - 1)
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

    free(line);
    *p_n = n;
    *p_m = m;
    *p_parsed_arr = arr;
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