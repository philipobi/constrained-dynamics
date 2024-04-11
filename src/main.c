#include <stdio.h>


// sfloat mul_rows(const sparse_row *p_row_i, const sparse_row *p_row_j)
// {
//     int n, k, m=0;
//     sfloat 
//     sum = 0,
//     *data1, *data2;
//     for (n = 0; n < p_row_i->nnz; n += 3) {
//         for (; m < p_row_j->nnz; m += 3) {
//             if (p_row_i->col[n] < p_row_j->col[m]) break;
//             if (p_row_i->col[n] > p_row_j->col[m]) continue;
//             data1 = p_row_i->data + n;
//             data2 = p_row_j->data + m;
//             for(k = 0; k < 3; k++) sum += (*data1++) * (*data2++);
//         }
//     }
// }

int main()
{
    int i=2;
    int out = i++ - i;
    printf("%d\n", out);
    return 0;
}