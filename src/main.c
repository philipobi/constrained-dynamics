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
    int row1[] = {
        0,1,2,
        6,7,8,
        9,10,11
    };
    int row2[] = {
        3,4,5,
        6,7,8,
        12,13,14
    };

    int 
    m=0,
    *data1, *data2;

    for (int n=0; n < 9; n+=3) {
        for (; m<9; m+=3) {
            if(row1[n] < row2[m]) break;
            if(row1[n] > row2[m]) continue;
            data1 = row1 + n;
            data2 = row2 + m;
            for (int k=0; k<3; k++) printf("(%i,%i)\n",*data1++,*data2++);
        }
    }
    return 0;
}