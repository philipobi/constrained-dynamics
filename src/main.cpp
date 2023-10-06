#include <iostream>
#include <utils.h>


void grid(double *positions, Constraint *constraints, int n, int m, double l)
{
    double *p_pos = positions;
    
    for(int i=0; i<n; i++) for(int j=0; j<m; j++) 
    {
        *p_pos++ = i*l; //x
        *p_pos++ = j*l; //y
        *p_pos++ = 0; //z
    }

    Constraint *p_constraints = constraints;
    //horizontal bonds
    for(int i=0; i<n; i++) for(int j=0; j<(m-1); j++) *p_constraints++ = Constraint(l, positions+3*(i*n+j), positions+3*(i*n+j+1));

    //vertical bonds
    for(int j=0; j<m; j++) for(int i=0; i<(n-1); i++) *p_constraints++ = Constraint(l, positions+3*(i*n+j), positions+3*((i+1)*n+j));
}


int main()
{
    int n,m,nConstraints;
    n=3;
    m=3;
    nConstraints = n*(m-1) + m*(n-1);

    double *positions = new double[3*n*m];
    Constraint *constraints = new Constraint[nConstraints];

    grid(positions,constraints,n,m,3.);
    
    Jacobian J(n,m,nConstraints,constraints,positions);
    printMatrix(nConstraints,3*n*m,J.eval());

    delete[] positions; delete[] constraints;
}